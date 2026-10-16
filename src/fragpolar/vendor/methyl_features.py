#!/usr/bin/env python3
"""Extract read-level methylation features from BAM or MethylDackel bedGraph.

The BAM path can compute AMF, PDR, PM, MHL, FDRP, and qFDRP because read-level
methylation haplotypes are available. A MethylDackel bedGraph is aggregate
site-level output, so it supports exact AMF only when methylated/unmethylated
counts are present.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import pysam
except ImportError:  # pragma: no cover - handled by CLI at runtime
    pysam = None


@dataclass
class Region:
    chrom: str
    start: int
    end: int


@dataclass
class ReadMethylCalls:
    read_id: str
    chrom: str
    start: int
    end: int
    calls: Dict[int, int]
    conversion_context: str = "unknown"
    alignment_count: int = 1
    is_paired_fragment: bool = False
    is_complete_pair: bool = False
    mate_overlap_cpgs: int = 0
    mate_conflict_cpgs: int = 0


@dataclass
class FeatureResult:
    summary: Dict[str, float]
    counts: Dict[str, int]
    pm_windows: List[Dict[str, object]]
    mhl_by_length: List[Dict[str, object]]
    pairwise: List[Dict[str, object]]


def parse_region(text: str) -> Region:
    chrom, rest = text.replace(",", "").split(":", 1)
    start_text, end_text = rest.split("-", 1)
    return Region(chrom=chrom, start=int(start_text), end=int(end_text))


def load_cpg_sites(path: str, region: Region) -> List[int]:
    """Load 0-based CpG cytosine positions from BED-like or TSV input."""
    positions: List[int] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom = fields[0]
            if chrom != region.chrom:
                continue
            if len(fields) >= 3:
                pos = int(fields[1])
            else:
                pos = int(fields[1])
            if region.start <= pos < region.end:
                positions.append(pos)
    return sorted(set(positions))


def parse_methyl_tag(tag_value: str, region: Region) -> Dict[int, int]:
    """Parse a simple methyl-call tag: ``100:M,110:U`` or ``100:1,110:0``."""
    calls: Dict[int, int] = {}
    if not tag_value:
        return calls
    for item in tag_value.replace(";", ",").split(","):
        item = item.strip()
        if not item:
            continue
        pos_text, state_text = item.split(":", 1)
        pos = int(pos_text)
        state = state_text.strip().upper()
        if not (region.start <= pos < region.end):
            continue
        if state in {"M", "1", "C", "Z"}:
            calls[pos] = 1
        elif state in {"U", "0", "T", "X"}:
            calls[pos] = 0
        else:
            raise ValueError("Unsupported methylation state in tag: %s" % item)
    return calls


def _call_from_base(base: str, conversion_context: str) -> Optional[int]:
    base = base.upper()
    if conversion_context == "top_ct":
        if base == "C":
            return 1
        if base == "T":
            return 0
    elif conversion_context == "bottom_ga":
        if base == "G":
            return 1
        if base == "A":
            return 0
    return None


def _alignment_conversion_hint(aln) -> str:
    """Infer conversion mode from paired-end orientation.

    - 99/147-style fragments use CG -> TG, so call C/T at the CpG C coordinate.
    - 163/83-style fragments use CG -> CA, so call G/A at the CpG G coordinate.
    """
    if aln.is_paired:
        if (aln.is_read1 and not aln.is_reverse) or (aln.is_read2 and aln.is_reverse):
            return "top_ct"
        if (aln.is_read1 and aln.is_reverse) or (aln.is_read2 and not aln.is_reverse):
            return "bottom_ga"
    return "top_ct" if not aln.is_reverse else "bottom_ga"


def _fragment_conversion_context(alns: Sequence[object]) -> str:
    hints = Counter(_alignment_conversion_hint(aln) for aln in alns)
    if not hints:
        return "unknown"
    if len(hints) == 1:
        return next(iter(hints))
    return hints.most_common(1)[0][0]


def _calls_from_alignment_by_context(
    aln,
    cpg_positions: Sequence[int],
    conversion_context: str,
    min_base_quality: int = 0,
    trim_5prime: int = 0,
    trim_3prime: int = 0,
) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Call methylation from aligned bases with optional quality and M-bias filters.

    Returns (calls, call_qualities) where call_qualities maps CpG position to
    the base quality at that position, used for quality-weighted overlap
    resolution.

    ``trim_5prime``/``trim_3prime`` exclude positions near the sequencing 5'/3'
    ends to mitigate M-bias artifacts common in EM-seq and bisulfite libraries.
    The trim is applied in *sequencing* coordinates: for reverse-mapped reads the
    BAM query position is converted back to the original sequencing cycle before
    checking.
    """
    calls: Dict[int, int] = {}
    call_quals: Dict[int, int] = {}
    query = aln.query_sequence or ""
    qualities = aln.query_qualities
    read_len = len(query)
    if read_len == 0:
        return calls, call_quals
    ref_to_query = {
        ref_pos: query_pos
        for query_pos, ref_pos in aln.get_aligned_pairs(matches_only=True)
        if query_pos is not None and ref_pos is not None
    }
    for cpg_c_pos in cpg_positions:
        ref_pos = cpg_c_pos if conversion_context == "top_ct" else cpg_c_pos + 1
        query_pos = ref_to_query.get(ref_pos)
        if query_pos is None or query_pos >= read_len:
            continue
        bq = qualities[query_pos] if qualities is not None else 0
        if min_base_quality > 0 and bq < min_base_quality:
            continue
        if trim_5prime > 0 or trim_3prime > 0:
            seq_pos = (read_len - 1 - query_pos) if aln.is_reverse else query_pos
            if seq_pos < trim_5prime or seq_pos >= read_len - trim_3prime:
                continue
        call = _call_from_base(query[query_pos], conversion_context)
        if call is not None:
            calls[cpg_c_pos] = call
            call_quals[cpg_c_pos] = bq
    return calls, call_quals


def read_calls_from_bam(
    bam_path: str,
    region: Region,
    cpg_positions: Optional[Sequence[int]] = None,
    call_source: str = "tag",
    methyl_tag: str = "XM",
    merge_paired_end: bool = True,
    overlap_policy: str = "quality",
    skip_duplicates: bool = True,
    min_mapping_quality: int = 10,
    min_base_quality: int = 5,
    trim_5prime: int = 0,
    trim_3prime: int = 0,
    trim_5prime_r1: Optional[int] = None,
    trim_3prime_r1: Optional[int] = None,
    trim_5prime_r2: Optional[int] = None,
    trim_3prime_r2: Optional[int] = None,
    mbias_trim: Optional[Dict[Tuple[str, str], Tuple[int, int]]] = None,
    return_qc: bool = False,
):
    """Extract per-fragment methylation calls from BAM.

    ``call_source=tag`` reads a simple synthetic/fixture tag such as
    ``XM=100:M,110:U``. ``call_source=pe_strand`` calls methylation from read
    bases using paired-end conversion orientation: 99/147-style fragments use
    C/T at the CpG C coordinate, and 163/83-style fragments use G/A at the CpG G
    coordinate. ``call_source=c_to_t`` is kept as a backwards-compatible alias
    for this PE-aware base mode.

    Assumes a directional EM-seq (or directional bisulfite) library.  Non-
    directional and PBAT libraries use different strand conventions and are not
    currently supported.

    For paired-end BAMs, alignments with the same query name are merged into one
    fragment haplotype before feature calculation. Overlapping mate CpGs are
    counted once; conflicting mate calls are resolved by the chosen
    ``overlap_policy``.

    Quality filters (``min_mapping_quality``, ``min_base_quality``,
    ``trim_5prime``, ``trim_3prime``) are applied during extraction.  MAPQ
    filtering applies to all call sources; base quality and M-bias trim only
    apply to base-level modes (``pe_strand``/``c_to_t``).

    M-bias trim can be specified per read in a pair via ``trim_5prime_r1``,
    ``trim_3prime_r1``, ``trim_5prime_r2``, ``trim_3prime_r2``.  When a
    per-read value is ``None`` (the default), the corresponding global
    ``trim_5prime``/``trim_3prime`` is used as a fallback.

    For strand-specific M-bias trim, pass ``mbias_trim`` — a dict mapping
    ``(conversion_context, read_key)`` to ``(trim_5prime, trim_3prime)``.
    Example for EM-seq directional libraries::

        mbias_trim = {
            ("top_ct", "R1"):     (3, 0),
            ("top_ct", "R2"):     (0, 20),
            ("bottom_ga", "R1"):  (5, 3),
            ("bottom_ga", "R2"):  (20, 0),
        }

    When ``mbias_trim`` is provided and contains the matching key, it takes
    precedence over per-read and global trim parameters.  ``mbias_trim``
    only affects CpG methylation call extraction — it does NOT discard
    reads, so non-methylation features (coverage, fragment length, etc.)
    are unaffected.

    Overlap policies:
      - ``quality``: MethylDackel-style quality-weighted resolution.  When
        mates agree the higher-quality call is kept.  When they disagree the
        higher-quality call wins with its quality reduced by the opposing
        quality; equal-quality conflicts are dropped.
      - ``drop_conflict``: any mate disagreement drops the CpG from the
        fragment.
      - ``first``: keep the first mate's call on conflict.
    """
    if pysam is None:
        raise RuntimeError("pysam is required for BAM extraction")
    if call_source not in {"tag", "c_to_t", "pe_strand"}:
        raise ValueError("call_source must be 'tag', 'c_to_t', or 'pe_strand'")
    if overlap_policy not in {"quality", "drop_conflict", "first"}:
        raise ValueError("overlap_policy must be 'quality', 'drop_conflict', or 'first'")
    cpg_set = set(cpg_positions or [])
    if call_source in {"c_to_t", "pe_strand"} and not cpg_set:
        raise ValueError("cpg_positions are required for base-based call_source")

    t5_r1 = trim_5prime_r1 if trim_5prime_r1 is not None else trim_5prime
    t3_r1 = trim_3prime_r1 if trim_3prime_r1 is not None else trim_3prime
    t5_r2 = trim_5prime_r2 if trim_5prime_r2 is not None else trim_5prime
    t3_r2 = trim_3prime_r2 if trim_3prime_r2 is not None else trim_3prime

    # Records carry (alignment, calls, per-position-qualities).
    # For tag mode qualities are empty since BQ is not available.
    grouped: Dict[str, List[Tuple[object, Dict[int, int], Dict[int, int]]]] = defaultdict(list)
    qc: Dict[str, int] = defaultdict(int)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(region.chrom, region.start, region.end):
            qc["n_alignments_seen"] += 1
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                qc["n_alignments_skipped_nonprimary"] += 1
                continue
            if aln.is_qcfail:
                qc["n_alignments_skipped_qcfail"] += 1
                continue
            if skip_duplicates and aln.is_duplicate:
                qc["n_duplicate_alignments_skipped"] += 1
                continue
            if aln.mapping_quality < min_mapping_quality:
                qc["n_alignments_skipped_low_mapq"] += 1
                continue
            if call_source == "tag":
                if not aln.has_tag(methyl_tag):
                    qc["n_alignments_missing_methyl_tag"] += 1
                    continue
                calls = parse_methyl_tag(aln.get_tag(methyl_tag), region)
            else:
                calls = {}
            if call_source != "tag" or calls:
                if calls:
                    qc["n_alignments_with_calls"] += 1
                grouped[aln.query_name].append((aln, calls, {}))

    if call_source in {"c_to_t", "pe_strand"}:
        grouped_with_calls: Dict[str, List[Tuple[object, Dict[int, int], Dict[int, int]]]] = defaultdict(list)
        for query_name, records in grouped.items():
            alns = [r[0] for r in records]
            conversion_context = _fragment_conversion_context(alns)
            hint_counts = Counter(_alignment_conversion_hint(aln) for aln in alns)
            if len(hint_counts) > 1:
                qc["n_mixed_conversion_hint_fragments"] += 1
            for aln, _, _ in records:
                read_key = "R2" if aln.is_read2 else "R1"
                if mbias_trim and (conversion_context, read_key) in mbias_trim:
                    t5, t3 = mbias_trim[(conversion_context, read_key)]
                elif aln.is_read2:
                    t5, t3 = t5_r2, t3_r2
                else:
                    t5, t3 = t5_r1, t3_r1
                calls, quals = _calls_from_alignment_by_context(
                    aln,
                    cpg_positions or [],
                    conversion_context,
                    min_base_quality=min_base_quality,
                    trim_5prime=t5,
                    trim_3prime=t3,
                )
                if calls:
                    qc["n_alignments_with_calls"] += 1
                    grouped_with_calls[query_name].append((aln, calls, quals))
        grouped = grouped_with_calls
    else:
        grouped = {
            query_name: records
            for query_name, records in grouped.items()
            if any(calls for _, calls, _ in records)
        }
        qc["n_alignments_with_calls"] = sum(len(records) for records in grouped.values())

    fragments: List[ReadMethylCalls] = []
    for query_name, records in grouped.items():
        qc["n_query_names_with_calls"] += 1
        if not merge_paired_end:
            for idx, (aln, calls, _quals) in enumerate(records):
                fragments.append(
                    ReadMethylCalls(
                        read_id="%s/%s" % (query_name, idx + 1),
                        chrom=region.chrom,
                        start=aln.reference_start,
                        end=aln.reference_end,
                        calls=dict(calls),
                        conversion_context=_fragment_conversion_context([aln]),
                        alignment_count=1,
                        is_paired_fragment=bool(aln.is_paired),
                        is_complete_pair=False,
                    )
                )
            continue

        is_paired = any(r[0].is_paired for r in records)
        has_read1 = any(r[0].is_read1 for r in records)
        has_read2 = any(r[0].is_read2 for r in records)
        is_complete_pair = is_paired and has_read1 and has_read2
        conversion_context = _fragment_conversion_context([r[0] for r in records])
        if conversion_context == "top_ct":
            qc["n_top_ct_fragments"] += 1
        elif conversion_context == "bottom_ga":
            qc["n_bottom_ga_fragments"] += 1
        if is_complete_pair:
            qc["n_complete_pairs_with_calls"] += 1
        elif is_paired:
            qc["n_singletons_with_calls"] += 1
        else:
            qc["n_unpaired_reads_with_calls"] += 1

        merged: Dict[int, int] = {}
        merged_qual: Dict[int, int] = {}
        conflicts = set()
        mate_overlap_cpgs = 0
        mate_conflict_cpgs = 0
        for _aln, calls, quals in records:
            for pos, state in calls.items():
                q = quals.get(pos, 0)
                if pos not in merged:
                    merged[pos] = state
                    merged_qual[pos] = q
                    continue
                mate_overlap_cpgs += 1
                if merged[pos] == state:
                    if overlap_policy == "quality":
                        merged_qual[pos] = max(merged_qual.get(pos, 0), q)
                    continue
                mate_conflict_cpgs += 1
                if overlap_policy == "quality":
                    existing_q = merged_qual.get(pos, 0)
                    if q > existing_q:
                        merged[pos] = state
                        merged_qual[pos] = q - existing_q
                    elif existing_q > q:
                        merged_qual[pos] = existing_q - q
                    else:
                        conflicts.add(pos)
                elif overlap_policy == "drop_conflict":
                    conflicts.add(pos)
        for pos in conflicts:
            merged.pop(pos, None)
            merged_qual.pop(pos, None)
        if mate_overlap_cpgs:
            qc["n_overlapping_mate_cpgs"] += mate_overlap_cpgs
        if mate_conflict_cpgs:
            qc["n_mate_conflict_cpgs"] += mate_conflict_cpgs

        if not merged:
            qc["n_fragments_dropped_after_merge"] += 1
            continue
        starts = [r[0].reference_start for r in records if r[0].reference_start is not None]
        ends = [r[0].reference_end for r in records if r[0].reference_end is not None]
        fragments.append(
            ReadMethylCalls(
                read_id=query_name,
                chrom=region.chrom,
                start=min(starts) if starts else region.start,
                end=max(ends) if ends else region.end,
                calls=merged,
                conversion_context=conversion_context,
                alignment_count=len(records),
                is_paired_fragment=is_paired,
                is_complete_pair=is_complete_pair,
                mate_overlap_cpgs=mate_overlap_cpgs,
                mate_conflict_cpgs=mate_conflict_cpgs,
            )
        )

    qc["n_fragments_with_calls"] = len(fragments)
    if return_qc:
        return fragments, dict(qc)
    return fragments


def amf_from_methyldackel_bedgraph(path: str, region: Region) -> Dict[str, float]:
    """Compute AMF from MethylDackel-like bedGraph.

    Exact AMF requires six columns:
    chrom, start, end, methylation_percent, methylated_count, unmethylated_count.
    Four-column bedGraph can only provide an unweighted site mean.
    """
    meth = 0.0
    unmeth = 0.0
    site_values: List[float] = []
    covered_sites = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("track") or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if chrom != region.chrom or end <= region.start or start >= region.end:
                continue
            covered_sites += 1
            if len(fields) >= 6:
                meth += float(fields[4])
                unmeth += float(fields[5])
            else:
                site_values.append(float(fields[3]) / 100.0)

    total = meth + unmeth
    result = {
        "bedgraph_covered_sites": float(covered_sites),
        "bedgraph_methylated_count": meth,
        "bedgraph_unmethylated_count": unmeth,
        "bedgraph_total_count": total,
        "bedgraph_amf": meth / total if total else math.nan,
        "bedgraph_site_mean_amf": (
            sum(site_values) / len(site_values) if site_values else math.nan
        ),
    }
    return result


def _safe_div(num: float, den: float) -> float:
    return num / den if den else math.nan


def _pattern_for_read(read: ReadMethylCalls, positions: Sequence[int]) -> Optional[str]:
    if all(pos in read.calls for pos in positions):
        return "".join("1" if read.calls[pos] else "0" for pos in positions)
    return None


def compute_methylation_features(
    reads: Sequence[ReadMethylCalls],
    region: Region,
    cpg_positions: Sequence[int],
    k: int = 4,
    min_cpgs_per_read_for_pdr: int = 4,
    min_covering_reads_per_cpg: int = 1,
    min_reads_for_pm: int = 1,
    min_reads_for_mhl: int = 1,
    min_possible_read_pairs: int = 1,
    min_common_cpgs_per_pair: int = 2,
) -> FeatureResult:
    positions = [pos for pos in sorted(set(cpg_positions)) if region.start <= pos < region.end]
    region_pos_set = set(positions)

    total_calls = 0
    methylated_calls = 0
    read_valid_for_pdr: List[bool] = []
    read_discordant: List[bool] = []
    for read in reads:
        values = [state for pos, state in read.calls.items() if pos in region_pos_set]
        if not values:
            read_valid_for_pdr.append(False)
            read_discordant.append(False)
            continue
        total_calls += len(values)
        methylated_calls += sum(values)
        valid = len(values) >= min_cpgs_per_read_for_pdr
        read_valid_for_pdr.append(valid)
        read_discordant.append(
            valid and any(v == 1 for v in values) and any(v == 0 for v in values)
        )

    valid_reads_for_pdr = sum(read_valid_for_pdr)
    discordant_reads = sum(read_discordant)
    cpg_pdr_values: List[float] = []
    for pos in positions:
        n_covering = 0
        n_discordant_covering = 0
        for idx, read in enumerate(reads):
            if pos not in read.calls:
                continue
            if not read_valid_for_pdr[idx]:
                continue
            n_covering += 1
            if read_discordant[idx]:
                n_discordant_covering += 1
        if n_covering >= min_covering_reads_per_cpg:
            cpg_pdr_values.append(n_discordant_covering / n_covering)

    pm_windows: List[Dict[str, object]] = []
    weighted_pm_sum = 0.0
    weighted_pm_den = 0
    weighted_me_sum = 0.0
    weighted_me_den = 0
    if len(positions) >= k:
        for idx in range(0, len(positions) - k + 1):
            window = positions[idx : idx + k]
            patterns = Counter()
            for read in reads:
                pattern = _pattern_for_read(read, window)
                if pattern is not None:
                    patterns[pattern] += 1
            n_valid = sum(patterns.values())
            if n_valid >= min_reads_for_pm and n_valid:
                pm = 1.0 - sum((count / n_valid) ** 2 for count in patterns.values())
                me = (-1.0 / k) * sum(
                    (count / n_valid) * math.log2(count / n_valid)
                    for count in patterns.values()
                    if count > 0
                )
                weighted_pm_sum += pm * n_valid
                weighted_pm_den += n_valid
                weighted_me_sum += me * n_valid
                weighted_me_den += n_valid
            else:
                pm = math.nan
                me = math.nan
            pm_windows.append(
                {
                    "window_index": idx,
                    "cpg_start": window[0],
                    "cpg_end": window[-1],
                    "k": k,
                    "n_valid_reads": n_valid,
                    "pm": pm,
                    "me": me,
                    "pattern_counts": ";".join(
                        "%s:%s" % (pattern, count)
                        for pattern, count in sorted(patterns.items())
                    ),
                }
            )

    possible_by_l: Dict[int, int] = defaultdict(int)
    full_methyl_by_l: Dict[int, int] = defaultdict(int)
    full_unmethyl_by_l: Dict[int, int] = defaultdict(int)
    if len(positions) >= k:
        for idx in range(0, len(positions) - k + 1):
            window = positions[idx : idx + k]
            for read in reads:
                pattern = _pattern_for_read(read, window)
                if pattern is None:
                    continue
                for length in range(1, k + 1):
                    for offset in range(0, k - length + 1):
                        substring = pattern[offset : offset + length]
                        possible_by_l[length] += 1
                        if substring == "1" * length:
                            full_methyl_by_l[length] += 1
                        if substring == "0" * length:
                            full_unmethyl_by_l[length] += 1

    mhl_by_length: List[Dict[str, object]] = []
    weighted_mhl_sum = 0.0
    weighted_mhl_den = 0
    weighted_uhl_sum = 0.0
    weighted_uhl_den = 0
    for length in range(1, k + 1):
        possible = possible_by_l[length]
        full_m = full_methyl_by_l[length]
        full_u = full_unmethyl_by_l[length]
        p_l = _safe_div(full_m, possible)
        p_l_unmethyl = _safe_div(full_u, possible)
        if possible >= min_reads_for_mhl:
            weighted_mhl_sum += length * p_l
            weighted_mhl_den += length
            weighted_uhl_sum += length * p_l_unmethyl
            weighted_uhl_den += length
        mhl_by_length.append(
            {
                "substring_length": length,
                "n_possible_substrings": possible,
                "n_fully_methylated_substrings": full_m,
                "n_fully_unmethylated_substrings": full_u,
                "p_l": p_l,
                "p_l_unmethyl": p_l_unmethyl,
            }
        )

    pairwise: List[Dict[str, object]] = []
    possible_pairs = 0
    discordant_pairs = 0
    hamming_sum = 0.0
    reads_with_region_calls = [
        read for read in reads if any(pos in region_pos_set for pos in read.calls)
    ]

    # CpG별 누적기 — O(n_positions × n_pairs) 이중 루프 제거
    cpg_pos_pairs: Dict[int, int] = defaultdict(int)
    cpg_pos_discordant: Dict[int, int] = defaultdict(int)
    cpg_pos_hamming_sum: Dict[int, float] = defaultdict(float)

    for read_a, read_b in itertools.combinations(reads_with_region_calls, 2):
        common = sorted(
            set(read_a.calls.keys()) & set(read_b.calls.keys()) & region_pos_set
        )
        if len(common) < min_common_cpgs_per_pair:
            continue
        distance = sum(1 for pos in common if read_a.calls[pos] != read_b.calls[pos])
        normalized = distance / len(common)
        is_discordant = distance > 0
        possible_pairs += 1
        if is_discordant:
            discordant_pairs += 1
        hamming_sum += normalized

        # 해당 쌍이 공유하는 CpG에 즉시 누적
        for pos in common:
            cpg_pos_pairs[pos] += 1
            if is_discordant:
                cpg_pos_discordant[pos] += 1
            cpg_pos_hamming_sum[pos] += normalized

        pairwise.append(
            {
                "read_a": read_a.read_id,
                "read_b": read_b.read_id,
                "n_common_cpgs": len(common),
                "n_discordant_cpgs": distance,
                "normalized_hamming": normalized,
            }
        )

    cpg_fdrp_values: List[float] = []
    cpg_qfdrp_values: List[float] = []
    for pos in positions:
        pairs_count = cpg_pos_pairs[pos]
        if pairs_count >= min_possible_read_pairs:
            cpg_fdrp_values.append(cpg_pos_discordant[pos] / pairs_count)
            cpg_qfdrp_values.append(cpg_pos_hamming_sum[pos] / pairs_count)

    fdrp = sum(cpg_fdrp_values) / len(cpg_fdrp_values) if cpg_fdrp_values else math.nan
    qfdrp = sum(cpg_qfdrp_values) / len(cpg_qfdrp_values) if cpg_qfdrp_values else math.nan

    summary = {
        "amf": _safe_div(methylated_calls, total_calls),
        "pdr": sum(cpg_pdr_values) / len(cpg_pdr_values) if cpg_pdr_values else math.nan,
        "pm": _safe_div(weighted_pm_sum, weighted_pm_den),
        "me": _safe_div(weighted_me_sum, weighted_me_den),
        "mhl": _safe_div(weighted_mhl_sum, weighted_mhl_den),
        "uhl": _safe_div(weighted_uhl_sum, weighted_uhl_den),
        "fdrp": fdrp,
        "qfdrp": qfdrp,
    }
    counts = {
        "n_reads_with_calls": len(reads_with_region_calls),
        "n_fragments_with_calls": len(reads_with_region_calls),
        "n_source_alignments": sum(read.alignment_count for read in reads_with_region_calls),
        "n_paired_fragments": sum(1 for read in reads_with_region_calls if read.is_paired_fragment),
        "n_complete_pair_fragments": sum(1 for read in reads_with_region_calls if read.is_complete_pair),
        "n_singleton_fragments": sum(
            1
            for read in reads_with_region_calls
            if read.is_paired_fragment and not read.is_complete_pair
        ),
        "n_overlapping_mate_cpgs": sum(read.mate_overlap_cpgs for read in reads_with_region_calls),
        "n_mate_conflict_cpgs": sum(read.mate_conflict_cpgs for read in reads_with_region_calls),
        "n_cpg_sites": len(positions),
        "n_methylated_cpg_observations": methylated_calls,
        "n_total_cpg_observations": total_calls,
        "n_valid_reads_for_pdr": valid_reads_for_pdr,
        "n_discordant_reads": discordant_reads,
        "n_cpg_sites_with_pdr": len(cpg_pdr_values),
        "n_pm_windows": len(pm_windows),
        "n_possible_read_pairs": possible_pairs,
        "n_discordant_read_pairs": discordant_pairs,
        "n_cpg_sites_with_fdrp": len(cpg_fdrp_values),
    }
    return FeatureResult(
        summary=summary,
        counts=counts,
        pm_windows=pm_windows,
        mhl_by_length=mhl_by_length,
        pairwise=pairwise,
    )


def read_calls_to_rows(reads: Sequence[ReadMethylCalls], cpg_positions: Sequence[int]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for read in reads:
        row: Dict[str, object] = {
            "read_id": read.read_id,
            "chrom": read.chrom,
            "start": read.start,
            "end": read.end,
            "conversion_context": read.conversion_context,
            "alignment_count": read.alignment_count,
            "is_paired_fragment": int(read.is_paired_fragment),
            "is_complete_pair": int(read.is_complete_pair),
            "mate_overlap_cpgs": read.mate_overlap_cpgs,
            "mate_conflict_cpgs": read.mate_conflict_cpgs,
        }
        for pos in cpg_positions:
            value = read.calls.get(pos)
            row["cpg_%s" % pos] = "" if value is None else value
        rows.append(row)
    return rows


def write_tsv(path: str, rows: Sequence[Dict[str, object]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not rows:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n")
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _summary_rows(result: FeatureResult) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for key, value in result.summary.items():
        rows.append({"name": key, "kind": "feature", "value": value})
    for key, value in result.counts.items():
        rows.append({"name": key, "kind": "count", "value": value})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", help="Input BAM/CRAM for read-level features")
    parser.add_argument("--bedgraph", help="MethylDackel bedGraph for AMF")
    parser.add_argument("--region", required=True, help="Region, e.g. chr1:100-200")
    parser.add_argument("--cpg-sites", help="BED/TSV with 0-based CpG cytosine positions")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument(
        "--call-source",
        choices=["tag", "c_to_t", "pe_strand"],
        default="tag",
        help=(
            "Methylation call source. pe_strand uses PE orientation: 99/147 "
            "calls C/T at CpG C, 163/83 calls G/A at CpG G. c_to_t is an alias."
        ),
    )
    parser.add_argument("--methyl-tag", default="XM")
    parser.add_argument(
        "--no-merge-paired-end",
        action="store_true",
        help="Keep each alignment as a separate haplotype instead of merging mates.",
    )
    parser.add_argument(
        "--overlap-policy",
        choices=["quality", "drop_conflict", "first"],
        default="quality",
        help=(
            "How to handle conflicting methylation calls in overlapping mates. "
            "'quality' uses MethylDackel-style quality-weighted resolution "
            "(default)."
        ),
    )
    parser.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Keep duplicate-marked alignments.",
    )
    parser.add_argument("--min-mapping-quality", type=int, default=10,
                        help="Minimum mapping quality (default 10).")
    parser.add_argument("--min-base-quality", type=int, default=5,
                        help="Minimum base quality for pe_strand/c_to_t mode (default 5).")
    parser.add_argument("--trim-5prime", type=int, default=10,
                        help="Trim N bases from sequencing 5' end for M-bias (default 10). "
                             "Used as fallback when per-read trims are not set.")
    parser.add_argument("--trim-3prime", type=int, default=10,
                        help="Trim N bases from sequencing 3' end for M-bias (default 10). "
                             "Used as fallback when per-read trims are not set.")
    parser.add_argument("--trim-5prime-r1", type=int, default=None,
                        help="Read1 5' M-bias trim (overrides --trim-5prime for R1).")
    parser.add_argument("--trim-3prime-r1", type=int, default=None,
                        help="Read1 3' M-bias trim (overrides --trim-3prime for R1).")
    parser.add_argument("--trim-5prime-r2", type=int, default=None,
                        help="Read2 5' M-bias trim (overrides --trim-5prime for R2).")
    parser.add_argument("--trim-3prime-r2", type=int, default=None,
                        help="Read2 3' M-bias trim (overrides --trim-3prime for R2).")
    parser.add_argument("--k", type=int, default=4)
    args = parser.parse_args()

    region = parse_region(args.region)
    cpg_positions = load_cpg_sites(args.cpg_sites, region) if args.cpg_sites else []

    if args.bam:
        reads, bam_qc = read_calls_from_bam(
            args.bam,
            region,
            cpg_positions=cpg_positions,
            call_source=args.call_source,
            methyl_tag=args.methyl_tag,
            merge_paired_end=not args.no_merge_paired_end,
            overlap_policy=args.overlap_policy,
            skip_duplicates=not args.keep_duplicates,
            min_mapping_quality=args.min_mapping_quality,
            min_base_quality=args.min_base_quality,
            trim_5prime=args.trim_5prime,
            trim_3prime=args.trim_3prime,
            trim_5prime_r1=args.trim_5prime_r1,
            trim_3prime_r1=args.trim_3prime_r1,
            trim_5prime_r2=args.trim_5prime_r2,
            trim_3prime_r2=args.trim_3prime_r2,
            return_qc=True,
        )
        if not cpg_positions:
            cpg_positions = sorted({pos for read in reads for pos in read.calls})
        result = compute_methylation_features(reads, region, cpg_positions, k=args.k)
        write_tsv(args.out_prefix + ".summary.tsv", _summary_rows(result))
        write_tsv(
            args.out_prefix + ".bam_merge_qc.tsv",
            [{"name": key, "value": value} for key, value in sorted(bam_qc.items())],
        )
        write_tsv(args.out_prefix + ".read_calls.tsv", read_calls_to_rows(reads, cpg_positions))
        write_tsv(args.out_prefix + ".pm_windows.tsv", result.pm_windows)
        write_tsv(args.out_prefix + ".mhl_by_length.tsv", result.mhl_by_length)
        write_tsv(args.out_prefix + ".pairwise.tsv", result.pairwise)

    if args.bedgraph:
        bedgraph = amf_from_methyldackel_bedgraph(args.bedgraph, region)
        write_tsv(
            args.out_prefix + ".bedgraph_amf.tsv",
            [{"name": key, "value": value} for key, value in bedgraph.items()],
        )


if __name__ == "__main__":
    main()

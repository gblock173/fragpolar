#!/usr/bin/env python3
"""Extract CpG methylation calls from BAM to Parquet.

Phase 1 of the two-phase pipeline: reads BAM once per chromosome with
batched streaming, extracts all CpG calls for marker regions, merges
mates, and writes a per-sample Parquet file.

The output Parquet can be used by the fast feature computation without
re-reading the BAM. This decouples I/O from computation and enables
fast re-runs with different parameters.

Usage:
  python extract_cpg_calls.py \
    --bam /path/to/sample.bam \
    --ref /path/to/hg38.fa \
    --markers /path/to/marker_table.tsv \
    --out /path/to/output/sample_id.cpg_calls.parquet \
    [--call-source pe_strand] \
    [--chromosomes chr1,chr2,...] \
    [--threads 4]
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:
    pa = None
    pq = None

try:
    import pysam
except ImportError:
    pysam = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from methyl_features import (
    Region,
    _fragment_conversion_context,
    _calls_from_alignment_by_context,
)


def cpg_positions_from_fasta(
    fasta: "pysam.FastaFile", chrom: str, start: int, end: int
) -> List[int]:
    """Find 0-based CpG cytosine positions in a region from the reference."""
    seq = fasta.fetch(chrom, start, end).upper()
    positions = []
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            positions.append(start + i)
    return positions


def load_markers(path: str, chromosomes: Optional[set] = None) -> Dict[str, List[dict]]:
    """Load markers grouped by chromosome, sorted by start."""
    import csv

    markers_by_chr: Dict[str, List[dict]] = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = row["chr"]
            if chromosomes and chrom not in chromosomes:
                continue
            markers_by_chr[chrom].append(
                {
                    "chrom": chrom,
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "marker_type": row["marker_type"],
                    "cpg_count": int(row["cpg_count"]),
                }
            )
    for chrom in markers_by_chr:
        markers_by_chr[chrom].sort(key=lambda m: m["start"])
    return dict(markers_by_chr)


def extract_chromosome(
    bam_path: str,
    ref_path: str,
    chrom: str,
    markers: List[dict],
    call_source: str = "pe_strand",
    min_mapping_quality: int = 10,
    min_base_quality: int = 10,
    skip_duplicates: bool = True,
    trim_5prime: int = 10,
    trim_3prime: int = 10,
    mbias_trim: Optional[dict] = None,
    threads: int = 1,
) -> List[dict]:
    """Stream BAM for one chromosome, extract CpG calls for all markers.

    Uses a sliding window over sorted markers to assign reads to
    overlapping marker regions without random seeks.

    Returns list of dicts (one per fragment per marker), ready for
    Parquet conversion.
    """
    if not markers:
        return []

    results = []
    marker_idx = 0
    n_markers = len(markers)
    active: List[Tuple[int, dict, Dict[str, list], List[int]]] = []

    t5_r1 = trim_5prime
    t3_r1 = trim_3prime
    t5_r2 = trim_5prime
    t3_r2 = trim_3prime

    fasta = pysam.FastaFile(ref_path)
    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        for aln in bam.fetch(chrom):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.is_qcfail:
                continue
            if skip_duplicates and aln.is_duplicate:
                continue
            if aln.mapping_quality < min_mapping_quality:
                continue

            ref_start = aln.reference_start
            ref_end = aln.reference_end or ref_start

            while marker_idx < n_markers and markers[marker_idx]["start"] <= ref_end:
                m = markers[marker_idx]
                cpg_pos = cpg_positions_from_fasta(
                    fasta, m["chrom"], m["start"], m["end"]
                )
                active.append((marker_idx, m, defaultdict(list), cpg_pos))
                marker_idx += 1

            new_active = []
            for midx, m, grouped, cpg_pos in active:
                if m["end"] <= ref_start:
                    _finalize_marker(
                        m,
                        grouped,
                        cpg_pos,
                        call_source,
                        min_base_quality,
                        t5_r1,
                        t3_r1,
                        t5_r2,
                        t3_r2,
                        mbias_trim,
                        results,
                    )
                else:
                    if ref_start < m["end"] and ref_end > m["start"]:
                        grouped[aln.query_name].append(aln)
                    new_active.append((midx, m, grouped, cpg_pos))
            active = new_active

    for midx, m, grouped, cpg_pos in active:
        _finalize_marker(
            m,
            grouped,
            cpg_pos,
            call_source,
            min_base_quality,
            t5_r1,
            t3_r1,
            t5_r2,
            t3_r2,
            mbias_trim,
            results,
        )

    fasta.close()
    return results


def _finalize_marker(
    marker: dict,
    grouped: Dict[str, list],
    cpg_positions: List[int],
    call_source: str,
    min_base_quality: int,
    t5_r1: int,
    t3_r1: int,
    t5_r2: int,
    t3_r2: int,
    mbias_trim: Optional[dict],
    results: List[dict],
):
    """Merge mates, extract CpG calls, and append to results."""
    region = Region(
        chrom=marker["chrom"], start=marker["start"], end=marker["end"]
    )
    marker_id = "%s:%d-%d" % (marker["chrom"], marker["start"], marker["end"])

    for query_name, alns in grouped.items():
        conversion_context = _fragment_conversion_context(alns)

        if call_source in {"pe_strand", "c_to_t"}:
            merged: Dict[int, int] = {}
            merged_qual: Dict[int, int] = {}
            conflicts = set()

            for a in alns:
                read_key = "R2" if a.is_read2 else "R1"
                if mbias_trim and (conversion_context, read_key) in mbias_trim:
                    t5, t3 = mbias_trim[(conversion_context, read_key)]
                elif a.is_read2:
                    t5, t3 = t5_r2, t3_r2
                else:
                    t5, t3 = t5_r1, t3_r1

                calls, quals = _calls_from_alignment_by_context(
                    a,
                    cpg_positions,
                    conversion_context,
                    min_base_quality=min_base_quality,
                    trim_5prime=t5,
                    trim_3prime=t3,
                )
                for pos, state in calls.items():
                    if not (region.start <= pos < region.end):
                        continue
                    q = quals.get(pos, 0)
                    if pos not in merged:
                        merged[pos] = state
                        merged_qual[pos] = q
                    else:
                        if merged[pos] == state:
                            merged_qual[pos] = max(merged_qual.get(pos, 0), q)
                        else:
                            existing_q = merged_qual.get(pos, 0)
                            if q > existing_q:
                                merged[pos] = state
                                merged_qual[pos] = q - existing_q
                            elif existing_q > q:
                                merged_qual[pos] = existing_q - q
                            else:
                                conflicts.add(pos)
            for pos in conflicts:
                merged.pop(pos, None)

            if not merged:
                continue
        else:
            continue

        is_paired = any(a.is_paired for a in alns)
        has_r1 = any(a.is_read1 for a in alns)
        has_r2 = any(a.is_read2 for a in alns)
        starts = [a.reference_start for a in alns if a.reference_start is not None]
        ends = [a.reference_end for a in alns if a.reference_end is not None]

        for pos, state in sorted(merged.items()):
            results.append(
                {
                    "marker_id": marker_id,
                    "fragment_id": query_name,
                    "chrom": marker["chrom"],
                    "pos": pos,
                    "call": state,
                    "fragment_start": min(starts) if starts else marker["start"],
                    "fragment_end": max(ends) if ends else marker["end"],
                    "conversion_context": conversion_context,
                    "is_paired": is_paired,
                    "is_complete_pair": is_paired and has_r1 and has_r2,
                    "n_alignments": len(alns),
                }
            )


PARQUET_SCHEMA = pa.schema(
    [
        ("marker_id", pa.string()),
        ("fragment_id", pa.string()),
        ("chrom", pa.string()),
        ("pos", pa.int32()),
        ("call", pa.int8()),
        ("fragment_start", pa.int32()),
        ("fragment_end", pa.int32()),
        ("conversion_context", pa.string()),
        ("is_paired", pa.bool_()),
        ("is_complete_pair", pa.bool_()),
        ("n_alignments", pa.int16()),
    ]
) if pa else None


class ParquetWriter:
    """Write CpG calls in per-chromosome row groups for efficient streaming."""

    def __init__(self, output_path: str):
        if pa is None:
            raise RuntimeError("pyarrow is required for Parquet output")
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        self.output_path = output_path
        self.tmp_path = output_path + ".tmp"
        self.writer = pq.ParquetWriter(
            self.tmp_path, PARQUET_SCHEMA, compression="snappy",
        )
        self.total_rows = 0

    def write_chromosome(self, records: List[dict]):
        if not records:
            return
        table = pa.Table.from_pylist(records, schema=PARQUET_SCHEMA)
        self.writer.write_table(table)
        self.total_rows += len(records)

    def close(self):
        self.writer.close()
        os.rename(self.tmp_path, self.output_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Input BAM/CRAM")
    parser.add_argument("--ref", required=True, help="Reference FASTA (hg38)")
    parser.add_argument("--markers", required=True, help="Marker table TSV")
    parser.add_argument("--out", required=True, help="Output Parquet path")
    parser.add_argument(
        "--call-source",
        choices=["pe_strand", "c_to_t"],
        default="pe_strand",
    )
    parser.add_argument(
        "--chromosomes",
        help="Comma-separated chromosomes (default: chr1-chr22)",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-bq", type=int, default=10)
    parser.add_argument("--trim-5prime", type=int, default=10)
    parser.add_argument("--trim-3prime", type=int, default=10)
    args = parser.parse_args()

    if args.chromosomes:
        chroms = set(args.chromosomes.split(","))
    else:
        chroms = {"chr%d" % i for i in range(1, 23)}

    mbias_trim = {
        ("top_ct", "R1"): (3, 0),
        ("top_ct", "R2"): (0, 20),
        ("bottom_ga", "R1"): (5, 3),
        ("bottom_ga", "R2"): (20, 0),
    }

    print("Loading markers from %s" % args.markers)
    markers_by_chr = load_markers(args.markers, chroms)
    total_markers = sum(len(v) for v in markers_by_chr.values())
    print("  %d markers across %d chromosomes" % (total_markers, len(markers_by_chr)))

    pq_writer = ParquetWriter(args.out)
    total_calls = 0
    t_start = time.time()

    for chrom in sorted(markers_by_chr.keys()):
        markers = markers_by_chr[chrom]
        print(
            "Processing %s: %d markers ..." % (chrom, len(markers)),
            flush=True,
        )
        t0 = time.time()
        records = extract_chromosome(
            args.bam,
            args.ref,
            chrom,
            markers,
            call_source=args.call_source,
            min_mapping_quality=args.min_mapq,
            min_base_quality=args.min_bq,
            skip_duplicates=True,
            trim_5prime=args.trim_5prime,
            trim_3prime=args.trim_3prime,
            mbias_trim=mbias_trim,
            threads=args.threads,
        )
        pq_writer.write_chromosome(records)
        total_calls += len(records)
        elapsed = time.time() - t0
        print(
            "  %s: %d CpG calls extracted (%.1fs)"
            % (chrom, len(records), elapsed)
        )

    pq_writer.close()
    total_elapsed = time.time() - t_start
    file_size = os.path.getsize(args.out) / (1024 * 1024)
    print(
        "\nTotal: %d CpG calls from %d markers in %.1f min (%.1f MB)"
        % (total_calls, total_markers, total_elapsed / 60, file_size)
    )


if __name__ == "__main__":
    main()

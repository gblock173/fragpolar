# fragpolar — fragment-level methylation states of plasma cfDNA in HCC

Analysis code for the manuscript on fragment polarisation (opposite strict-uniform fragments, BFV/WFV variance split,
FPI) of plasma cell-free DNA in hepatocellular carcinoma (HCC), chronic hepatitis and ageing.
원고(간세포암 혈장 cfDNA의 fragment 상태 양극화)의 결과와 그림을 만든 분석 코드이다.

> **Data are not included.** This repository contains code only: no sample identifiers, clinical variables, CpG calls,
> BAM files, result tables or source data. All inputs are given through `config/config.yaml`.
> 표본 ID, 임상정보, CpG call, BAM, 결과표, source data는 포함하지 않는다. 모든 입력 경로는 config로 지정한다.

## Contents / 구성

```
code_release/
├── README.md
├── environment.yml            main analysis environment (Python 3.11)
├── environment_py38.yml       BAM -> CpG-call extraction and PBMC step (Python 3.8, pysam 0.16)
├── requirements.txt, pyproject.toml
├── config/config.example.yaml paths, input files, column names, sample lists, parameters (placeholders only)
├── src/fragpolar/
│   ├── config.py              configuration loader (FRAGPOLAR_CONFIG or config/config.yaml)
│   ├── extract.py             CpG-call cleaning + per-marker sufficient statistics (K >= 4 fragments)
│   ├── engine.py              one-pass streaming engine (sufficient statistics, 10k-panel calls, production summaries)
│   ├── feature_definitions.py reference implementation of the 12 summaries on toy fragments (Fig. 1b)
│   ├── metrics.py             test-facing restatement of the core formulas
│   └── vendor/                fixed copy of the BAM -> CpG-call extractor
├── pipelines/                 analysis steps, run in numeric order
│   ├── 01_group_effects       marker-matrix group comparison, BFV decomposition, K/length adjustment, direction margin
│   ├── 02_train_rebuild       re-aggregation of Train samples from raw CpG calls; Train re-analysis
│   ├── 03_fpi_k               K-corrected FPI (FPI_K) and its approximation error
│   ├── 04_generators          synthetic generators G1/G2/G3, cosine similarity, polarisation rule, matched pairs
│   ├── 05_direct_fragments    sOUF removal, matched-dose controls, read split
│   ├── 06_tumour_hepatitis_age tumour-fraction dose-response, hepatitis, ageing
│   ├── 07_pbmc                PBMC genomic DNA extraction and comparison
│   ├── 08_component_origin    component-origin vectors vs tissue / normal-cell atlas / PBMC
│   └── 09_mixture_and_length  mixture-expected BFV, fragment length and control composition
├── figures/main/              main figures 1-6 (prep_* write source data, fig* draw; style_v8.py)
├── figures/supp/              supplementary figures S1-S18
├── scheduler/job_template.sh  generic batch-job template (PBS syntax; adapt to your scheduler)
├── tests/                     unit tests on synthetic input (pytest)
└── docs/CODE_STRUCTURE_KO.md  script-by-script description (Korean)
```

Each step writes to `<work_dir>/<package>/{results,data,work,figures,logs}` with package = group_effects, train_rebuild,
fpi_k, generators, direct_fragments, tumour_hepatitis_age, pbmc, component_origin, mixture_and_length.

## Installation / 설치

```bash
conda env create -f environment.yml && conda activate fragpolar
pip install -e .                             # optional; scripts add src/ to sys.path themselves
conda env create -f environment_py38.yml     # BAM -> CpG calls (vendor/extract_cpg_calls.py) and 07_pbmc
```
Figures use Liberation Sans (metric-compatible with Arial).

## Configuration / 설정

```bash
cp config/config.example.yaml config/config.yaml   # fill in every <PLACEHOLDER>
export FRAGPOLAR_CONFIG=$PWD/config/config.yaml
```
* `paths`: input roots and the output roots `work_dir`, `figure_dir`.
* `inputs`: every input file/folder by a generic name (relative to `data_dir`, or `{reference_dir}/...`, `{extraction_dir}/...`).
* `columns`: maps generic column names (e.g. `healthy_direction`) to the names in your files, so existing files need not be renamed.
* `sample_lists`: hepatitis samples analysed (the only hepatitis samples used anywhere), same-site HCC/healthy subsets, generator base samples.
* `params`: collection-year code of the same-site subset.

## Run order / 실행 순서

**[batch]** = heavy step; submit with `scheduler/job_template.sh` (`SCRIPT=...`, `ARGS=...`, array index `$ARRAY_INDEX`).

| Step | Command | Main outputs |
|---|---|---|
| 0 | `src/fragpolar/vendor/extract_cpg_calls.py --bam B --ref FASTA --markers MARKERS --out <cpg_calls_dir>/<id>.cpg_calls.parquet` **[batch]** | CpG calls |
| 1 | `01_group_effects/ge1_aggregate.py full 8` **[batch]** → `ge2_analysis.py` → `ge3_strict_subset_check.py` → `ge4_draft_figures.py` | ge_* tables (Test) |
| 2 | `02_train_rebuild/tr0_make_ids.py`; `tr1_extract.py --sid ID --force-raw` **[batch]**; `tr2_aggregate.py full 8` → `tr3_analysis.py` → `tr4_strict_subset_check.py` → `tr5_compare_legacy_raw.py`; parity `tr6_parity.py ID…`, `tr7_compare_strict_recompute.py` (after step 6) | ge_* tables (Train) |
| 3 | `03_fpi_k/fk1_aggregate.py 8`; `fk2_raw_invK.py ID` **[batch]**; `fk3_khist.py`; `fk5_jensen_raw.py ID ID ID`; `fk4_analysis.py`; `fk6_draft_figures.py` | fpik_* tables |
| 4 | `04_generators/onepass_engine.py --sid ID --path <cpg_calls> --outroot <work>/generators/data/engine` and `gen_panel_job.py --sid ID --engdir ... --b INDEX --out <work>/generators/data/panel/ID.parquet` **[batch, 20 base samples]**; `gen_observed.py 4`, `gen_observed_rebuilt.py 4` **[batch]**; `gen_analysis.py` (main), `gen_analysis_initial.py`, `gen_seed_consistency.py`, `gen_draft_figures.py` | gen_* tables, `gen_matched_pairs_rebuilt.tsv` |
| 5 | `05_direct_fragments/df0_initial_pairs.py`; `df1_extract.py TASK`, `df2_sample.py TASK` **[batch]**; `df0b_pairs_from_generators.py`; extra samples: `SAMPLES=<work>/direct_fragments/results_genpairs/df_extra_samples.tsv` with `df1_extract.py`/`df2_sample.py` **[batch]**; then `SAMPLES=.../df_union_samples.tsv PKL=<work>/direct_fragments/work/stratum_values_gen.pkl PAIRS=.../results_genpairs/df_pairs.tsv OUTD=<work>/direct_fragments/results_genpairs FIGD=<work>/direct_fragments/figures_genpairs` with `df3_stratum.py` → `df4_stats.py` → `df5_draft_figures.py` | df_* tables (main: results_genpairs) |
| 6 | `06_tumour_hepatitis_age/tha0_strict_extract.py ID` **[batch]**; `tha1_extract.py` → `tha2_tumour_fraction.py` → `tha3_hepatitis.py` → `tha3b_hepatitis_contrasts.py` → `tha4_ageing.py` → `tha5_draft_figures.py` | tf_*, hep_*, age_* tables |
| 7 | `07_pbmc/pb0_manifest.py`; `pb1_extract.py --id ID --bam BAM --outdir <work>/pbmc/data` **[batch, Python 3.8]**; `pb2_plasma_strict.py CHUNK N` **[batch]**; `pb3_qc_strict.py` → `pb4_summarize.py` → `pb5_analyze.py` → `pb5_analyze.py common` → `pb6_draft_figure.py`; parity `pb7_parity_check.py ID` | PBMC effects |
| 8 | `08_component_origin/co0_pat_tasks.py`, `co1_marker_cpg_ranges.py`; `co2_pat_extract.py PAT OUT`, `co3_tissue_suff.py ID` **[batch]**; `co4_marker_obs.py` → `co5_analysis.py full` → `co5_analysis.py t167` → `co6_draft_figures.py` | component cosines |
| 9 | `09_mixture_and_length/ml1_mixture_model.py 300` **[batch]** → `ml2_mixture_draft_figures.py`; `ml3_length_manifest.py`; `ml4_length.py TASK` **[batch]**; `ml5_length_analysis.py` → `ml6_length_draft_figures.py`; `ml7_length_tumour_fraction.py` (after step 6) | mix_*, len_* tables |
| 10 | `figures/main`: `prep_fig1.py`, `prep_fig1b_cards.py`, `prep_fig2.py`, `prep_fig3.py`, `prep_fig3_raster.py`, `fig456_prep.py`, `prep_fig6d.py`, `prep_fig4a_pairs.py`, `prep_fig4c_stats.py`, `fig5a_auc_make.py`, then `fig1.py` … `fig6.py`; `figures/supp`: `s1.py` … `s18.py` | figures and source data |

User-prepared sample tables read by batch steps (`pbmc/data/plasma_missing_strict.tsv`, `component_origin/data/tissue_tasks.tsv`,
identifier arguments) are described in docs/CODE_STRUCTURE_KO.md.

## Results and figures / 결과·그림 대응

| Manuscript item | Script(s) |
|---|---|
| Fig. 1a, 1c | `figures/main/prep_fig1.py`, `fig1.py` |
| Fig. 1b | `figures/main/prep_fig1b_cards.py` (`fragpolar.feature_definitions`) |
| Fig. 1d | `04_generators/gen_panel_job.py`, `gen_analysis.py` → `prep_fig1.py` |
| Fig. 2a–d | `01_group_effects/ge2_analysis.py` (Test), `02_train_rebuild/tr3_analysis.py` (Train), `04_generators/gen_analysis.py`, `03_fpi_k/fk4_analysis.py` → `prep_fig2.py` |
| Fig. 2e | `06_tumour_hepatitis_age/tha1_extract.py`, `tha2_tumour_fraction.py` |
| Fig. 3a | `05_direct_fragments/df2_sample.py` → `prep_fig3_raster.py` |
| Fig. 3b | `09_mixture_and_length/ml1_mixture_model.py` → `prep_fig3.py` |
| Fig. 3c–e | `05_direct_fragments/df1_extract.py` … `df4_stats.py` → `prep_fig3.py` |
| Fig. 3f | `04_generators/gen_analysis.py` → `prep_fig3.py` |
| Fig. 4a | `fig456_prep.py`, `prep_fig4a_pairs.py` (`08_component_origin/co3_tissue_suff.py`) |
| Fig. 4b | `fig456_prep.py` (earlier analysis table, `inputs.prior_cnv_effect_summary`) |
| Fig. 4c | `06_tumour_hepatitis_age/tha2_tumour_fraction.py` → `fig456_prep.py`, `prep_fig4c_stats.py` |
| Fig. 4d | `08_component_origin/co5_analysis.py` → `fig456_prep.py` |
| Fig. 5a | `figures/main/fig5a_auc_make.py` |
| Fig. 5b–c | `06_tumour_hepatitis_age/tha1_extract.py`, `tha3b_hepatitis_contrasts.py` → `fig456_prep.py`, `fig5.py` |
| Fig. 6 | `06_tumour_hepatitis_age/tha4_ageing.py`, `tha3b_hepatitis_contrasts.py` → `fig456_prep.py` (a–c), `prep_fig6d.py` (d), `fig6.py` |
| Supp. Fig. 1, 2, 11 | `s1.py`, `s2.py`, `s11.py` (stored or earlier-analysis tables) |
| Supp. Fig. 3, 5, 6 | `gen_analysis.py`, `ge2_analysis.py` / `tr3_analysis.py` → `s3.py`, `s5.py`, `s6.py` |
| Supp. Fig. 4 | `03_fpi_k/*` → `s4.py` |
| Supp. Fig. 7, 8 | `ml1_mixture_model.py`, `ml4_length.py`, `ml5_length_analysis.py` → `s7.py`, `s8.py` |
| Supp. Fig. 9 | `07_pbmc/*` → `s9.py` |
| Supp. Fig. 10 | `co5_analysis.py`, `tha3b_hepatitis_contrasts.py` → `s10.py` |
| Supp. Fig. 12, 13 | `tha4_ageing.py` → `s12.py`; `s13.py` |
| Supp. Fig. 14 | `gen_analysis.py` → `s14.py` (full comparison, agreement with matched pairs) |
| Supp. Fig. 15 | `ge2_analysis.py` / `tr3_analysis.py`, `ml5_length_analysis.py` → `s15.py` |
| Supp. Fig. 16 | `ge2_analysis.py` / `tr3_analysis.py`, `fk4_analysis.py`, `ml5_length_analysis.py`, `ml7_length_tumour_fraction.py` → `s16.py` |
| Supp. Fig. 17 | `tha2_tumour_fraction.py` → `s17.py` |
| Supp. Fig. 18 | `tha4_ageing.py` → `s18.py` |

## Tests / 테스트

```bash
python -m pytest -q tests      # synthetic input only
```

## Repository

https://github.com/gblock173/fragpolar

## Citation and licence

Please cite the accompanying article (citation and DOI to be added after publication).

The code is distributed under the [MIT License](LICENSE). Copyright (c) 2026 fragpolar contributors.
External software and datasets retain their own licences and access conditions.

"""fragpolar: fragment-level methylation state analysis of plasma cfDNA (analysis code accompanying the manuscript).

Modules
-------
config               configuration loader (paths, sample lists, parameters)
extract              CpG-call cleaning and per-marker sufficient statistics (raw_extract / aggregate / check_metrics)
engine               one-pass streaming engine over a sample's CpG calls (sufficient statistics, 10k-panel calls,
                     production feature summaries via gold_summary)
feature_definitions  reference implementation of the 12 fragment summaries on toy fragments (Fig. 1b cards)
metrics              small pure functions shared by the tests (variance decomposition, FPI, sOUF, removal, mixture)
vendor/              fixed copy of the BAM -> CpG-call extractor (extract_cpg_calls.py, methyl_features.py)
"""
__version__ = '1.0.0'

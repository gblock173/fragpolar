"""Configuration loader: placeholders count as unset, sample lists are read from files, outputs go under work_dir."""
import importlib
import pytest


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    ids = tmp_path / 'ids.txt'; ids.write_text('# comment\nS001\n\nS002\n')
    f = tmp_path / 'c.yaml'
    f.write_text(f"paths:\n  data_dir: {tmp_path}/data/\n  reference_dir: {tmp_path}/ref\n  work_dir: {tmp_path}/work\n  cpg_calls_dir: <PATH_TO_CPG_CALLS>\n"
                 f"inputs:\n  samples_table: cohort_tables/samples.tsv\n  healthy_direction_npy: \"{{reference_dir}}/hd.npy\"\n"
                 f"columns:\n  healthy_direction: orig_direction_col\n  hof_feature: orig_hof\n"
                 f"sample_lists:\n  hepatitis_included: {ids}\nparams:\n  same_site_era: 22\n")
    monkeypatch.setenv('FRAGPOLAR_CONFIG', str(f))
    from fragpolar import config
    importlib.reload(config)
    return config, tmp_path


def test_paths_and_lists(cfg):
    C, t = cfg
    assert C.DATA == f'{t}/data' and C.work('group_effects', 'results') == f'{t}/work/group_effects/results'
    assert C.ids('hepatitis_included') == ['S001', 'S002'] and C.param('same_site_era') == 22


def test_placeholder_is_unset(cfg):
    C, _ = cfg
    with pytest.raises(KeyError):
        C.CPG_CALLS
    assert C.path('gallery_values_tsv', required=False) is None


def test_pipeline_dir(cfg):
    C, _ = cfg
    assert C.pipeline_dir('group_effects').endswith('pipelines/01_group_effects')


def test_inputs_columns_features(cfg):
    import pandas as pd
    C, t = cfg
    assert C.inp('samples_table') == f'{t}/data/cohort_tables/samples.tsv'
    assert C.inp('healthy_direction_npy') == f'{t}/ref/hd.npy'
    assert C.col('healthy_direction') == 'orig_direction_col' and C.col('unknown_generic') == 'unknown_generic'
    assert C.featname('hof') == 'orig_hof' and C.featname('pdr') == 'pdr'
    d = C.unmap_features(pd.DataFrame({'feature': ['orig_hof', 'pdr']}))
    assert d.feature.tolist() == ['hof', 'pdr']

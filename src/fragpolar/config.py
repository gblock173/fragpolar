"""Central configuration loader for the fragpolar analysis code.

Every path, sample list and site-specific parameter used by the pipelines is read from one YAML file:

    1. the file named by the environment variable FRAGPOLAR_CONFIG, or
    2. <repo>/config/config.yaml (copy config/config.example.yaml and fill in the placeholders).

The module is imported by every script as ``from fragpolar import config as C``.
Attributes are resolved lazily (PEP 562), so importing the module never fails; a missing or unfilled
key raises a clear error only when a script actually needs it.

Main attributes
---------------
C.DATA          analysis data directory (marker matrices, cohort tables, strata definition, ...)
C.REFERENCE     reference directory (marker z-score references, healthy direction, strata markers, cohort matrices)
C.EXTRACTION    extraction directory (marker regions, extraction manifests)
C.CPG_CALLS     directory of per-sample <sample>.cpg_calls.parquet files
C.WORK          root of analysis outputs; each pipeline writes to C.work('<package>')
C.FIG           root of figure outputs (Figure_*.pdf, source_data/, supp/)
C.path(key)     any key of the ``paths`` section
C.inp(name)     an input file/folder of the ``inputs`` section (relative to data_dir unless absolute; {reference_dir} etc. expand)
C.col(name)     a column name of the ``columns`` section (maps generic names to the column names of the user's files)
C.featname(f)   file stem / table label of a feature (only 'hof' is mapped, via columns.hof_feature)
C.ids(key)      list of sample identifiers read from a file of the ``sample_lists`` section (one ID per line)
C.param(key)    any key of the ``params`` section
C.pipeline_dir(pkg)  directory of a pipeline package inside this repository (for cross-package imports)
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
_CFG = None

# analysis package name -> folder under pipelines/ (output sub-folder names under work_dir keep the package name)
PIPELINES = {
    'group_effects': '01_group_effects', 'train_rebuild': '02_train_rebuild', 'fpi_k': '03_fpi_k', 'generators': '04_generators',
    'direct_fragments': '05_direct_fragments', 'tumour_hepatitis_age': '06_tumour_hepatitis_age', 'pbmc': '07_pbmc',
    'component_origin': '08_component_origin', 'mixture_and_length': '09_mixture_and_length',
}

# attribute name -> key in the 'paths' section
_ALIASES = {
    'DATA': 'data_dir', 'REFERENCE': 'reference_dir', 'EXTRACTION': 'extraction_dir', 'CPG_CALLS': 'cpg_calls_dir',
    'HEP_CPG_CALLS': 'hepatitis_cpg_calls_dir', 'WORK': 'work_dir', 'FIG': 'figure_dir',
}


def config_file():
    return Path(os.environ.get('FRAGPOLAR_CONFIG', REPO / 'config' / 'config.yaml'))


def load():
    """Return the parsed configuration dictionary (cached)."""
    global _CFG
    if _CFG is None:
        f = config_file()
        if not f.exists():
            raise FileNotFoundError(f'configuration file not found: {f}. Copy config/config.example.yaml to config/config.yaml '
                                    f'(or set FRAGPOLAR_CONFIG) and fill in the placeholders.')
        import yaml
        with open(f) as fh:
            _CFG = yaml.safe_load(fh) or {}
    return _CFG


def _get(section, key, required=True):
    v = (load().get(section) or {}).get(key)
    if isinstance(v, str) and v.startswith('<') and v.endswith('>'):
        v = None                                   # unfilled placeholder
    if v is None and required:
        raise KeyError(f"config: '{section}.{key}' is not set in {config_file()}")
    return v


def path(key, required=True):
    """Value of paths.<key> as a string without a trailing slash (None if optional and unset)."""
    v = _get('paths', key, required)
    return None if v is None else str(v).rstrip('/')


def param(key, default=None):
    v = _get('params', key, required=default is None)
    return default if v is None else v


def ids(key):
    """Sample identifiers listed one per line in the file sample_lists.<key> (blank lines and '#' comments ignored)."""
    f = _get('sample_lists', key)
    with open(f) as fh:
        return [l.strip() for l in fh if l.strip() and not l.lstrip().startswith('#')]


def ids_file(key):
    return str(_get('sample_lists', key))


def _expand(v):
    paths = {k: str(x).rstrip('/') for k, x in (load().get('paths') or {}).items() if x is not None}
    return str(v).format(**paths)


def inp(name):
    """Input file or folder inputs.<name>; relative values are taken under paths.data_dir."""
    v = _expand(_get('inputs', name)).rstrip('/')
    return v if os.path.isabs(v) else os.path.join(path('data_dir'), v)


def col(name):
    """Column name columns.<name> (defaults to the generic name itself)."""
    v = (load().get('columns') or {}).get(name)
    return name if v is None else str(v)


def featname(f):
    """Feature stem used in reference/cohort-matrix file names and earlier tables ('hof' -> columns.hof_feature)."""
    return col('hof_feature') if f == 'hof' else f


def unmap_features(df, column='feature'):
    """Rename the user's HOF feature label back to the generic 'hof' in a table read from an earlier analysis."""
    if column in df.columns and featname('hof') != 'hof':
        df = df.assign(**{column: df[column].replace({featname('hof'): 'hof'})})
    return df


def work(pkg, *sub):
    """Output directory of an analysis package (string), e.g. work('group_effects', 'results')."""
    return os.path.join(path('work_dir'), pkg, *sub)


def cpg_calls(sid):
    """Path of the per-sample CpG-call parquet in cpg_calls_dir."""
    return os.path.join(path('cpg_calls_dir'), f'{sid}.cpg_calls.parquet')


def pipeline_dir(pkg):
    return str(REPO / 'pipelines' / PIPELINES[pkg])


def __getattr__(name):                             # lazy module attributes (Python >= 3.7)
    if name in _ALIASES:
        return path(_ALIASES[name])
    raise AttributeError(name)

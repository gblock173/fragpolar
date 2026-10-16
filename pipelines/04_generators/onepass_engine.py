"""generator step entry point of the one-pass engine (cfDNA, explicit input path / output root).

The engine code itself lives in src/fragpolar/engine.py (single copy shared with FPI_K step and CD); this file keeps the original
command line:  python onepass_engine.py --sid SAMPLE --path <cpg_calls.parquet> --outroot <dir> [--limit N]
and the names imported by gen_lib (gold_summary, trim_keep, marker_table).
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar.engine import *               # noqa: E402,F401,F403
from fragpolar.engine import main            # noqa: E402

if __name__ == '__main__':
    main()

from pkg_resources import get_distribution

try:
    __version__ = get_distribution("cte").version
except:
    __version__ = "local"


__all__ = [
    "amplicon_scheme",
    "built_in_data",
    "iupac",
    "msa",
    "one_run_evaluator",
    "tasks",
    "utils",
    "variants",
    "varifier_tools",
]

from cte import *

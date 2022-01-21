from pkg_resources import get_distribution

try:
    __version__ = get_distribution("cte").version
except:
    __version__ = "local"


__all__ = [
    "built_in_data",
    "one_run_evaluator",
    "primers",
    "tasks",
    "utils",
    "varifier_tools",
]

from cte import *

from cte import multi_run_evaluator

import os


def run(options):
    cte.multi_run_evaluator.MultiRunEvaluator(options.manifest_tsv)
    cte.run()

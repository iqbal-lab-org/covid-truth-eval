import csv

from cte import one_run_evaluator


class MultiRunEvaluator:
    def __init__(self, data_tsv):
        self.data_tsv = data_tsv

    def run(self):
        raise NotImplementedError()

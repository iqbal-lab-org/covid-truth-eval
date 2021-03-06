import logging
import os

from cte import one_run_evaluator, utils


def run(options):
    one_run_evaluator.eval_one_fasta(
        options.outdir,
        options.fasta_to_eval,
        options.ref_fasta,
        options.primers,
        truth_vcf=options.truth_vcf,
        truth_fasta=options.truth_fasta,
        force=options.force,
        debug=options.debug,
    )

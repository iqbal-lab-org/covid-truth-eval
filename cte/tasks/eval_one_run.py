import logging
import os

from cte import one_run_evaluator, utils


def run(options):
    if options.force and os.path.exists(options.outdir):
        logging.info(
            f"Option --force used. Deleting existing output directory {options.outdir}"
        )
        utils.syscall(f"rm -rf {options.outdir}")
    one_run_evaluator.eval_one_fasta(
        options.outdir,
        options.fasta_to_eval,
        options.ref_fasta,
        options.truth_vcf,
        options.primers,
    )

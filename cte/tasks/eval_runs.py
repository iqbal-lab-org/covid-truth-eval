from cte import multi_run_evaluator, utils

import logging
import os


def run(options):
    if options.force and os.path.exists(options.outdir):
        logging.info(
            f"Option --force used. Deleting existing output directory {options.outdir}"
        )
        utils.syscall(f"rm -rf {options.outdir}")
    multi_run_evaluator.evaluate_runs(
        options.manifest_tsv, options.outdir, options.ref_fasta, debug=options.debug
    )

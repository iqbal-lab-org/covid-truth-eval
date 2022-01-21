from cte import one_run_evaluator, utils

import os


def run(options):
    if options.force:
        utils.syscall(f"rm -rf {options.outdir}")
    os.mkdir(options.outdir)
    truth_files_dir = os.path.join(options.outdir, "Truth_files")
    evaluator = one_run_evaluator.OneSampleEvaluator(
        options.ref_fasta, options.truth_vcf, truth_files_dir
    )
    eval_dir = os.path.join(options.outdir, "Results")
    evaluator.evaluate_one_fasta(
        eval_dir,
        options.fasta_to_eval,
        options.primers,
    )

import logging
import os

from cte import amplicon_scheme, msa, utils


def eval_one_fasta(
    outdir,
    fasta_to_eval,
    ref_fasta,
    primers_name_or_tsv,
    truth_vcf=None,
    truth_fasta=None,
    debug=False,
    force=False,
):
    if truth_vcf is None and truth_fasta is None:
        raise Exception("Must provide truth_vcf or truth_fasta. Cannot continue")

    if force:
        logging.info(
            f"--force option used, deleting output directory if it already exists {outdir}"
        )
        utils.syscall(f"rm -rf {outdir}")
    elif os.path.exists(outdir):
        raise Exception(f"Ouptut directory already exists. Cannot continue. {outdir}")
    logging.info(f"Start evaluating {fasta_to_eval}")
    logging.info(f"Using reference FASTA: {ref_fasta}")

    stats_summary_tsv_out = os.path.join(outdir, "results.tsv")
    stats_summary_json_out = os.path.join(outdir, "results.json")
    per_position_tsv_out = os.path.join(outdir, "per_position.tsv")
    msa_dir = os.path.join(outdir, "MSA")

    amp_scheme = amplicon_scheme.AmpliconScheme(primers_name_or_tsv)
    os.mkdir(outdir)
    if truth_fasta is None:
        logging.info(f"Using truth VCF to make truth sequence: {truth_vcf}")
        truth_fasta = os.path.join(outdir, "truth.fasta")
        utils.apply_variants_to_genome(truth_vcf, truth_fasta, ref_fasta=ref_fasta)
        user_provided_truth_fasta = False
    else:
        logging.info(f"Using provided truth FASTA: {truth_fasta}")
        user_provided_truth_fasta = True
    multi_aln = msa.Msa(ref_fasta, truth_fasta, fasta_to_eval)
    logging.info("Making multiple sequence alignment (ref vs truth vs seq to evaluate)")
    multi_aln.run_msa(msa_dir)
    logging.info("Gathering stats from MSA")
    multi_aln.add_eval_ends_missing()
    multi_aln.make_coords_lookups()
    if truth_vcf is not None:
        multi_aln.add_truth_dropped_amps(truth_vcf)
    multi_aln.add_eval_dropped_amps(amp_scheme)
    multi_aln.gather_stats(amp_scheme, per_pos_tsv=per_position_tsv_out)
    multi_aln.write_stats_summary_tsv(stats_summary_tsv_out)
    multi_aln.write_stats_summary_json(stats_summary_json_out)
    if debug:
        logging.info("Debug option used, so not cleaning up intermediate files")
    else:
        logging.info("Deleting intermediate files")
        utils.syscall(f"rm -r {msa_dir}")
        if not user_provided_truth_fasta:
            utils.syscall(f"rm {truth_fasta}")
    logging.info(f"Finished evaluating {fasta_to_eval}")
    return multi_aln.stats

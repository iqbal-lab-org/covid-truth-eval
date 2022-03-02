import logging
import os

from cte import amplicon_scheme, msa, utils


def eval_one_fasta(
    outdir,
    fasta_to_eval,
    ref_fasta,
    truth_vcf,
    primers_name_or_tsv,
    debug=False,
    force=False,
):
    if force:
        logging.info(
            f"--force option used, deleting output directory if it already exists {outdir}"
        )
        utils.syscall(f"rm -rf {outdir}")
    elif os.path.exists(outdir):
        raise Exception(f"Ouptut directory already exists. Cannot continue. {outdir}")
    logging.info(f"Start evaluating {fasta_to_eval}")
    logging.info(f"Using reference FASTA: {ref_fasta}")
    logging.info(f"Using truth VCF: {truth_vcf}")

    stats_summary_tsv_out = os.path.join(outdir, "results.tsv")
    stats_summary_json_out = os.path.join(outdir, "results.json")
    per_position_tsv_out = os.path.join(outdir, "per_position.tsv")
    msa_dir = os.path.join(outdir, "MSA")
    truth_fasta = os.path.join(outdir, "truth.fasta")

    amp_scheme = amplicon_scheme.AmpliconScheme(primers_name_or_tsv)
    os.mkdir(outdir)
    utils.apply_variants_to_genome(truth_vcf, truth_fasta, ref_fasta=ref_fasta)
    multi_aln = msa.Msa(ref_fasta, truth_fasta, fasta_to_eval)
    logging.info("Making multiple sequence alignment (ref vs truth vs seq to evaluate)")
    multi_aln.run_msa(msa_dir)
    logging.info("Gathering stats from MSA")
    multi_aln.add_eval_ends_missing()
    multi_aln.make_coords_lookups()
    multi_aln.add_truth_dropped_amps(truth_vcf)
    multi_aln.add_eval_dropped_amps(amp_scheme)
    multi_aln.gather_stats(amp_scheme, per_pos_tsv=per_position_tsv_out)
    multi_aln.write_stats_summary_tsv(stats_summary_tsv_out)
    multi_aln.write_stats_summary_json(stats_summary_json_out)
    logging.info(f"Finished evaluating {fasta_to_eval}")
    return multi_aln.stats

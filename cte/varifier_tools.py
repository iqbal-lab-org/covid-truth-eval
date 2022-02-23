import logging
import os

import cluster_vcf_records

from cte import variants, utils


def varifier_vcf_to_counts_and_errors(
    infile, primers, first_primer_start, last_primer_end, fp_or_fn
):
    errors = []
    stats = {
        "snp": {"TP": 0, fp_or_fn: 0},
        "indel": {"TP": 0, fp_or_fn: 0},
        "primer_snp": {"TP": 0, fp_or_fn: 0},
        "primer_indel": {"TP": 0, fp_or_fn: 0},
    }
    header, records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(infile)
    for record in records:
        if not first_primer_start <= record.POS <= last_primer_end:
            logging.warning(f"Skipping VCF record because outside primer range: {record}")
            continue
        if record.FORMAT["VFR_IN_MASK"] == "1":
            logging.info(f"Skipping record because is in mask: {record}")
            continue
        assert record.FORMAT["VFR_IN_MASK"] == "0"

        if record.FORMAT["VFR_RESULT"] == "TP":
            tp_or_fp = "TP"
        else:
            tp_or_fp = fp_or_fn
            errors.append(str(record))

        assert len(record.ALT) == 1
        if len(record.REF) == 1 == len(record.ALT[0]):
            snp_or_indel = "snp"
        else:
            snp_or_indel = "indel"

        stats[snp_or_indel][tp_or_fp] += 1

        for i in range(record.POS, record.ref_end_pos() + 1):
            if i in primers:
                stats["primer_" + snp_or_indel][tp_or_fp] += 1
                break

    return stats, errors


def varifier_outdir_to_stats(dirname, primers, first_primer_start, last_primer_end):
    precision_vcf = os.path.join(dirname, "precision.vcf")
    stats, fp_errors = varifier_vcf_to_counts_and_errors(
        precision_vcf, primers, first_primer_start, last_primer_end, "FP"
    )
    recall_vcf = os.path.join(dirname, "recall", "recall.vcf.masked.vcf")
    recall_stats, fn_errors = varifier_vcf_to_counts_and_errors(
        recall_vcf, primers, first_primer_start, last_primer_end, "FN"
    )
    for key, d in stats.items():
        d["FN"] = recall_stats[key]["FN"]
    stats["Errors"] = {"FP": fp_errors, "FN": fn_errors}
    return stats


def load_msa_file(msa_file):
    with open(msa_file) as f:
        msa = [x.rstrip() for x in f]
        if len(msa) != 2 or len(msa[0]) != len(msa[1]):
            raise Exception(f"Error loading MSA file {self.msa_file}")

    return msa[0], msa[1]


def make_msa_and_vcf(ref_fasta, query_fasta, min_ref_coord, max_ref_coord, outdir, debug=False):
    """Makes a multiple alignment between the ref and query fasta,
    and a VCF of inferred variants bewteen the two.
    returns them as a tuple: (ref_msa, qry_msa, vcf)"""
    # Uses varifier to compare the fasta to evaluate against the reference
    # genome. Runs varifier in the mode that makes a global alignment and
    # writes the MSA to a file. Also use the option --use_non_acgt, which
    # means the file 01.merged.vcf has ALL variants between the sequences,
    # including those involving Ns and other ambiguous IUPAC codes like
    # R (=A or G), Y (=C or T) ... etc.
    utils.syscall(
        f"varifier make_truth_vcf --use_non_acgt --global_align --global_align_min_coord {min_ref_coord} --global_align_max_coord {max_ref_coord} {query_fasta} {ref_fasta} {outdir}"
    )
    vcf = os.path.join(outdir, "01.merged.vcf")
    msa_file = os.path.join(outdir ,"04.msa")
    if not (os.path.exists(vcf) and os.path.exists(msa_file)):
        raise Exception("Error making VCF and/or MSA. Cannot continue")

    msa_ref, msa_qry = load_msa_file(msa_file)
    return msa_ref, msa_qry, vcf

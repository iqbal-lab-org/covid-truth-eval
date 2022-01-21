import logging
import os

import cluster_vcf_records


def varifier_vcf_to_counts(
    infile, primers, first_primer_start, last_primer_end, fp_or_fn
):
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

    return stats


def varifier_outdir_to_stats(dirname, primers, first_primer_start, last_primer_end):
    precision_vcf = os.path.join(dirname, "precision.vcf")
    stats = varifier_vcf_to_counts(
        precision_vcf, primers, first_primer_start, last_primer_end, "FP"
    )
    recall_vcf = os.path.join(dirname, "recall", "recall.vcf.masked.vcf")
    recall_stats = varifier_vcf_to_counts(
        recall_vcf, primers, first_primer_start, last_primer_end, "FN"
    )
    for key, d in stats.items():
        d["FN"] = recall_stats[key]["FN"]
    return stats

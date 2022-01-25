import copy
import json
import os
import pytest
import random

from cte import one_run_evaluator, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "one_run_evaluator")


def load_vcf_except_headers(filename):
    with open(filename) as f:
        return [x for x in f if not x.startswith("#")]


def vcfs_same_exclude_headers(file1, file2):
    lines1 = load_vcf_except_headers(file1)
    lines2 = load_vcf_except_headers(file2)
    return lines1 == lines2


def test_eval_one_fasta():
    outdir = "tmp.eval_one_fasta"
    ref_fasta = f"{outdir}.ref.fa"
    fasta_to_eval = f"{outdir}.to_eval.fa"
    truth_vcf = f"{outdir}.truth.vcf"
    primers_tsv = f"{outdir}.primers.tsv"
    utils.syscall(f"rm -rf {outdir}*")

    # Making these variants:
    # 11 C -> T. TP, but before first amplicon so should get ignored
    # 101 A -> G. TP SNP
    # 201 G -> GA. TP indel
    # 301 CA -> C. TP indel
    # 601 T -> A. SNP in a primer. But test fasta has T -> G, so FP and FN
    # 801 A -> C, but is a HET call. test fasta has A -> G, which should get
    # ignored instead of being called as wrong. Position should get added
    # to the mask BED file.
    # Adding those up, we get expected counts:
    expect = {
        "snp": {"TP": 1, "FP": 1, "FN": 1},
        "indel": {"TP": 2, "FP": 0, "FN": 0},
        "primer_snp": {"TP": 0, "FP": 1, "FN": 1},
        "primer_indel": {"TP": 0, "FP": 0, "FN": 0},
        "Errors": {
            "FP": [
                "ref\t601\t4\tT\tG\t.\tPASS\t.\tGT:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT\t1/1:1:1:1:1:0:0.0:0:FP"
            ],
            "FN": [
                "ref\t601\t.\tT\tA\t42.42\tPASS\t.\tGT:VFR_ED_RA:VFR_ED_TR:VFR_ED_TA:VFR_ALLELE_LEN:VFR_ALLELE_MATCH_COUNT:VFR_ALLELE_MATCH_FRAC:VFR_IN_MASK:VFR_RESULT\t1/1:1:1:1:1:0:0.0:0:FP"
            ],
        },
    }

    with open(truth_vcf, "w") as f:
        print(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1",
            "ref\t11\t.\tC\tT\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t101\t.\tA\tG\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t201\t.\tG\tGA\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t301\t.\tCA\tC\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t601\t.\tT\tA\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t801\t.\tA\tC\t42.42\tHET\t.\tGT\t0/1",
            sep="\n",
            file=f,
        )

    random.seed(42)
    ref_seq = random.choices(["A", "C", "G", "T"], k=1000)
    eval_seq = copy.copy(ref_seq)
    ref_seq[10] = "C"
    ref_seq[100] = "A"
    ref_seq[200] = "G"
    ref_seq[300] = "C"
    ref_seq[301] = "A"
    ref_seq[600] = "T"
    ref_seq[800] = "A"
    with open(ref_fasta, "w") as f:
        print(">ref", file=f)
        print("".join(ref_seq), file=f)
    eval_seq[10] = "T"
    eval_seq[100] = "G"
    eval_seq[200] = "GA"
    eval_seq[300] = "C"
    eval_seq[301] = ""
    eval_seq[600] = "G"
    eval_seq[800] = "G"
    with open(fasta_to_eval, "w") as f:
        print(">to_eval", file=f)
        print("".join(eval_seq), file=f)

    with open(primers_tsv, "w") as f:
        print(
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
            sep="\t",
            file=f,
        )
        print("A1", "A1_left", "left", "".join(ref_seq[20:30]), 20, sep="\t", file=f)
        print(
            "A1", "A1_right", "right", "".join(ref_seq[590:620]), 590, sep="\t", file=f
        )
        print("A2", "A2_left", "left", "".join(ref_seq[610:630]), 630, sep="\t", file=f)
        print(
            "A2", "A2_right", "right", "".join(ref_seq[950:970]), 950, sep="\t", file=f
        )

    got = one_run_evaluator.eval_one_fasta(
        outdir,
        fasta_to_eval,
        ref_fasta,
        truth_vcf,
        primers_tsv,
    )
    print(got)
    assert got == expect
    with open(os.path.join(outdir, "Results", "results.json")) as f:
        got = json.load(f)
    assert got == expect

    with open(os.path.join(outdir, "Truth_files", "truth.mask.bed")) as f:
        mask_lines = [x.rstrip() for x in f]
    assert mask_lines == ["ref\t800\t801"]

    expect_recall_vcf = os.path.join(data_dir, "eval_one_fasta.expect_recall.vcf")
    got_recall_vcf = os.path.join(
        outdir, "Results", "Varifier", "recall", "recall.vcf.masked.vcf"
    )
    assert vcfs_same_exclude_headers(got_recall_vcf, expect_recall_vcf)
    expect_precision_vcf = os.path.join(data_dir, "eval_one_fasta.expect_precision.vcf")
    got_precision_vcf = os.path.join(outdir, "Results", "Varifier", "precision.vcf")
    assert vcfs_same_exclude_headers(got_precision_vcf, expect_precision_vcf)
    utils.syscall(f"rm -r {outdir}*")

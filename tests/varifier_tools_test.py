import os
import pytest

from cte import varifier_tools

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "varifier_tools")


def test_varifier_vcf_to_counts_and_errors():
    vcf_file = os.path.join(data_dir, "varifier_vcf_to_counts_and_errors.vcf")
    # primers should be a dict with keys = positions, and values = primer info,
    # but the funtion we're testing only uses the keys. We can just make a set
    # of numbers instead.
    primers = {4, 5, 6, 7, 14, 15, 16}
    first_primer_start = min(primers)
    last_primer_end = max(primers)
    got_counts, got_errors = varifier_tools.varifier_vcf_to_counts_and_errors(vcf_file, primers, first_primer_start, last_primer_end, "FP")
    expect_counts = {
        "snp": {"TP": 3, "FP": 2},
        "indel": {"TP": 2, "FP": 2},
        "primer_snp": {"TP": 1, "FP": 1},
        "primer_indel": {"TP": 1, "FP": 1},
    }
    assert got_counts == expect_counts
    expect_errors = [
        "ref\t5\t0\tG\tTT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
        "ref\t7\t0\tG\tT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
        "ref\t10\t0\tG\tTT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
        "ref\t11\t0\tG\tT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
    ]
    assert got_errors == expect_errors

    got_counts, got_errors = varifier_tools.varifier_vcf_to_counts_and_errors(vcf_file, primers, first_primer_start, last_primer_end, "FN")
    expect_counts = {
        "snp": {"TP": 3, "FN": 2},
        "indel": {"TP": 2, "FN": 2},
        "primer_snp": {"TP": 1, "FN": 1},
        "primer_indel": {"TP": 1, "FN": 1},
    }
    assert got_counts == expect_counts
    assert got_errors == expect_errors


def test_varifier_outdir_to_stats():
    indir = os.path.join(data_dir, "varifier_outdir_to_stats")
    primers = {4, 5, 6, 7, 14, 15, 16}
    first_primer_start = min(primers)
    last_primer_end = max(primers)
    got = varifier_tools.varifier_outdir_to_stats(indir, primers, first_primer_start, last_primer_end)
    expect = {
        "snp": {"TP": 3, "FP": 2, "FN": 2},
        "indel": {"TP": 3, "FP": 1, "FN": 1},
        "primer_snp": {"TP": 1, "FP": 1, "FN": 1},
        "primer_indel": {"TP": 1, "FP": 1, "FN": 1},
        "Errors": {
            "FP": [
                "ref\t5\t0\tG\tTT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
                "ref\t7\t0\tG\tT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
                "ref\t11\t0\tG\tA\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
            ],
            "FN": [
                "ref\t5\t0\tG\tTA\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
                "ref\t7\t0\tG\tC\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
                "ref\t11\t0\tG\tT\t.\tPASS\t.\tGT:VFR_IN_MASK:VFR_RESULT\t1/1:0:FP",
            ]
        },
    }
    assert got == expect

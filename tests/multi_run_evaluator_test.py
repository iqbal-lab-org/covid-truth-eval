import copy
import filecmp
import json
import os
import pytest
import random

from cte import multi_run_evaluator, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "multi_run_evaluator")


def test_load_manifest_tsv():
    infile = os.path.join(data_dir, "load_manifest_tsv.bad.tsv")
    with pytest.raises(Exception):
        multi_run_evaluator.load_manifest_tsv(infile)

    infile = os.path.join(data_dir, "load_manifest_tsv.tsv")
    got = multi_run_evaluator.load_manifest_tsv(infile)
    expect = {
        "run1": {
            "name": "run1",
            "eval_fasta": "run1.fa",
            "primers": "COVID-ARTIC-V3",
            "truth_vcf": "truth_a.vcf",
        },
        "run2": {
            "name": "run2",
            "eval_fasta": "run2.fa",
            "primers": "COVID-ARTIC-V4",
            "truth_vcf": "truth_a.vcf",
        },
        "run3": {
            "name": "run3",
            "eval_fasta": "run3.fa",
            "primers": "COVID-ARTIC-V4",
            "truth_vcf": "truth_b.vcf",
        },
    }
    assert got == expect


def test_evaluate_runs():
    # Make minimal test data for 3 runs. We're more interested in the overall
    # pipeline running on all runs, not the quality of results - that is tested
    # in one_run_evaluator_test. Spike in a couple of SNPs so we have something
    # in results other than all ref
    tmp_data_root = "tmp.evaluate_runs.data"
    utils.syscall(f"rm -rf {tmp_data_root}")
    os.mkdir(tmp_data_root)
    ref_fasta = os.path.join(tmp_data_root, "ref.fa")
    manifest_tsv = os.path.join(tmp_data_root, "manifest.tsv")

    random.seed(42)
    ref_seq = random.choices(["A", "C", "G", "T"], k=1000)
    ref_seq[500] = "A"
    ref_seq[600] = "T"
    with open(ref_fasta, "w") as f:
        print(">ref", file=f)
        print("".join(ref_seq), file=f)

    truth_vcfs = []
    primers = []
    eval_fastas = []

    truth_vcfs.append(os.path.join(tmp_data_root, "truth.a.vcf"))
    with open(truth_vcfs[0], "w") as f:
        print(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1",
            "ref\t501\t.\tA\tT\t42.42\tPASS\t.\tGT\t1/1",
            sep="\n",
            file=f,
        )
    truth_vcfs.append(truth_vcfs[0])

    truth_vcfs.append(os.path.join(tmp_data_root, "truth.b.vcf"))
    with open(truth_vcfs[2], "w") as f:
        print(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1",
            "ref\t501\t.\tA\tT\t42.42\tPASS\t.\tGT\t1/1",
            "ref\t601\t.\tT\tG\t42.42\tPASS\t.\tGT\t1/1",
            sep="\n",
            file=f,
        )

    eval_seq = copy.copy(ref_seq)
    eval_seq[500] = "T"
    eval_fastas.append(os.path.join(tmp_data_root, "truth.a.eval1.fa"))
    with open(eval_fastas[0], "w") as f:
        print(">eval1", "".join(eval_seq), sep="\n", file=f)

    eval_seq[500] = "G"
    eval_fastas.append(os.path.join(tmp_data_root, "truth.a.eval2.fa"))
    with open(eval_fastas[1], "w") as f:
        print(">eval2", "".join(eval_seq), sep="\n", file=f)

    eval_seq[600] = "G"
    eval_fastas.append(os.path.join(tmp_data_root, "truth.b.eval1.fa"))
    with open(eval_fastas[2], "w") as f:
        print(">eval3", "".join(eval_seq), sep="\n", file=f)

    primers_tsv = os.path.join(tmp_data_root, "primers.tsv")
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

    with open(manifest_tsv, "w") as f:
        print("name", "truth_vcf", "eval_fasta", "primers", sep="\t", file=f)
        for i in range(3):
            print(
                f"name{i}", truth_vcfs[i], eval_fastas[i], primers_tsv, sep="\t", file=f
            )

    outdir = "tmp.evaluate_runs"
    utils.syscall(f"rm -rf {outdir}")
    multi_run_evaluator.evaluate_runs(manifest_tsv, outdir, ref_fasta)
    with open(os.path.join(outdir, "results.json")) as f:
        got_results = json.load(f)
    with open(os.path.join(data_dir, "evaluate_runs.expect_results.json")) as f:
        expect_results = json.load(f)
    assert got_results == expect_results

    expect_tsv = os.path.join(data_dir, "evaluate_runs.expect_results.tsv")
    got_tsv = os.path.join(outdir, "results.tsv")
    assert filecmp.cmp(got_tsv, expect_tsv, shallow=False)

    expect_tsv = os.path.join(data_dir, "evaluate_runs.expect_results_per_run.tsv")
    got_tsv = os.path.join(outdir, "results_per_run.tsv")
    assert filecmp.cmp(got_tsv, expect_tsv, shallow=False)

    utils.syscall(f"rm -r {outdir} {tmp_data_root}")

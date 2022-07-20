import copy
from unittest import mock
import os
import pytest
import random

from cte import amplicon_scheme, msa, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "msa_test")


def test_aln_bases_to_stats_row_and_col():
    f = msa.aln_bases_to_stats_row_and_col
    row = msa.StatRow
    col = msa.StatCol
    assert f("A", "Z", "A") == (row.Dropped_amplicon, col.Called_ref)
    assert f("A", "Z", "N") == (row.Dropped_amplicon, col.Called_N)
    assert f("A", "Z", "Z") == (row.Dropped_amplicon, col.Dropped_amplicon)
    assert f("A", "Z", "C") == (row.Dropped_amplicon, col.Called_other)
    assert f("A", "Z", "e") == (row.Dropped_amplicon, col.No_call_genome_ends)
    assert f("T", "T", "-") == (row.True_ref, col.Called_wrong_indel)
    assert f("G", "T", "-") == (row.SNP_true_alt, col.Called_wrong_indel)
    assert f("-", "A", "A") == (row.True_indel, col.Called_correct_alt)
    assert f("-", "A", "T") == (row.True_indel, col.Called_wrong_indel)
    assert f("A", "-", "A") == (row.True_indel, col.Called_wrong_indel)
    assert f("A", "-", "G") == (row.True_indel, col.Called_wrong_indel)
    assert f("-", "-", "A") == (row.True_ref, col.Called_wrong_indel)
    assert f("-", "A", "-") == (row.True_indel, col.Called_wrong_indel)
    assert f("A", "-", "-") == (row.True_indel, col.Called_correct_alt)


def test_msa():
    # Basic test of Msa class. We test in more detail elsewhere with end-to-end
    # tests of the whole pipeline
    random.seed(42)
    root_out = "tmp.msa"
    utils.syscall(f"rm -rf {root_out}")
    os.mkdir(root_out)
    ref_fasta = os.path.join(root_out, "ref.fa")
    truth_fasta = os.path.join(root_out, "truth.fa")
    eval_fasta = os.path.join(root_out, "eval.fa")
    ref_seq = random.choices(["A", "C", "G", "T"], k=50)
    truth_seq = copy.copy(ref_seq)
    eval_seq = copy.copy(ref_seq)
    # add a few variants. Throw in IUPAC codes to check doesn't crash aligner
    ref_seq[14] = "A"
    ref_seq[20] = "A"
    original_ref_seq = "".join(ref_seq)
    truth_seq[20] = "Y"
    eval_seq[20] = "N"
    eval_seq[26:35] = "N" * 9
    truth_seq.insert(46, "A")
    truth_seq.pop(3)
    with open(ref_fasta, "w") as f:
        print(">ref", "".join(ref_seq), sep="\n", file=f)
    with open(truth_fasta, "w") as f:
        print(">truth", "".join(truth_seq), sep="\n", file=f)
    with open(eval_fasta, "w") as f:
        print(">eval", "".join(eval_seq), sep="\n", file=f)

    outdir = os.path.join(root_out, "msa")
    seqaln = msa.Msa(ref_fasta, truth_fasta, eval_fasta)
    seqaln.run_msa(outdir)
    truth_seq.insert(3, "-")
    ref_seq.insert(44, "-")
    eval_seq.insert(44, "-")
    assert seqaln.ref_aln == ref_seq
    assert seqaln.truth_aln == truth_seq
    assert seqaln.eval_aln == eval_seq

    seqaln.make_coords_lookups()
    expect_coords = list(range(44)) + list(range(45, 51))
    assert seqaln.ref_coords == expect_coords
    assert seqaln.eval_coords == expect_coords

    vcf_file = os.path.join(root_out, "truth_vars.vcf")
    with open(vcf_file, "w") as f:
        print(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "sample",
            sep="\t",
            file=f,
        )
        print(
            "ref",
            15,
            ".",
            ref_seq[14],
            "G",
            ".",
            "PASS",
            ".",
            "GT",
            "1/1",
            sep="\t",
            file=f,
        )
        print(
            "ref",
            27,
            ".",
            ref_seq[25],
            "N",
            ".",
            "DROPPED_AMP",
            "AMP_START=26;AMP_END=35",
            "GT",
            "1/1",
            sep="\t",
            file=f,
        )
        print(
            "ref",
            33,
            ".",
            ref_seq[31],
            "N",
            ".",
            "DROPPED_AMP",
            "AMP_START=32;AMP_END=40",
            "GT",
            "1/1",
            sep="\t",
            file=f,
        )
    seqaln.add_truth_dropped_amps(vcf_file)
    truth_seq[26:41] = "Z" * 15
    assert seqaln.truth_aln == truth_seq

    primers_file = os.path.join(root_out, "primers.tsv")
    with open(primers_file, "w") as f:
        print(
            "Amplicon_name",
            "Primer_name",
            "Left_or_right",
            "Sequence",
            "Position",
            sep="\t",
            file=f,
        )
        print("amp1", "p1_l", "left", original_ref_seq[10:13], 10, sep="\t", file=f)
        print("amp1", "p1_r", "right", original_ref_seq[27:30], 27, sep="\t", file=f)
        print("amp2", "p2_l", "left", original_ref_seq[26:29], 26, sep="\t", file=f)
        print("amp2", "p2_r", "right", original_ref_seq[33:36], 33, sep="\t", file=f)
        print("amp3", "p3_l", "left", original_ref_seq[32:35], 32, sep="\t", file=f)
        print("amp3", "p3_r", "right", original_ref_seq[38:41], 38, sep="\t", file=f)
        print("amp4", "p4_l", "left", original_ref_seq[40:42], 40, sep="\t", file=f)
        print("amp4", "p4_r", "right", original_ref_seq[47:49], 47, sep="\t", file=f)
    amp_scheme = amplicon_scheme.AmpliconScheme(primers_file)

    seqaln.add_eval_dropped_amps(amp_scheme, min_gap_length=5)
    eval_seq[26:35] = "Z" * 9
    assert seqaln.eval_aln == eval_seq
    utils.syscall(f"rm -r {root_out}")


def test_stats_from_three_way_aln():
    # Do a few basic tests here. We test the end-to-end one_run_evaluator
    # properly on a lot more test cases. We'd just be effectively repeating
    # them here
    amp_scheme = mock.Mock()
    primer_positions = {1, 17}
    amp_scheme.in_primer = lambda x: x in primer_positions
    # in the real scheme, these next two lines would match the primer
    # positions. Here we're just hacking by forcing the coords, to test
    # that we only calculate stats inside this range
    amp_scheme.first_primer_start = 0
    amp_scheme.last_primer_end = 3
    assert amp_scheme.in_primer(1)
    assert not amp_scheme.in_primer(2)

    ref = "ACGT"
    tru = "ACGT"
    con = "ACGT"
    stats = msa.make_empty_stats_dict()
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 4
    stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 1
    assert msa.stats_from_three_way_aln(ref, tru, con, amp_scheme) == stats
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 3
    amp_scheme.first_primer_start = 1
    assert msa.stats_from_three_way_aln(ref, tru, con, amp_scheme) == stats
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 2
    amp_scheme.last_primer_end = 2
    assert msa.stats_from_three_way_aln(ref, tru, con, amp_scheme) == stats

    # diffs marked with * so easier to spot!
    #       * * *  * *
    ref = "CATACTGA-AT"
    tru = "CATTCTGA-AA"
    con = "CTTYC-GACAT"
    stats = msa.make_empty_stats_dict()
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 6
    # SNP at position 1:
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_alt] = 1
    stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_alt] = 1
    # SNP at position 3:
    stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_wrong_IUPAC] = 1
    # deletion at position 5, and insertion after position 7:
    stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] = 2
    # SNP at the end:
    stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_ref] = 1
    tsv_out = "tmp.msa.stats_from_three_way_aln.tsv"
    utils.syscall(f"rm -f {tsv_out}")
    amp_scheme.first_primer_start = 0
    amp_scheme.last_primer_end = 100
    assert (
        msa.stats_from_three_way_aln(ref, tru, con, amp_scheme, tsv_out=tsv_out)
        == stats
    )
    assert os.path.exists(tsv_out)
    utils.syscall(f"rm {tsv_out}")

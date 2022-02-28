import filecmp
import os
import pytest

from cte import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_load_single_seq_fasta():
    tmp_fasta = "tmp.test_load_single_seq_fasta.fa"
    name = "test_seq"
    seq = "ACGT"
    with open(tmp_fasta, "w") as f:
        print(f">{name}", file=f)
        print(seq, file=f)

    got = utils.load_single_seq_fasta(tmp_fasta)
    assert got.id == name
    assert got.seq == seq
    os.unlink(tmp_fasta)


def test_vcf_file_to_sorted_list():
    vcf_file = os.path.join(data_dir, "vcf_file_to_sorted_list.vcf")
    got = utils.vcf_file_to_sorted_list(vcf_file)
    assert len(got) == 2
    assert got[0].POS < got[1].POS

    bad_vcf = os.path.join(data_dir, "vcf_file_to_sorted_list.bad_two_refs.vcf")
    with pytest.raises(Exception):
        got = utils.vcf_file_to_sorted_list(bad_vcf)


def test_apply_variants_to_genome():
    ref_fasta = os.path.join(data_dir, "apply_variants_to_genome.ref.fa")
    vcf_file = os.path.join(data_dir, "apply_variants_to_genome.vcf")
    expect_fasta = os.path.join(data_dir, "apply_variants_to_genome.expect.fa")
    out_fasta = "tmp.test_apply_variants_to_genome.fa"
    utils.syscall(f"rm -f {out_fasta}")
    utils.apply_variants_to_genome(vcf_file, out_fasta, ref_fasta=ref_fasta)
    assert filecmp.cmp(expect_fasta, out_fasta, shallow=False)
    os.unlink(out_fasta)

    vcf_file = os.path.join(data_dir, "apply_variants_to_genome.bad.vcf")
    with pytest.raises(Exception):
        utils.apply_variants_to_genome(vcf_file, out_fasta, ref_fasta=ref_fasta)


def test_string_to_gaps():
    assert utils.string_to_gaps("A", 1) == []
    assert utils.string_to_gaps("AN", 1) == [(1, 1)]
    assert utils.string_to_gaps("AN", 2) == []
    assert utils.string_to_gaps("ANN", 2) == [(1, 2)]
    assert utils.string_to_gaps("ANN", 3) == []
    assert utils.string_to_gaps("ANA", 1) == [(1, 1)]
    assert utils.string_to_gaps("NANA", 1) == [(0, 0), (2, 2)]
    assert utils.string_to_gaps("NAANNNNANN", 1) == [(0, 0), (3, 6), (8, 9)]
    assert utils.string_to_gaps("NAANNNNANN", 2) == [(3, 6), (8, 9)]
    assert utils.string_to_gaps("NAANNNNANN", 3) == [(3, 6)]
    assert utils.string_to_gaps("NAANNNNANN", 4) == [(3, 6)]
    assert utils.string_to_gaps("NAANNNNANN", 5) == []

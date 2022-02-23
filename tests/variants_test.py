import os
import pytest

import cluster_vcf_records as cvr

from cte import variants

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "variants")


#def test_qry_to_ref_coords_from_msa():
#    #            01234 567
#    ref_aln = "--ACGTC-AGT"
#    qry_aln = "CGAC-TNTAGT"
#    #          0123 456789
#    expect = [0, 0, 0, 1, 3, 4, 5, 5, 6, 7]
#    assert variants.qry_to_ref_coords_from_msa(ref_aln, qry_aln) == expect


def test_variant_init():
    primers = {0: "primer1", 1: "primer1", 100: "primer2"}
    record = cvr.vcf_record.VcfRecord("ref\t3\t.\tA\tC\t.\tPASS\t.\tGT\t1/1")
    variant = variants.Variant(record, primers, True)
    assert variant.is_snp
    assert not variant.in_primer
    assert not variant.is_het
    assert not variant.is_unsure


    primers[2] = "primer1"
    variant = variants.Variant(record, primers, False)
    assert variant.is_snp
    assert variant.in_primer
    assert not variant.is_het
    assert not variant.is_unsure

    record = cvr.vcf_record.VcfRecord("ref\t3\t.\tA\tC\t.\tHET\t.\tGT\t0/1")
    variant = variants.Variant(record, primers, True)
    assert variant.is_snp
    assert variant.in_primer
    assert variant.is_het
    assert not variant.is_unsure

    # TODO - test more cases



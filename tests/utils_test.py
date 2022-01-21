import os
import pytest

from cte import varifier_tools

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_apply_variants_to_genome():
    #utils.apply_variants_to_genome(ref_fasta, vcf_file, out_fasta)
    # FIXME
    pass

import os
import pytest

from cte import varifier_tools

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "varifier_tools")


def test_varfifier_vcf_to_counts():
    # FIXME
    # varifier_tools.varfifier_vcf_to_counts(infile, primers, first_primer_start, last_primer_end, fp_or_fn)
    pass


def test_varifier_outdir_to_stats():
    # FIXME
    # varifier_tools.varifier_outdir_to_stats(dirname, primers, first_primer_start, last_primer_end)
    pass

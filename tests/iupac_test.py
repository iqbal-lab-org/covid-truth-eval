import pytest

from cte import iupac

def test_expand_nucleotide_list():
    assert {"A"} == iupac.expand_nucleotide_list(["A"])
    assert {"A"} == iupac.expand_nucleotide_list(["A", "A"])
    assert {"A", "G"} == iupac.expand_nucleotide_list(["R"])
    assert {"A", "G"} == iupac.expand_nucleotide_list(["A", "R"])
    assert {"A", "C", "G", "T"} == iupac.expand_nucleotide_list(["R", "Y"])


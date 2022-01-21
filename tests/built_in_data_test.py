import os
import pytest

from cte import built_in_data


def test_covid_files_exist():
    assert os.path.exists(built_in_data.COVID_REF)
    for scheme_name in built_in_data.COVID_SCHEME_NAMES:
        assert os.path.exists(built_in_data.COVID_PRIMER_TSVS[scheme_name])

import json
import os

this_dir = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(this_dir, "data")

manifest_json = os.path.join(DATA_DIR, "manifest.json")
with open(manifest_json) as f:
    data_files = json.load(f)

COVID_REF = os.path.join(DATA_DIR, data_files["covid.MN908947"]["ref_fasta"])
COVID_PRIMER_TSVS = {
    k: os.path.join(DATA_DIR, v)
    for k, v in data_files["covid.MN908947"]["primer_tsvs"].items()
}
COVID_SCHEME_NAMES = sorted(list(COVID_PRIMER_TSVS.keys()))

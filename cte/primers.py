import csv
import os

from cte import built_in_data


def load_viridian_workflow_primers_tsv(infile):
    """Loads TSV primers file in 'viridian workflow' format.
    Returns tuple:
      - a lookup of genome position (zero based) to set of (amplicon, primer)
        tuples that position is in
      - zero-based coord of start of first primer
      - zero-based coord of end of last primer"""
    primer_lookup = {}
    first_primer_start = float("inf")
    last_primer_end = 0
    if not os.path.exists(infile):
        if infile in built_in_data.COVID_PRIMER_TSVS:
            infile = built_in_data.COVID_PRIMER_TSVS[infile]
        else:
            raise Exception(f"Primers '{infile}' not found")

    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            start = int(d["Position"])
            end = start + len(d["Sequence"]) - 1
            first_primer_start = min(first_primer_start, start)
            last_primer_end = max(last_primer_end, end)

            for i in range(start, end + 1):
                if i not in primer_lookup:
                    primer_lookup[i] = set()
                primer_lookup[i].add((d["Amplicon_name"], d["Primer_name"]))

    return primer_lookup, first_primer_start, last_primer_end

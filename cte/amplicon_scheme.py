import csv
import logging
import os

from intervaltree import IntervalTree

from cte import built_in_data


class AmpliconScheme:
    def __init__(self, infile):
        """Loads data from TSV primers file in 'viridian workflow' format"""
        self.primers = IntervalTree()
        first_primer_start = float("inf")
        last_primer_end = 0
        amplicons = {}
        if not os.path.exists(infile):
            if infile in built_in_data.COVID_PRIMER_TSVS:
                infile = built_in_data.COVID_PRIMER_TSVS[infile]
            else:
                raise Exception(f"Primers '{infile}' not found")

        logging.info(f"Loading primers from file {infile}")
        with open(infile) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for d in reader:
                start = int(d["Position"])
                end = start + len(d["Sequence"]) - 1
                if d["Amplicon_name"] not in amplicons:
                    amplicons[d["Amplicon_name"]] = {"start": float("inf"), "end": -1}
                amplicons[d["Amplicon_name"]]["start"] = min(amplicons[d["Amplicon_name"]]["start"], start)
                amplicons[d["Amplicon_name"]]["end"] = max(amplicons[d["Amplicon_name"]]["end"], end)
                first_primer_start = min(first_primer_start, start)
                last_primer_end = max(last_primer_end, end)
                self.primers[start:end+1] = d["Primer_name"]

        self.amplicons = IntervalTree()
        for name, d in amplicons.items():
            self.amplicons[d["start"]:d["end"] + 1] = name

        self.first_primer_start = min(self.amplicons).begin
        self.last_primer_end = max(self.amplicons).end - 1


    def overlaps_primer(self, start, end):
        return len(self.primers[start:end+1]) > 0


    def overlapping_amplicons(self, start, end):
        matches = self.amplicons[start:end+1]
        if len(matches) > 0:
            return [{"name": x.data, "start": x.begin, "end":x.end - 1} for x in sorted(matches)]
        else:
            return None


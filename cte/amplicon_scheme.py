import csv
import logging
import os

from intervaltree import IntervalTree

from cte import built_in_data


class AmpliconScheme:
    def __init__(self, infile):
        """Loads data from TSV primers file in 'viridian workflow' format"""
        if not os.path.exists(infile):
            if infile in built_in_data.COVID_PRIMER_TSVS:
                infile = built_in_data.COVID_PRIMER_TSVS[infile]
            else:
                raise Exception(f"Primers '{infile}' not found")

        self.primers, self.amplicons = self._load_primers_file(infile)
        self.amplicon_tree = IntervalTree()
        for name, d in self.amplicons.items():
            self.amplicon_tree[d["start"] : d["end"] + 1] = name

        self.first_primer_start = min(self.amplicon_tree).begin
        self.last_primer_end = max(self.amplicon_tree).end - 1

    @classmethod
    def _load_primers_file(cls, infile):
        logging.info(f"Loading primers from file {infile}")
        primers = IntervalTree()
        amplicons = {}
        with open(infile) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for d in reader:
                start = int(d["Position"])
                end = start + len(d["Sequence"]) - 1
                name = d["Amplicon_name"]
                if name not in amplicons:
                    amplicons[name] = {
                        "start": float("inf"),
                        "end": -1,
                        "left_primers": [],
                        "right_primers": [],
                    }
                amplicons[name]["start"] = min(amplicons[name]["start"], start)
                amplicons[name]["end"] = max(amplicons[name]["end"], end)
                l_or_r = d["Left_or_right"].lower()
                amplicons[name][f"{l_or_r}_primers"].append((start, end))
                primers[start : end + 1] = d["Primer_name"]

        amps_in_order = sorted(amplicons.keys(), key=lambda k: amplicons[k]["start"])
        for i, amp in enumerate(amps_in_order):
            if i == 0:
                amplicons[amp]["previous"] = None
            else:
                amplicons[amp]["previous"] = amps_in_order[i - 1]

            if i < len(amps_in_order) - 1:
                amplicons[amp]["next"] = amps_in_order[i + 1]
            else:
                amplicons[amp]["next"] = None

        return primers, amplicons

    def in_primer(self, position):
        return len(self.primers[position : position + 1]) > 0

    def overlaps_primer(self, start, end):
        return len(self.primers[start : end + 1]) > 0

    def overlapping_amplicons(self, start, end):
        matches = self.amplicon_tree[start : end + 1]
        if len(matches) > 0:
            return [x.data for x in sorted(matches)]
        else:
            return None

    def amp_coords(self, amp_name):
        return self.amplicons[amp_name]["start"], self.amplicons[amp_name]["end"]

    def inner_amp_coords(self, amp_name):
        """Returns the "inner coords", which is end of previous amplicon to
        start of next amplicon. This is the minimum we would expect to get
        masked if a consensus has a dropped amplicon"""
        amp_dict = self.amplicons[amp_name]
        if amp_dict["previous"] is None:
            start = max([x[1] for x in amp_dict["left_primers"]]) + 1
        else:
            start = self.amplicons[amp_dict["previous"]]["end"] + 1
        if amp_dict["next"] is None:
            end = min([x[0] for x in amp_dict["right_primers"]]) - 1
        else:
            end = self.amplicons[amp_dict["next"]]["start"] - 1
        assert start < end
        return start, end

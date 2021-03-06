ambiguous_codes = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

rev_ambiguous_codes = {"".join(sorted(list(v))): k for k, v in ambiguous_codes.items()}


def expand_nucleotide_list(l):
    expanded = set()
    for x in l:
        expanded.update(ambiguous_codes.get(x, {x}))
    return expanded

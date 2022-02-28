import copy
from enum import Enum, auto
import json
import os

import cluster_vcf_records
import pyfastaq

from cte import iupac, utils


class StatRow(Enum):
    True_ref = auto()
    SNP_true_alt = auto()
    SNP_true_mixed = auto()
    Unknown_truth = auto()
    True_indel = auto()
    Dropped_amplicon = auto()


class StatCol(Enum):
    Called_ref = auto()
    Called_correct_alt = auto()
    Called_wrong_alt = auto()
    Called_N = auto()
    Called_correct_IUPAC = auto()
    Called_wrong_IUPAC = auto()
    Called_wrong_indel = auto()
    Dropped_amplicon = auto()
    Called_other = auto()


ACGT = {"A", "C", "G", "T"}
AMBIG_NON_N = {x for x in iupac.ambiguous_codes if x != "N"}


def make_empty_stats_dict():
    stats = {row: {} for row in StatRow}
    for row in StatRow:
        for col in StatCol:
            stats[row][col] = 0
    return {"All": copy.deepcopy(stats), "Primer": copy.deepcopy(stats)}


def write_stats_summary_tsv(stats, outfile):
    with open(outfile, "w") as f:
        print("Truth", *[x.name for x in StatCol], sep="\t", file=f)
        for row in StatRow:
            print(row.name, *[stats["All"][row][x] for x in StatCol], sep="\t", file=f)
        for row in StatRow:
            print(
                f"P_{row}",
                *[stats["Primer"][row][x] for x in StatCol],
                sep="\t",
                file=f,
            )


def write_stats_summary_json(stats, outfile):
    to_print = {"All": {}, "Primer": {}}
    for all_or_primer in to_print:
        for row, d in stats[all_or_primer].items():
            to_print[all_or_primer][row.name] = {k.name: v for k, v in d.items()}


    with open(outfile, "w") as f:
        json.dump(to_print, f, indent=2)


# Some notes on working out which rows and columns. There are a lot
# of cases to consider!
# Reminder - using "Z" to denote dropped amplicon.

# All indel possibilities
# (where the "x" could be any A/C/G/T/ambiguous IUPAC but not Z):
# ref  truth cons   row              col
# -    x     x      True_indel       Called_correct_alt
# -    x     Z      True_indel       Dropped_amplicon
# -    -     x      True_ref         Called_wrong_indel
# -    -     Z      True_ref         Dropped_amplicon
# -    x     -      True_indel       Called_wrong_indel
# x    -     x      True_indel       Called_wrong_indel
# x    -     Z      True_indel       Dropped_amplicon
# x    -     -      True_indel       Called_correct_alt
# x    x     -      True_ref         Called_wrong_indel
# x    Z     -      Dropped_amplicon Called_wrong_indel
# Should not happen:
# ref  truth cons
# -    Z     (anything)


# Non-indel possibilities
# ref  truth cons   row               col
# C    C     C      True_ref          Called_ref
# C    C     G      True_ref          Called_wrong_alt
# C    C     Y      True_ref          Called_wrong_IUPAC
# C    C     N      True_ref          Called_N
# C    C     Z      True_ref          Dropped_amplicon
# C    T     C      SNP_true_alt      Called_ref
# C    T     T      SNP_true_alt      Called_correct_alt
# C    T     G      SNP_true_alt      Called_wrong_alt
# ...etc ... end then dropped amlicons:
# C    Z     Z      Dropped_amplicon  Dropped_amplicon
# C    Z     C      Dropped_amplicon  Called_ref
# C    Z     T      Dropped_amplicon  Called_wrong_alt
# ...etc ...
# C    C     Z      True_ref          Dropped_amplicon
# C    T     Z      SNP_true_alt      Dropped_amplicon
# C    N     C      Unknown_truth     Called_ref
# ...etc
#
# The point here is that the row only depends on ref/truth.
# If ref==truth, it's True_ref, otherwise only depends on truth, so we can
# build a dictionary of truth to row:
TRUTH_TO_ROW = {x: StatRow.SNP_true_alt for x in ACGT}
for x in AMBIG_NON_N:
    TRUTH_TO_ROW[x] = StatRow.SNP_true_mixed
TRUTH_TO_ROW["N"] = StatRow.Unknown_truth
TRUTH_TO_ROW["Z"] = StatRow.Dropped_amplicon

# When ref==truth, there should be no variant. And then cons!=ref/truth, then
# column depends only on the consensus, so can build a dictionary.
CONS_FP_TO_COL = {x: StatCol.Called_wrong_alt for x in ACGT}
for x in AMBIG_NON_N:
    CONS_FP_TO_COL[x] = StatCol.Called_wrong_IUPAC
CONS_FP_TO_COL["N"] = StatCol.Called_N
CONS_FP_TO_COL["Z"] = StatCol.Dropped_amplicon


def aln_bases_to_stats_row_and_col(ref, truth, cons):
    assert ref == "-" or ref in ACGT
    if ref == "-" or truth == "-" or cons == "-":
        if ref == truth:
            row = StatRow.True_ref
        elif truth == "Z":
            row = StatRow.Dropped_amplicon
        else:
            row = StatRow.True_indel
        if cons == truth:
            col = StatCol.Called_correct_alt
        elif cons == "Z":
            col = StatCol.Dropped_amplicon
        else:
            col = StatCol.Called_wrong_indel
    else:
        if ref == truth:
            row = StatRow.True_ref
            col = StatCol.Called_ref if cons == ref else CONS_FP_TO_COL[cons]
        else:
            # if we're here then ref != truth, so there really is a variant.
            # Means we need to consider all three ref, truth, cons...
            row = TRUTH_TO_ROW[truth]

            if cons == "Z":
                col = StatCol.Dropped_amplicon
            elif cons == ref:
                col = StatCol.Called_ref
            elif cons == "N":
                col = StatCol.Called_N
            elif row == StatRow.Unknown_truth:
                col = StatCol.Called_other
            elif row == StatRow.SNP_true_alt:
                if cons == truth:
                    col = StatCol.Called_correct_alt
                elif cons in ACGT:
                    col = StatCol.Called_wrong_alt
                else:
                    assert cons in AMBIG_NON_N
                    col = StatCol.Called_wrong_IUPAC
            elif row == StatRow.SNP_true_mixed:
                if cons == truth:
                    col = StatCol.Called_correct_IUPAC
                elif cons in ACGT:
                    col = StatCol.Called_wrong_alt
                else:
                    assert cons in AMBIG_NON_N
                    col = StatCol.Called_wrong_IUPAC
            else:
                raise NotImplementedError(
                    f"Case not found for ref,truth,cons bases: {ref},{truth},{cons}"
                )

    return row, col


def stats_from_three_way_aln(ref_aln, truth_aln, eval_aln, amp_scheme, tsv_out=None):
    ref_pos = 0
    stats = make_empty_stats_dict()
    tsv_lines_out = [
        ("Ref_pos", "Ref", "Truth", "Consensus", "Truth_category", "Consensus_category")
    ]

    for aln_index, (ref_base, truth_base, eval_base) in enumerate(
        zip(ref_aln, truth_aln, eval_aln)
    ):
        if ref_pos < amp_scheme.first_primer_start:
            if ref_base != "-":
                ref_pos += 1
            continue
        elif ref_pos > amp_scheme.last_primer_end:
            break

        row, col = aln_bases_to_stats_row_and_col(ref_base, truth_base, eval_base)
        tsv_lines_out.append(
            (ref_pos + 1, ref_base, truth_base, eval_base, row.name, col.name)
        )
        stats["All"][row][col] += 1
        if amp_scheme.in_primer(ref_pos):
            stats["Primer"][row][col] += 1
        if ref_base != "-":
            ref_pos += 1

    if tsv_out is not None:
        with open(tsv_out, "w") as f:
            for t in tsv_lines_out:
                print(*t, sep="\t", file=f)

    return stats


class Msa:
    def __init__(self, ref_fasta, truth_fasta, eval_fasta):
        assert os.path.exists(ref_fasta)
        assert os.path.exists(truth_fasta)
        assert os.path.exists(eval_fasta)
        self.ref_fasta = ref_fasta
        self.truth_fasta = truth_fasta
        self.eval_fasta = eval_fasta
        self.dropped_amplicons = None
        self.ref_aln = None
        self.truth_aln = None
        self.eval_aln = None
        self.ref_coords = None
        self.eval_coords = None

    def run_msa(self, outdir):
        os.mkdir(outdir)
        all_seqs_fasta = os.path.join(outdir, "to_align.fa")
        with open(all_seqs_fasta, "w") as f:
            for fname in self.ref_fasta, self.truth_fasta, self.eval_fasta:
                reader = pyfastaq.sequences.file_reader(fname)
                for seq in reader:
                    print(seq, file=f)

        msa_fasta = os.path.join(outdir, "msa.fasta")
        command = f"mafft --anysymbol --nuc --quiet --maxiterate 1000 --thread 1 {all_seqs_fasta} > {msa_fasta}"
        utils.syscall(command)
        reader = pyfastaq.sequences.file_reader(msa_fasta)
        seqs = [list(x.seq.upper()) for x in reader]
        assert len(seqs) == 3
        assert len({len(x) for x in seqs}) == 1
        self.ref_aln, self.truth_aln, self.eval_aln = seqs

    def make_primer_mask(self, amp_scheme):
        assert self.ref_coords is not None
        raise NotImplementedError()

    def make_coords_lookups(self):
        # make lookup of original ref coord -> coord in msa
        self.ref_coords = [i for i, base in enumerate(self.ref_aln) if base != "-"]
        assert self.ref_coords[-1] + 1 == len(self.ref_aln)

        # lookup of original eval coord -> coord in msa
        self.eval_coords = [i for i, base in enumerate(self.eval_aln) if base != "-"]
        assert self.eval_coords[-1] + 1 == len(self.eval_aln)

    def add_truth_dropped_amps(self, truth_vcf_file):
        assert self.ref_coords is not None
        header, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            truth_vcf_file
        )
        for record in vcf_records:
            if "DROPPED_AMP" not in record.FILTER:
                continue

            assert "AMP_START" in record.INFO and "AMP_END" in record.INFO
            end = int(record.INFO["AMP_END"])
            assert int(record.INFO["AMP_START"]) == record.POS
            # for i in range(record.POS, end + 1):
            #    self.ref_aln[self.ref_coords[i]] = "D"
            for i in range(self.ref_coords[record.POS], self.ref_coords[end + 1]):
                self.truth_aln[i] = "Z"

    def add_eval_dropped_amps(self, amp_scheme, min_gap_length=50):
        eval_seq = "".join([x for x in self.eval_aln if x != "-"])
        gaps = utils.string_to_gaps(eval_seq, min_gap_length=min_gap_length)
        dropped_eval_ranges = []

        for (gap_start, gap_end) in gaps:
            amp_matches = amp_scheme.overlapping_amplicons(gap_start, gap_end)
            if amp_matches is None:
                continue

            for amp_name in amp_matches:
                inner_start, inner_end = amp_scheme.inner_amp_coords(amp_name)
                if gap_start <= inner_start <= inner_end <= gap_end:
                    amp_start, amp_end = amp_scheme.amp_coords(amp_name)
                    drop_start = max(amp_start, gap_start)
                    drop_end = min(amp_end, gap_end)
                    dropped_eval_ranges.append((drop_start, drop_end + 1))

        for (start, end) in dropped_eval_ranges:
            for i in range(start, end):
                self.eval_aln[self.eval_coords[i]] = "Z"

    def gather_stats(self, amp_scheme, per_pos_tsv=None):
        self.stats = stats_from_three_way_aln(
            self.ref_aln,
            self.truth_aln,
            self.eval_aln,
            amp_scheme,
            tsv_out=per_pos_tsv,
        )

    def write_stats_summary_tsv(self, outfile):
        write_stats_summary_tsv(self.stats, outfile)


    def write_stats_summary_json(self, outfile):
        write_stats_summary_json(self.stats, outfile)



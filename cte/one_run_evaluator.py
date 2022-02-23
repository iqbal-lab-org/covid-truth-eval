import copy
import json
import logging
import os

import cluster_vcf_records

from cte import amplicon_scheme, utils, variants, varifier_tools

STATS_ROWS = ["true_ref", "SNP_true_alt", "SNP_true_mixed", "unknown_truth", "True_indel", "Dropped_amplicon"]
STATS_COLS = ["Called_ref", "Called_correct_alt", "Called_wrong_alt", "Called_N", "Called_correct_IUPAC", "Called_wrong_IUPAC", "Called_wrong_indel", "Dropped_amplicon"]


def make_empty_stats_dict():
    stats = {row: {} for row in STATS_ROWS}
    for row in STATS_ROWS:
        for col in STATS_COLS:
            stats[row][col] = 0
    return {"All": copy.deepcopy(stats), "Primer": copy.deepcopy(stats)}


def write_stats_counts(stats, outfile):
    with open(outfile, "w") as f:
        print("Truth", *STATS_COLS, sep="\t", file=f)
        for row in STATS_ROWS:
            print(row, *[stats["All"][row][x] for x in STATS_COLS], sep="\t", file=f)
        for row in STATS_ROWS:
            print(f"P_{row}", *[stats["Primer"][row][x] for x in STATS_COLS], sep="\t", file=f)


class OneSampleEvaluator:
    def __init__(self, ref_fasta, truth_vcf, eval_fasta, outdir, primers_tsv, debug=False):
        self.ref_fasta = ref_fasta
        self.truth_vcf = truth_vcf
        self.eval_fasta = eval_fasta
        self.outdir = outdir
        self.debug = debug
        self.primers_tsv = primers_tsv
        self.stats = None
        self.inferred_vcf_dir = os.path.join(outdir, "Inferred_VCF_and_MSA")
        self.stats_tsv_out = os.path.join(self.outdir, "results.tsv")
        self.stats_json_out = os.path.join(self.outdir, "results.json")
        self.per_record_tsv_out = os.path.join(self.outdir, "per_record_breakdown.tsv")
        self.ref_seq = utils.load_single_seq_fasta(self.ref_fasta)
        self.used_ref_positions = None
        self.per_record_output = []
        self.site_count = None


    def init_stats(self):
        self.stats = make_empty_stats_dict()


    def handle_truth_variants(self):
        # Find the FNs: for each truth variant, see if it is in the called
        # variants. We'll be finding called variants that match truth variants.
        # Need to remember which called variants we have used, so that they are
        # not counted again when we check them all later for TP/FPs.
        self.used_called_variants = set()
        self.used_ref_positions = set()
        self.dropped_amp_ref_positions = set()
        self.site_count = 1

        for truth_var in self.truth_variants.vars.values():
            self.used_ref_positions.update(range(truth_var.ref_start(), truth_var.ref_end() + 1))
            row_key = truth_var.truth_type_for_reporting()
            called_matches = self.called_variants.get_overlapping_variants(truth_var)
            if len(called_matches) == 0:
                col_key = "Called_ref"
                self.per_record_output.append((self.site_count, "Check_truth", truth_var.in_primer, row_key, col_key, "Truth", truth_var.vcf_record))
                self.stats["All"][row_key][col_key] += truth_var.ref_length()
                if truth_var.in_primer:
                    self.stats["Primer"][row_key][col_key] += truth_var.ref_length()
                self.site_count += 1
                continue

            for called_match in called_matches:
                stat_increment = 1
                self.used_ref_positions.update(range(called_match.ref_start(), called_match.ref_end() + 1))

                if row_key != "unknown_truth" and called_match.is_same_variant(truth_var, self.ref_seq):
                    col_key = "Called_correct_IUPAC" if truth_var.is_het else "Called_correct_alt"
                    self.used_called_variants.add(called_match)
                elif row_key == "Dropped_amplicon":
                    all_amp_matches = called_match.amplicon_dropout_matches()
                    if all_amp_matches is not None:
                        amp_matches = [x for x in all_amp_matches if x["start"] == int(truth_var.vcf_record.INFO["AMP_START"])]
                    else:
                        amp_matches = []
                    assert len(amp_matches) <= 1
                    if len(amp_matches) == 0:
                        raise NotImplementedError()
                    else:
                        col_key =


                else:
                    # We have exactly one variant at the same position, and it is wrong.
                    # Need to figure out in what way it is wrong
                    self.used_called_variants.add(called_match)
                    if row_key == "True_indel" or not called_match.is_snp:
                        col_key = "Called_wrong_indel"
                    elif called_match.alt_is_all_N:
                        stat_increment = called_match.ref_length()
                        col_key = "Called_N"
                    elif called_match.is_het:
                        col_key = "Called_wrong_IUPAC"
                    elif called_match.is_snp:
                        col_key = "Called_wrong_alt"
                    else:
                        raise NotImplementedError(f"error getting result from truth_var:\n{truth_var}\nand called variant:\n{called}")

                self.per_record_output.append((self.site_count, "Check_truth", truth_var.in_primer, row_key, col_key, "Truth", truth_var.vcf_record))
                if called_match is not None:
                    self.per_record_output.append((self.site_count, "Check_truth", truth_var.in_primer, row_key, col_key, "Called", called_match.vcf_record))
                self.site_count += 1

                self.stats["All"][row_key][col_key] += stat_increment
                if row_key != "Dropped_amplicon" and truth_var.in_primer:
                    self.stats["Primer"][row_key][col_key] += stat_increment


    def handle_called_variants(self):
        for called_var in self.called_variants.vars.values():
            # We already went through the truth variants and matched each one to
            # a called variant (if possible). We don't want to count those called
            # variants twice
            if called_var in self.used_called_variants:
                continue

            # We should not get any truth matches, because those would have been
            # matched to this called variant when we ran self.handle_truth_variants()
            truth_matches = self.truth_variants.get_overlapping_variants(called_var)
            if len(truth_matches) > 0:
                variants = "\n".join([str(x.vcf_record) for x in truth_matches])
                raise Exception(f"Error detecting overlapping variants, possible counting twice error\nCalled variant:{called_var.vcf_record}\nOverlapping truth variants:\n{variants}")

            self.used_ref_positions.update(range(called_var.ref_start(), called_var.ref_end() + 1))
            stat_increment = 1

            if called_var.alt_is_all_N:
                stat_increment = called_var.ref_length()
                col_key = "Called_N"
            elif not called_var.is_snp:
                col_key = "Called_wrong_indel"
            elif called_var.is_het:
                col_key = "Called_wrong_IUPAC"
            elif called_var.is_snp:
                col_key = "Called_wrong_alt"
            else:
                raise Exception(f"Error determining error type of false positive: {called_var.vcf_record}")

            row_key = "true_ref"
            self.per_record_output.append((self.site_count, "Check_called", called_var.in_primer, row_key, col_key, "Called", called_var.vcf_record))
            self.site_count += 1

            self.stats["All"][row_key][col_key] += stat_increment
            if called_var.in_primer:
                self.stats["Primer"][row_key][col_key] += stat_increment


    def update_ref_position_counts(self):
        r = range(self.amplicon_scheme.first_primer_start, self.amplicon_scheme.last_primer_end + 1)
        unused_all = len([i for i in r if i not in self.used_ref_positions])
        unused_in_primer = len([i for i in r if i not in self.used_ref_positions and self.amplicon_scheme.overlaps_primer(i, i)])
        self.stats["All"]["true_ref"]["Called_ref"] += unused_all
        self.stats["Primer"]["true_ref"]["Called_ref"] += unused_in_primer


    def run(self):
        assert not os.path.exists(self.outdir)
        self.amplicon_scheme = amplicon_scheme.AmpliconScheme(self.primers_tsv)
        os.mkdir(self.outdir)
        self.ref_msa, self.qry_msa, called_vcf = varifier_tools.make_msa_and_vcf(
            self.ref_fasta,
            self.eval_fasta,
            self.amplicon_scheme.first_primer_start + 1,
            self.amplicon_scheme.last_primer_end + 1,
            self.inferred_vcf_dir,
            debug=self.debug
        )
        self.truth_variants = variants.VariantSet(self.truth_vcf, self.amplicon_scheme, True)
        self.called_variants = variants.VariantSet(called_vcf, self.amplicon_scheme, False)
        self.init_stats()
        self.per_record_output = []
        self.handle_truth_variants()
        self.handle_called_variants()
        self.update_ref_position_counts()
        write_stats_counts(self.stats, self.stats_tsv_out)
        with open(self.stats_json_out, "w") as f:
            json.dump(self.stats, f, indent=2)

        with open(self.per_record_tsv_out, "w") as f:
            site = 1
            print("Site", "Process", "In_primer", "Truth_category", "Called_category", "Variant_source",
                    "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample",
                    sep="\t", file=f)
            for fields in self.per_record_output:
                if fields[-1] is not None:
                    print(*fields, sep="\t", file=f)


def eval_one_fasta(outdir, fasta_to_eval, ref_fasta, truth_vcf, primers_name_or_tsv):
    logging.info(f"Start evaluating {fasta_to_eval}")
    logging.info(f"Using reference FASTA: {ref_fasta}")
    logging.info(f"Using truth VCF: {truth_vcf}")
    evaluator = OneSampleEvaluator(ref_fasta, truth_vcf, fasta_to_eval, outdir, primers_name_or_tsv, debug=True)
    evaluator.run()
    logging.info(f"Finished evaluating {fasta_to_eval}")
    return evaluator.stats

import logging
import os

import cluster_vcf_records
import json

from cte import primers, varifier_tools, utils


class OneSampleEvaluator:
    def __init__(self, ref_fasta, truth_vcf, truth_files_dir):
        self.ref_fasta = ref_fasta
        self.truth_vcf = truth_vcf
        self.truth_files_dir = truth_files_dir
        os.mkdir(self.truth_files_dir)
        self.filtered_truth_vcf = os.path.join(
            self.truth_files_dir, "truth.filtered.vcf"
        )
        self.truth_mask_bed = os.path.join(self.truth_files_dir, "truth.mask.bed")
        self.truth_fasta = os.path.join(self.truth_files_dir, "truth.fa")
        self._make_filtered_truth_vcf_and_truth_mask()
        utils.apply_variants_to_genome(
            self.filtered_truth_vcf, self.truth_fasta, ref_fasta=self.ref_fasta
        )

    def _make_filtered_truth_vcf_and_truth_mask(self):
        header, records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
            self.truth_vcf
        )
        mask_lines = []
        with open(self.filtered_truth_vcf, "w") as f:
            print(*header, sep="\n", file=f)
            for record in records:
                if "PASS" in record.FILTER:
                    print(record, file=f)
                elif "HET" in record.FILTER or "UNSURE" in record.FILTER:
                    mask_lines.append(
                        (record.CHROM, record.POS, record.ref_end_pos() + 1)
                    )

        with open(self.truth_mask_bed, "w") as f:
            for line in mask_lines:
                print(*line, sep="\t", file=f)

    def evaluate_one_fasta(self, outdir, fasta_to_eval, primers_tsv):
        assert not os.path.exists(outdir)
        os.mkdir(outdir)
        inferred_vcf_dir = os.path.join(outdir, "Inferred_from_consensus")
        varifier_dir = os.path.join(outdir, "Varifier")
        logging.info(
            f"Making VCF of variants inferred from FASTA being evaluated {fasta_to_eval}"
        )
        utils.syscall(
            f"varifier make_truth_vcf --global_align {fasta_to_eval} {self.ref_fasta} {inferred_vcf_dir}"
        )
        vcf_to_eval = os.path.join(inferred_vcf_dir, "04.truth.vcf")
        logging.info("Evaluating inferred variants")
        utils.syscall(
            f"varifier --debug vcf_eval --ref_mask {self.truth_mask_bed} --truth_vcf {self.filtered_truth_vcf} {self.truth_fasta} {self.ref_fasta} {vcf_to_eval} {varifier_dir}"
        )
        (
            primer_lookup,
            first_primer_start,
            last_primer_end,
        ) = primers.load_viridian_workflow_primers_tsv(primers_tsv)
        logging.info(
            f"Gathering summary stats from Varifier output directory {varifier_dir}"
        )
        varifier_stats = varifier_tools.varifier_outdir_to_stats(
            varifier_dir, primer_lookup, first_primer_start, last_primer_end
        )

        json_out = os.path.join(outdir, "results.json")
        with open(json_out, "w") as f:
            json.dump(varifier_stats, f, indent=2)
        return varifier_stats


def eval_one_fasta(outdir, fasta_to_eval, ref_fasta, truth_vcf, primers_name_or_tsv):
    logging.info(f"Start evaluating {fasta_to_eval}")
    logging.info(f"Using reference FASTA: {ref_fasta}")
    logging.info(f"Using truth VCF: {truth_vcf}")
    os.mkdir(outdir)
    truth_files_dir = os.path.join(outdir, "Truth_files")
    evaluator = OneSampleEvaluator(ref_fasta, truth_vcf, truth_files_dir)
    eval_dir = os.path.join(outdir, "Results")
    logging.info(f"Finished evaluating {fasta_to_eval}")
    return evaluator.evaluate_one_fasta(
        eval_dir,
        fasta_to_eval,
        primers_name_or_tsv,
    )

import copy
import csv
import json
import logging
import os

from cte import msa, one_run_evaluator


def load_manifest_tsv(infile):
    manifest_data = {}
    expect_cols = {"name", "truth_vcf", "eval_fasta", "primers"}

    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not expect_cols.issubset(set(reader.fieldnames)):
            raise Exception(
                f"manifest tsv file {infile} must contain these columns: {','.join(sorted(list(expect_cols)))}"
            )

        for d in reader:
            if d["name"] in manifest_data:
                raise Exception(
                    f"Name must be unique across all input. Got this more than once: {d['name']}"
                )
            manifest_data[d["name"]] = d

    return manifest_data


def evaluate_runs(manifest_tsv, outdir, ref_fasta, debug=False):
    os.mkdir(outdir)
    processing_dir = os.path.join(outdir, "Processing")
    os.mkdir(processing_dir)
    logging.info(f"Load runs from file {manifest_tsv}")
    manifest_data = load_manifest_tsv(manifest_tsv)
    per_run_results = {}
    results_totals = None
    completed = 0
    per_sample_tsv_lines = []

    for run_name, run_data in sorted(manifest_data.items()):
        logging.info("=" * 60)
        logging.info(
            f"Start processing run {run_name} {completed+1}/{len(manifest_data)}"
        )
        run_outdir = os.path.join(processing_dir, run_name)
        new_results = one_run_evaluator.eval_one_fasta(
            run_outdir,
            run_data["eval_fasta"],
            ref_fasta,
            run_data["truth_vcf"],
            run_data["primers"],
            debug=debug,
        )
        per_run_results[run_name] = msa.stats_to_json_friendly(new_results)

        tsv_lines = msa.stats_to_summary_tsv_lines(new_results, run_name)
        if len(per_sample_tsv_lines) == 0:
            per_sample_tsv_lines.extend(tsv_lines)
        else:
            per_sample_tsv_lines.extend(tsv_lines[1:])

        if results_totals is None:
            results_totals = copy.deepcopy(new_results)
        else:
            for all_or_primer in results_totals:
                for row in results_totals[all_or_primer]:
                    for col, value in new_results[all_or_primer][row].items():
                        results_totals[all_or_primer][row][col] += value
        logging.info(
            f"Finished processing run {run_name} {completed+1}/{len(manifest_data)}"
        )

    logging.info("=" * 60)
    logging.info("Finished processing all runs")
    per_run_results_tsv = os.path.join(outdir, "results_per_run.tsv")
    logging.info(f"Writing per run results TSV {per_run_results_tsv}")
    with open(per_run_results_tsv, "w") as f:
        for row in per_sample_tsv_lines:
            print(*row, sep="\t", file=f)
    results_tsv = os.path.join(outdir, "results.tsv")
    logging.info(f"Writing summary results TSV {results_tsv}")
    msa.write_stats_summary_tsv(results_totals, results_tsv)
    results_json = os.path.join(outdir, "results.json")
    logging.info(f"Writing summary results JSON {results_json}")
    with open(results_json, "w") as f:
        json.dump(per_run_results, f, indent=2)
    logging.info("Finished")

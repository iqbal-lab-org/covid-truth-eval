import csv
import json
import logging
import os

from cte import one_run_evaluator


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
            if d["truth_vcf"] not in manifest_data:
                manifest_data[d["truth_vcf"]] = {}
            vcf_d = manifest_data[d["truth_vcf"]]
            assert d["name"] not in vcf_d
            vcf_d[d["name"]] = {
                "fasta_to_eval": d["eval_fasta"],
                "primers": d["primers"],
            }

    return manifest_data


def all_results_dict_to_tsv(all_results, tsv_out):
    snp_indels = ["snp", "indel"]
    tp_etc = ["TP", "FP", "FN"]
    cols = [f"{x}_{y}" for x in snp_indels for y in tp_etc]
    cols = cols + [f"primer_{x}" for x in cols]

    with open(tsv_out, "w") as f:
        print("Name", "Truth_vcf", *cols, sep="\t", file=f)
        for truth_vcf in sorted(all_results):
            for name, d in sorted(all_results[truth_vcf].items()):
                print(
                    name,
                    truth_vcf,
                    *[d["snp"][x] for x in tp_etc],
                    *[d["indel"][x] for x in tp_etc],
                    *[d["primer_snp"][x] for x in tp_etc],
                    *[d["primer_indel"][x] for x in tp_etc],
                    sep="\t",
                    file=f,
                )


def evaluate_runs(manifest_tsv, outdir, ref_fasta):
    os.mkdir(outdir)
    processing_dir = os.path.join(outdir, "Processing")
    os.mkdir(processing_dir)
    manifest_data = load_manifest_tsv(manifest_tsv)
    truth_vcf_to_dir = {}
    all_results = {}

    for truth_vcf in manifest_data:
        logging.info("=" * 60)
        logging.info(f"Start processing truth VCF file {truth_vcf}")

        truth_root_dir = os.path.join(processing_dir, str(len(truth_vcf_to_dir)))
        truth_vcf_to_dir[truth_vcf] = truth_root_dir
        os.mkdir(truth_root_dir)
        truth_files_dir = os.path.join(truth_root_dir, "Truth_files")
        evaluator = one_run_evaluator.OneSampleEvaluator(
            ref_fasta, truth_vcf, truth_files_dir
        )
        all_results[truth_vcf] = {}
        logging.info(f"Made truth files. Processing each run for VCF file {truth_vcf}")

        for run_name in manifest_data[truth_vcf]:
            logging.info("_" * 40)
            logging.info(f"Start processing {run_name} with truth VCF {truth_vcf}")
            run_outdir = os.path.join(truth_root_dir, run_name)
            fasta_to_eval = manifest_data[truth_vcf][run_name]["fasta_to_eval"]
            primers = manifest_data[truth_vcf][run_name]["primers"]
            all_results[truth_vcf][run_name] = evaluator.evaluate_one_fasta(
                run_outdir,
                fasta_to_eval,
                primers,
            )
            logging.info(f"Finished processing {run_name} with truth VCF {truth_vcf}")
        logging.info(f"Finished processing runs for truth VCF file {truth_vcf}")

    results_json = os.path.join(outdir, "results.json")
    logging.info(f"Writing JSON file of results {results_json}")
    with open(results_json, "w") as f:
        json.dump(all_results, f, indent=2)

    results_tsv = os.path.join(outdir, "results.tsv")
    logging.info(f"Writing TSV file of results {results_tsv}")
    all_results_dict_to_tsv(all_results, results_tsv)

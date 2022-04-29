#!/usr/bin/env python3

import argparse
import logging
import sys
import cte


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="cte", usage="cte <command> <options>", description="cte: test-covid-eval"
    )

    parser.add_argument("--version", action="version", version=cte.__version__)

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------- general options common to all tasks ------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )
    common_parser.add_argument(
        "--ref_fasta",
        help="Reference FASTA file. Default is to use covid MN908947.3 [%(default)s]",
        default=cte.built_in_data.COVID_REF,
        metavar="FILENAME",
    )
    common_parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory if it already exists",
    )
    common_parser.add_argument(
        "--outdir",
        required=True,
        help="REQUIRED. Output directory (will be created, or see --force)",
        metavar="FILENAME",
    )

    # ------------------------ eval_one_run -----------------------------------
    subparser_eval_one_run = subparsers.add_parser(
        "eval_one_run",
        parents=[common_parser],
        help="Evaluate one run",
        usage="cte eval_one_run [options] <truth files options> --outdir out --fasta_to_eval to_eval.fa --primers name",
        description="Evaluate one consensus sequence",
    )
    subparser_eval_one_run.add_argument(
        "--truth_vcf",
        help="Truth VCF file (with respect to the reference genome). If not provided, must provide --truth_fasta. If this option and --truth_fasta is used, then only dropped amplicon entries are used from the truth_vcf file",
        metavar="FILENAME",
    )
    subparser_eval_one_run.add_argument(
        "--truth_fasta",
        help="Truth FASTA file. If not provided, must provide --truth_vcf",
        metavar="FILENAME",
    )
    subparser_eval_one_run.add_argument(
        "--fasta_to_eval",
        required=True,
        help="REQUIRED. FASTA file of consensus sequence to be evaluated",
        metavar="FILENAME",
    )
    scheme_names = ",".join(cte.built_in_data.COVID_SCHEME_NAMES)
    subparser_eval_one_run.add_argument(
        "--primers",
        required=True,
        help=f"REQUIRED. TSV file of primers (in 'viridian_workflow' format). Or use a built-in scheme by providing one of these names: {scheme_names}",
        metavar="SCHEME_NAME/FILENAME",
    )
    subparser_eval_one_run.set_defaults(func=cte.tasks.eval_one_run.run)

    # ------------------------ eval_runs --------------------------------------
    subparser_eval_runs = subparsers.add_parser(
        "eval_runs",
        parents=[common_parser],
        help="Evaluate multiple consensus sequences",
        usage="cte eval_runs [options] --outdir out manifest.tsv",
        description="Evaluate multiple consensus sequences",
    )
    subparser_eval_runs.add_argument(
        "manifest_tsv",
        help="TSV file containing files to be evaluated",
    )
    subparser_eval_runs.set_defaults(func=cte.tasks.eval_runs.run)

    args = parser.parse_args()

    logging.basicConfig(
        format="[%(asctime)s cte %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if hasattr(args, "debug") and args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        if (
            args.func == cte.tasks.eval_one_run.run
            and args.truth_fasta is None
            and args.truth_vcf is None
        ):
            print(
                "Must use --truth_fasta or --truth_vcf. Cannot continue",
                file=sys.stderr,
            )
            sys.exit(1)
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

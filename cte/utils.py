import logging
import operator
import subprocess
import sys

import cluster_vcf_records
import pyfastaq


def syscall(command):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code: {completed_process.returncode}")
    if completed_process.returncode != 0:
        print("Error running this command:", command, file=sys.stderr)
        print("Return code:", completed_process.returncode, file=sys.stderr)
        print(
            "Output from stdout:", completed_process.stdout, sep="\n", file=sys.stderr
        )
        print(
            "Output from stderr:", completed_process.stderr, sep="\n", file=sys.stderr
        )
        raise Exception("Error in system call. Cannot continue")

    logging.debug(f"stdout:\n{completed_process.stdout.rstrip()}")
    logging.debug(f"stderr:\n{completed_process.stderr.rstrip()}")
    return completed_process


def load_single_seq_fasta(infile):
    d = {}
    pyfastaq.tasks.file_to_dict(infile, d)
    if len(d) != 1:
        raise Exception(
            f"Expected exatcly 1 sequence in {infile} but got {len(d)} sequences"
        )
    ref = list(d.values())[0]
    ref.id = ref.id.split()[0]
    return ref


def vcf_file_to_sorted_list(vcf_file):
    """Loads VCF file. Assumption is only one reference sequence.
    Returns a list of variants sorted by POS"""
    records = []
    ref_name = None

    header_lines, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(
        vcf_file
    )
    for record in vcf_records:
        if ref_name is None:
            ref_name = record.CHROM
        elif ref_name != record.CHROM:
            raise Exception(
                f"Only one CHROM allowed in VCF file, found more than one: {ref_name}, {record.CHROM}"
            )
        records.append(record)

    records.sort(key=operator.attrgetter("POS"))
    return records


def apply_variants_to_genome(vcf_file, out_fasta, ref_seq=None, ref_fasta=None):
    """Takes the variants in vcf_file, and applies them to the associated
    reference genome in ref_fasta. Writes a new file out_fasta that has those
    variants applied"""
    if ref_seq is None:
        assert ref_fasta is not None
        ref_seq = load_single_seq_fasta(ref_fasta)
    else:
        assert ref_fasta is None
    vcf_records = vcf_file_to_sorted_list(vcf_file)
    new_seq = list(ref_seq.seq)
    previous_ref_start = None

    # Applying indels messes up the coords of any subsequent variant,
    # so start at the end and work backwards
    for vcf_record in reversed(vcf_records):
        genotype = set(vcf_record.FORMAT["GT"].split("/"))
        if len(genotype) != 1:
            # FIXME
            raise NotImplementedError()
        try:
            allele_index = int(genotype.pop())
        except:
            raise Exception(
                f"Genotype must be int(s). Cannot use this line of VCF file: {vcf_record}"
            )
        if allele_index == 0:
            logging.warn(
                f"Genotype of zero in truth VCF. Is that deliberate? Ignoring the line: {vcf_record}"
            )
            continue

        # If the current record overlaps the previous one, stop.
        # Should not happen in a nice truth VCF.
        if (
            previous_ref_start is not None
            and vcf_record.ref_end_pos() >= previous_ref_start
        ):
            raise Exception(
                f"VCF record overlaps another record, cannot continue: {vcf_record}"
            )

        previous_ref_start = vcf_record.POS
        allele = vcf_record.ALT[allele_index - 1]
        start, end = vcf_record.POS, vcf_record.ref_end_pos() + 1
        assert ref_seq[start:end] == "".join(new_seq[start:end])
        new_seq[start:end] = [allele]

    with open(out_fasta, "w") as f:
        new_seq = pyfastaq.sequences.Fasta(f"{ref_seq.id}.mutated", "".join(new_seq))
        print(new_seq, file=f)

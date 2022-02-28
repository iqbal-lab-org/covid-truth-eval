import logging
import operator
import subprocess
import sys

import cluster_vcf_records
import pyfastaq

from cte import iupac


def syscall(command):
    logging.info(f"Run command: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    logging.info(f"Return code = {completed_process.returncode} from: {command}")
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
        if "DROPPED_AMP" in vcf_record.FILTER:
            continue

        try:
            genotypes = {int(x) for x in vcf_record.FORMAT["GT"].split("/")}
        except:
            raise Exception(
                f"Genotype(s) must be int(s). Cannot use this line of VCF file: {vcf_record}"
            )

        if "UNSURE" in vcf_record.FILTER:
            allele = "N" * len(vcf_record.REF)
        elif len(genotypes) == 1:
            allele_index = genotypes.pop()
            if allele_index == 0:
                logging.warning(
                    f"Genotype of zero in truth VCF. Is that deliberate? Ignoring the line: {vcf_record}"
                )
                continue
            allele = vcf_record.ALT[allele_index - 1]
        else:
            # We only allow HET SNPs, not indels.
            if len(vcf_record.REF) != 1 or {len(x) for x in vcf_record.ALT} != {1}:
                raise Exception(
                    f"Heterozygous SNPs are not allowed. Cannot continue. VCF line:\n{vcf_record}"
                )

            alleles = []
            if 0 in genotypes:
                alleles.append(vcf_record.REF)
            for i in genotypes:
                if i > 0:
                    alleles.append(vcf_record.ALT[i - 1])
            allele = iupac.rev_ambiguous_codes["".join(sorted(alleles))]

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
        start, end = vcf_record.POS, vcf_record.ref_end_pos() + 1
        assert ref_seq[start:end] == "".join(new_seq[start:end])
        new_seq[start:end] = [allele]

    with open(out_fasta, "w") as f:
        new_seq = pyfastaq.sequences.Fasta(f"{ref_seq.id}.mutated", "".join(new_seq))
        print(new_seq, file=f)


def string_to_gaps(string, min_gap_length):
    gaps = []
    gap_start = None
    for i, base in enumerate(string):
        if base == "N":
            if gap_start is None:
                gap_start = i
        elif gap_start is not None:
            if i - gap_start >= min_gap_length:
                gaps.append((gap_start, i - 1))
            gap_start = None

    if gap_start is not None and i - gap_start + 1 >= min_gap_length:
        gaps.append((gap_start, i))

    return gaps

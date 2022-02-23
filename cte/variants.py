import logging

import cluster_vcf_records

from cte import iupac


def vcf_record_overlaps_primer(record, amplicon_scheme):
    return amplicon_scheme.overlaps_primer(record.POS, record.ref_end_pos())


def vcf_record_is_snp(record):
    if len(record.REF) != 1:
        return False
    return {len(x) for x in record.ALT} == {1}


def vcf_record_is_truth_dropped_amplicon(record, from_truth):
    if not (from_truth and "DROPPED_AMP" in record.FILTER):
        return False

    assert "AMP_START" in record.INFO and "AMP_END" in record.INFO
    assert int(record.INFO["AMP_START"]) == record.POS
    return True


def vcf_record_is_het(record, from_truth):
    if from_truth:
        return "HET" in record.FILTER or "UNSURE" in record.FILTER
    else:
        # if it's not a truth VCF record, then it was made from MSA.
        # It can be het if alt allele has an ambiguous IUPAC code
        return any([x in iupac.ambiguous_codes for x in record.ALT[0]])


def vcf_record_alt_is_all_N(record):
    return set(record.ALT[0]) == {"N"}


def vcf_record_is_unsure(record, from_truth):
    return from_truth and "UNSURE" in record.FILTER


def vcf_record_amplicon_dropout_matches(record, amplicon_scheme, min_percent_overlap=90.0):
    """Returns list of names of matching amplicons"""
    if len(record.ALT) != 1 or not vcf_record_alt_is_all_N(record):
        return None

    if 100 * (abs(len(record.REF) - len(record.ALT[0])) / len(record.REF)) <= min_percent_overlap:
        return None

    matching_amplicons = amplicon_scheme.overlapping_amplicons(record.POS, record.ref_end_pos())
    if len(matching_amplicons) == 0:
        return None

    good_matches = []
    for d in matching_amplicons:
        overlap_start = max(d["start"], record.POS)
        overlap_end = min(d["end"], record.ref_end_pos())
        percent_overlap = 100 * (overlap_end - overlap_start + 1) / (d["end"] - d["start"] + 1)
        if percent_overlap >= min_percent_overlap:
            good_matches.append(d)
    if len(good_matches) > 0:
        return good_matches
    else:
        return None



#def qry_to_ref_coords_from_msa(ref_aln, qry_aln):
#    """Returns list of 0-based query position -> ref position"""
#    coords = []
#    assert len(ref_aln) == len(qry_aln)
#    ref_coord = 0
#    for i in range(len(ref_aln)):
#        if qry_aln[i] != "-":
#            coords.append(ref_coord)
#        if ref_aln[i] != "-":
#            ref_coord += 1
#    return coords



class Variant:
    def __init__(self, vcf_record, amplicon_scheme, from_truth):
        self.vcf_record = vcf_record
        self.amplicon_scheme = amplicon_scheme
        self.is_truth_dropped_amplicon = vcf_record_is_truth_dropped_amplicon(self.vcf_record, from_truth)
        self.is_snp = vcf_record_is_snp(self.vcf_record)
        self.in_primer = vcf_record_overlaps_primer(self.vcf_record, self.amplicon_scheme)
        self.is_het = vcf_record_is_het(self.vcf_record, from_truth)
        self.is_unsure = vcf_record_is_unsure(self.vcf_record, from_truth)
        self.alt_is_all_N = vcf_record_alt_is_all_N(self.vcf_record)
        self.amplicons_dropped = vcf_record_amplicon_dropout_matches(self.vcf_record, self.amplicon_scheme, min_percent_overlap=90.0)

    def __hash__(self):
        return hash(str(self.vcf_record))

    def __eq__(self, other):
        return type(other) is type(self) and self.vcf_record == other.vcf_record

    def ref_start(self):
        return self.vcf_record.POS

    def ref_end(self):
        if self.is_truth_dropped_amplicon:
            return int(self.vcf_record.INFO["AMP_END"])
        else:
            return self.vcf_record.ref_end_pos()

    def ref_length(self):
        return self.ref_end() - self.ref_start() + 1

    def truth_type_for_reporting(self):
        if self.is_truth_dropped_amplicon:
            return "Dropped_amplicon"
        elif self.is_unsure:
            return "unknown_truth"
        elif self.is_snp:
            if self.is_het:
                return "SNP_true_mixed"
            else:
                return "SNP_true_alt"
        else:
            return "True_indel"

    def expanded_alleles(self):
        return iupac.expand_nucleotide_list(self.vcf_record.ALT + [self.vcf_record.REF])

    def is_same_variant(self, other, ref_seq):
        if self.is_snp and other.is_snp:
            is_same = self.vcf_record.POS == other.vcf_record.POS and self.vcf_record.REF == other.vcf_record.REF
            if is_same:
                if self.is_het or other.is_het:
                    is_same = self.expanded_alleles() == other.expanded_alleles()
                else:
                    is_same = self.vcf_record.ALT == other.vcf_record.ALT
        else:
            is_same = self.vcf_record.is_the_same_indel(other.vcf_record, ref_seq)
        logging.debug(f"is same variant: {is_same}. Variants compared:")
        logging.debug(f"  {self.vcf_record}")
        logging.debug(f"  {other.vcf_record}")
        return is_same

    def amplicon_dropout_matches(self):
        return vcf_record_amplicon_dropout_matches(self.vcf_record, self.amplicon_scheme)


class VariantSet:
    def __init__(self, vcf_file, amplicon_scheme, from_truth):
        self.vars = {} # 0-based ref POS -> Variant
        header, vcf_records = cluster_vcf_records.vcf_file_read.vcf_file_to_list(vcf_file)
        self.ref_pos_to_var_pos = {}
        for vcf_record in vcf_records:
            assert vcf_record.POS not in self.vars
            self.vars[vcf_record.POS] = Variant(vcf_record, amplicon_scheme, from_truth)
            for i in range(vcf_record.POS, vcf_record.ref_end_pos() + 1):
                assert i not in self.ref_pos_to_var_pos
                self.ref_pos_to_var_pos[i] = vcf_record.POS


    def get_overlapping_variants(self, qry_var):
        """Returns a list of all variants in this variant set whose REF
        allele intersects the REF allele of qry_var"""
        match_keys = {self.ref_pos_to_var_pos[i] for i in range(qry_var.ref_start(), qry_var.ref_end()+1) if i in self.ref_pos_to_var_pos}
        return [self.vars[x] for x in match_keys]


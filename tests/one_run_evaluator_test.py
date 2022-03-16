import filecmp
import os
import pytest

from cluster_vcf_records import vcf_record

from cte import built_in_data, msa, one_run_evaluator, utils


def different_nucleotide(nuc_in):
    return {"A": "C", "C": "G", "G": "T", "T": "C"}[nuc_in]


def make_test_data(outdir, test_type):
    os.mkdir(outdir)
    truth_vcf_records = []
    files = {
        "ref_fasta": built_in_data.COVID_REF,
        "truth_vcf": os.path.join(outdir, "truth.vcf"),
        "fasta_to_eval": os.path.join(outdir, "to_eval.fa"),
        "expect_stats_tsv": os.path.join(outdir, "expect_stats.tsv"),
        "expect_stats_json": os.path.join(outdir, "expect_stats.json"),
        "expect_per_record_tsv": os.path.join(outdir, "expect_record_breakdown.tsv"),
        "primers_tsv": built_in_data.COVID_PRIMER_TSVS["COVID-ARTIC-V3"],
    }
    ref_seq = utils.load_single_seq_fasta(files["ref_fasta"])
    eval_seq = list(ref_seq)
    stats = msa.make_empty_stats_dict()
    total_ref = 29836
    total_in_primer = 4908

    # always add a SNP before the first amplicon and after the last one -
    # should get ignored
    eval_seq[19] = different_nucleotide(eval_seq[19])
    eval_seq[29899] = different_nucleotide(eval_seq[29899])

    if test_type == "No_vars":
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_in_primer
    elif test_type == "true_ref":
        # Add a false positive SNP non-het call.
        assert ref_seq[199] == "T"
        eval_seq[199] = "A"
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_alt] += 1

        # Add a false positive SNP non-N het call. Is in a primer
        eval_seq[399] = "Y"
        assert ref_seq[399] == "T"
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_IUPAC] += 1
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_IUPAC] += 1

        # Add a false positive N call.
        eval_seq[599] = "N"
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_N] += 1

        # Add a false positive two consecutive Ns
        eval_seq[699] = "N"
        eval_seq[700] = "N"
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_N] += 2

        # Add a false positve insertion call.
        eval_seq[799] = eval_seq[799] + "A"
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] += 1

        # Add a false positve deletion call.
        eval_seq[999] = ""
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] += 1

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref - 6
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_in_primer - 1
        )
    elif test_type == "SNP_true_alt":
        # A true SNP (not het) called as ref
        assert ref_seq[199] == "T"
        record = vcf_record.VcfRecord("MN908947.3\t200\t.\tT\tA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_ref] += 1
        truth_vcf_records.append(record)

        # A true SNP (not het) called correct alt. Also is in a primer
        assert ref_seq[399] == "T"
        eval_seq[399] = "A"
        record = vcf_record.VcfRecord("MN908947.3\t400\t.\tT\tA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_correct_alt] += 1
        stats["Primer"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_correct_alt] += 1
        truth_vcf_records.append(record)

        # A true SNP (not het) called wrong alt
        assert ref_seq[599] == "G"
        eval_seq[599] = "T"
        record = vcf_record.VcfRecord("MN908947.3\t600\t.\tG\tA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_wrong_alt] += 1
        truth_vcf_records.append(record)

        # A true SNP (not het) called as N
        assert ref_seq[799] == "G"
        eval_seq[799] = "N"
        record = vcf_record.VcfRecord("MN908947.3\t800\t.\tG\tA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_N] += 1
        truth_vcf_records.append(record)

        # A true SNP called as ambiguous code
        assert ref_seq[999] == "T"
        eval_seq[999] = "Y"
        record = vcf_record.VcfRecord("MN908947.3\t1000\t.\tT\tA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_wrong_IUPAC] += 1
        truth_vcf_records.append(record)

        # A true SNP called as indel
        assert ref_seq[1199] == "G"
        eval_seq[1199] = "GG"
        record = vcf_record.VcfRecord("MN908947.3\t1200\t.\tG\tA\t42.42\t.\t.\tGT\t1/1")
        # The alignment results in a bp indel, which is incorrect. And at the
        # SNP it's the same as the ref.
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] += 1
        stats["All"][msa.StatRow.SNP_true_alt][msa.StatCol.Called_ref] += 1
        truth_vcf_records.append(record)

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref - 6
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_in_primer - 1
        )
    elif test_type == "SNP_true_mixed":
        # True positive HET where ref is called
        assert ref_seq[199] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t200\t.\tT\tA\t42.42\tHET\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_ref] += 1
        truth_vcf_records.append(record)

        # A true het call called as an N. Is in primer
        assert ref_seq[399] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t400\t.\tT\tG\t42.42\tHET\t.\tGT\t0/1"
        )
        eval_seq[399] = "N"
        stats["All"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_N] += 1
        stats["Primer"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_N] += 1
        truth_vcf_records.append(record)

        # True positive HET where neither allele is the ref
        # And this is in a primer, so also adds to P_SNP_true_mixed
        assert ref_seq[599] == "G"
        eval_seq[599] = "Y"
        record = vcf_record.VcfRecord(
            "MN908947.3\t600\t.\tG\tC,T\t42.42\tHET\t.\tGT\t1/2"
        )
        stats["All"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_correct_IUPAC] += 1
        truth_vcf_records.append(record)

        # True positive HET called with wrong amibiguous code
        assert ref_seq[799] == "G"
        eval_seq[799] = "R"  # = A or G
        record = vcf_record.VcfRecord(
            "MN908947.3\t800\t.\tG\tC\t42.42\tHET\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_wrong_IUPAC] += 1
        truth_vcf_records.append(record)

        ## A true HET called as indel
        assert ref_seq[999] == "T"
        eval_seq[999] = "TA"
        record = vcf_record.VcfRecord(
            "MN908947.3\t1000\t.\tT\tA\t42.42\tHET\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.SNP_true_mixed][msa.StatCol.Called_ref] += 1
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] += 1
        truth_vcf_records.append(record)

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref - 5
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_in_primer - 1
        )
    elif test_type == "Unknown_truth":
        # Unknown position called as ref
        assert ref_seq[199] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t199\t.\tT\tA\t42.42\tUNSURE\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.Called_ref] += 1
        truth_vcf_records.append(record)

        # Unknown position called as snp
        assert ref_seq[399] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t400\t.\tT\tG\t42.42\tUNSURE\t.\tGT\t0/1"
        )
        eval_seq[399] = "G"
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.Called_other] += 1
        stats["Primer"][msa.StatRow.Unknown_truth][msa.StatCol.Called_other] += 1
        truth_vcf_records.append(record)

        # Unknown position called as N
        assert ref_seq[599] == "G"
        eval_seq[599] = "N"
        record = vcf_record.VcfRecord(
            "MN908947.3\t600\t.\tG\tC,T\t42.42\tUNSURE\t.\tGT\t1/2"
        )
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.Called_N] += 1
        truth_vcf_records.append(record)

        # Unknown position called as amibiguous code
        assert ref_seq[799] == "G"
        eval_seq[799] = "R"  # = A or G
        record = vcf_record.VcfRecord(
            "MN908947.3\t800\t.\tG\tC\t42.42\tUNSURE\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.Called_other] += 1
        truth_vcf_records.append(record)

        # Unknown position called as indel
        assert ref_seq[999] == "T"
        eval_seq[999] = "TA"
        record = vcf_record.VcfRecord(
            "MN908947.3\t1000\t.\tT\tA\t42.42\tUNSURE\t.\tGT\t0/1"
        )
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_indel] += 1
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.Called_ref] += 1
        truth_vcf_records.append(record)

        # Unknown position near then end, which doesn't get assembled
        assert ref_seq[29833] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t29834\t.\tT\tN\t42.42\tUNSURE\t.\tGT\t1/1"
        )
        truth_vcf_records.append(record)
        eval_seq[29830:] = "N" * len(eval_seq[29830:])
        stats["All"][msa.StatRow.Unknown_truth][msa.StatCol.No_call_genome_ends] += 1
        stats["All"][msa.StatRow.True_ref][msa.StatCol.No_call_genome_ends] += 35
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.No_call_genome_ends] += 30

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref - 41
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_in_primer - 31
        )
    elif test_type == "True_indel":
        # A true indel called as ref
        assert ref_seq[199] == "T"
        record = vcf_record.VcfRecord("MN908947.3\t200\t.\tT\tTA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_wrong_indel] += 1
        truth_vcf_records.append(record)

        # A true indel called correct alt. Also is in a primer
        assert ref_seq[399] == "T"
        eval_seq[399] = "TA"
        record = vcf_record.VcfRecord("MN908947.3\t400\t.\tT\tTA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_correct_alt] += 1
        stats["Primer"][msa.StatRow.True_indel][msa.StatCol.Called_correct_alt] += 1
        truth_vcf_records.append(record)

        # A true indel called as an incorrect SNP
        assert ref_seq[599] == "G"
        eval_seq[599] = "T"
        record = vcf_record.VcfRecord("MN908947.3\t600\t.\tG\tGA\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_alt] += 1
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_wrong_indel] += 1
        truth_vcf_records.append(record)

        # A true indel called as N
        assert ref_seq[799] == "G"
        assert ref_seq[800] == "G"
        eval_seq[799] = "N"
        eval_seq[800] = "N"
        record = vcf_record.VcfRecord("MN908947.3\t800\t.\tGG\tG\t42.42\t.\t.\tGT\t1/1")
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_N] += 1
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_N] += 1
        truth_vcf_records.append(record)

        # A true indel called as ambiguous code
        assert ref_seq[999] == "T"
        eval_seq[999] = "Y"
        record = vcf_record.VcfRecord(
            "MN908947.3\t1000\t.\tT\tTA\t42.42\t.\t.\tGT\t1/1"
        )
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_wrong_indel] += 1
        truth_vcf_records.append(record)

        # True indel called with wrong indel
        assert ref_seq[1199] == "G"
        eval_seq[1199] = "GA"
        record = vcf_record.VcfRecord(
            "MN908947.3\t1200\t.\tG\tGC\t42.42\t.\t.\tGT\t1/1"
        )
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_wrong_IUPAC] += 1
        stats["All"][msa.StatRow.True_indel][msa.StatCol.Called_wrong_indel] += 1
        truth_vcf_records.append(record)

        # Indel near then end, which doesn't get assembled
        assert ref_seq[29833] == "T"
        record = vcf_record.VcfRecord(
            "MN908947.3\t29834\t.\tT\tTA\t42.42\tPASS\t.\tGT\t1/1"
        )
        truth_vcf_records.append(record)
        eval_seq = eval_seq[:29830]
        stats["All"][msa.StatRow.True_indel][msa.StatCol.No_call_genome_ends] += 1
        stats["All"][msa.StatRow.True_ref][msa.StatCol.No_call_genome_ends] += 36
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.No_call_genome_ends] += 30

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = total_ref - 40
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_in_primer - 30
        )
    elif test_type == "Dropped_amplicon":
        # Drop amplicon 2. Is at 320-725 (0-based)
        record = vcf_record.VcfRecord(
            "MN908947.3\t321\t.\tG\tN\t42.42\tDROPPED_AMP\tAMP_START=320;AMP_END=725\tGT\t1/1"
        )
        truth_vcf_records.append(record)
        eval_seq[320:726] = "N" * (726 - 320)
        stats["All"][msa.StatRow.Dropped_amplicon][msa.StatCol.Dropped_amplicon] += (
            726 - 320
        )
        stats["Primer"][msa.StatRow.Dropped_amplicon][
            msa.StatCol.Dropped_amplicon
        ] += 91

        # Have consensus incorrectly drop most of amplicon 4. Is at 943-1336
        eval_seq[950:1300] = "N" * 350
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Dropped_amplicon] += 350
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Dropped_amplicon] += 61

        # Update the all ref counts
        stats["All"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = (
            total_ref - (726 - 320) - 350
        )
        stats["Primer"][msa.StatRow.True_ref][msa.StatCol.Called_ref] = 4756
    else:
        raise NotImplementedError(f"Not implemented: {test_type}")

    with open(files["truth_vcf"], "w") as f:
        print(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1",
            *truth_vcf_records,
            sep="\n",
            file=f,
        )
    with open(files["fasta_to_eval"], "w") as f:
        print(">to_eval", "".join(eval_seq), sep="\n", file=f)
    msa.write_stats_summary_tsv(stats, files["expect_stats_tsv"])
    return files


def eval_one_fasta_test(outdir, test_type):
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)
    data_dir = os.path.join(outdir, "data")
    results_dir = os.path.join(outdir, "results")
    files = make_test_data(
        data_dir,
        test_type,
    )
    one_run_evaluator.eval_one_fasta(
        results_dir,
        files["fasta_to_eval"],
        files["ref_fasta"],
        files["truth_vcf"],
        files["primers_tsv"],
        debug=True,
        force=False,
    )
    got_stats_tsv = os.path.join(results_dir, "results.tsv")
    assert filecmp.cmp(got_stats_tsv, files["expect_stats_tsv"], shallow=False)
    assert os.path.exists(os.path.join(results_dir, "per_position.tsv"))
    assert os.path.exists(os.path.join(results_dir, "results.json"))
    utils.syscall(f"rm -r {outdir}")


def test_eval_one_fasta_no_vars():
    eval_one_fasta_test("tmp.eval_one_fasta_no_vars", test_type="No_vars")


def test_eval_one_fasta_true_ref():
    eval_one_fasta_test("tmp.eval_one_fasta_true_ref", test_type="true_ref")


def test_eval_one_fasta_SNP_true_alt():
    eval_one_fasta_test("tmp.eval_one_fasta_SNP_true_alt", test_type="SNP_true_alt")


def test_eval_one_fasta_SNP_true_mixed():
    eval_one_fasta_test("tmp.eval_one_fasta_SNP_true_mixed", test_type="SNP_true_mixed")


def test_eval_one_fasta_unknown_truth():
    eval_one_fasta_test("tmp.eval_one_fasta_unknown_truth", test_type="Unknown_truth")


def test_eval_one_fasta_True_indel():
    eval_one_fasta_test("tmp.eval_one_fasta_True_indel", test_type="True_indel")


def test_eval_one_fasta_Dropped_amplicon():
    eval_one_fasta_test(
        "tmp.eval_one_fasta_Dropped_amplicon", test_type="Dropped_amplicon"
    )

# covid-truth-eval
Evaluate accuracy of covid assemblies where truth is available

## Installation

Minimal instructions are below. Please see the [wiki page](https://github.com/iqbal-lab-org/covid-truth-eval/wiki)
for more details.


### Docker

Get a Docker image of the latest release:

```
docker pull ghcr.io/iqbal-lab-org/cte:latest
```

All Docker images are listed in the
[packages page](https://github.com/iqbal-lab-org/covid-truth-eval/pkgs/container/covid-truth-eval).

Alternatively, build your own Docker image:

```
sudo docker build --network=host .
```

### Singularity

[Releases](https://github.com/iqbal-lab-org/covid-truth-eval/releases)
include a Singularity image to download.
Each release has a file called `cte_vX.Y.Z.img`, where `X.Y.Z` is the release version.

Alternatively, build your own Singularity image:

```
singularity build cte.simg Singularity.def
```


### From source

Dependencies:

* Python 3
* [mafft](https://mafft.cbrc.jp/alignment/software/) installed and in your `$PATH`.

Install by cloning this repository (or downloading the latest release), and
running:

```
python3 -m pip install .
```


## Quick start

To evaluate one SARS-CoV-2 consensus sequence, you will need:
1. A VCF file of the "truth" calls `truth.vcf`
2. The consensus sequence to evalaute in a FASTA file `cons.fa`
3. The primer scheme. Currently supported: COVID-ARTIC-V3, COVID-ARTIC-V4, COVID-MIDNIGHT-1200.
   Or use your own TSV file of primers in [Viridian Workflow format](https://github.com/iqbal-lab-org/viridian_workflow/wiki/Amplicon-schemes).
   
The program is called `cte` and is installed in the Docker and Singularity containers, and gets installed by `pip`.

Example, assuming primer scheme COVID-ARTIC-V4:
```
cte eval_one_run \
  --outdir OUT \
  --truth_vcf truth.vcf \
  --fasta_to_eval cons.fa \
  --primers COVID-ARTIC-V4
```

The output files are:
1. `results.tsv` - a TSV file of counts of the truth bases vs what was called in the consensus. The same information is also put in a JSON file `results.json`.
2. `per_position.tsv` - a TSV file, one line per reference position. It shows the multiple alignment of the reference, truth (inferred from the truth VCF file), and the sequence being evaluated. At each position the assigned category of truth and called bases is shown, where the categories are the same as those used in `results.tsv`.

The files are described in detail in the [output files](https://github.com/iqbal-lab-org/covid-truth-eval/wiki/Output-files) documentation.

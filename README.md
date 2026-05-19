# Plasmid Typing Pipeline

Pipeline to type the genome of plasmid fastas files using several tools and convert it into a single report file. This report can be put into the shiny_upload_template.csv to have consistent outputs with the [PlasmidNL Shiny application](https://apps.rivm.nl/bsr-ids-ienv/plasmidnl/), also available as [standalone R/Shiny](https://gitlab.com/gidjes/plasmidnl).

# Installation

Unfortunately, several tools used in this pipeline are incompatible in the same environment. So several seperate conda environments have to be created in order to get all parts to work. 

```
git clone https://gitlab.com/mmb-umcu/plasmidnl_typing.git
cd plasmidnl_typing

chmod +x scripts/init.sh
./scripts/init.sh

conda activate plasmidnl_typing
```

# Preperation

All fasta files should be placed in the input directory.

## Command-Line Arguments

The pipeline accepts the following command-line options.

### `--input` (required)

Path to the input directory containing plasmid FASTA files.

Example:

    --input fastas/plasmids/

------------------------------------------------------------------------

### `--card`

Enable CARD (Comprehensive Antibiotic Resistance Database) resistance
typing.

If provided, CARD-based resistance gene detection will be run.\
If omitted, CARD typing is skipped.

Example:

    --card

------------------------------------------------------------------------

### `--skip_amrfinder`

Disable resistance gene detection using AMRFinderPlus.

By default, AMRFinderPlus results are included.\
Use this flag to ignore AMRFinderPlus-derived resistance genes.

Example:

    --skip_amrfinder

------------------------------------------------------------------------

### `--jobs`

Number of parallel worker processes to use.

-   Type: integer\
-   Default: `1`

Example:

    --jobs 8

------------------------------------------------------------------------

### `--rerun-failed`

Rerun failed sequence from previous run. Only works if the pipeline has
been run previously.

Example:

    --rerun-failed

------------------------------------------------------------------------

## Example Usage

Run the pipeline with CARD typing enabled, AMRFinderPlus disabled, and 8
parallel jobs:

```bash
    type_plasmids \
      --input data/plasmids/ \
      --card \
      --skip_amrfinder \
      --jobs 8
```

# Components

## PlasmidFinder [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4068535/)]
Runs [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) on the input sequences and runs a subsequent script to extract replicon data.

## MOB-typer [[2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/), [3](https://pubmed.ncbi.nlm.nih.gov/32969786/)]
MOB-typer from [MOB-suite](https://github.com/phac-nml/mob-suite) is run to determine plasmid mobility.

## mge-cluster [[4](https://academic.oup.com/nargab/article/5/3/lqad066/7222077)]
Plasmids are clustered using [mge-cluster](https://gitlab.com/sirarredondo/mge-cluster). In basic usage, it applies a mge-cluster scheme constructed from Dutch plasmids obtained through the surveillance for [carbapenemase-producing organisms](http://dx.doi.org/10.5281/ZENODO.18920264). However, a custom mge-cluster scheme can be placed in data/mge_scheme_custom and used through the `--custom_mge_scheme` flag

## ResFinder [[5](https://pmc.ncbi.nlm.nih.gov/articles/PMC8914360/)]
Resistance genes are annotated using [ResFinder](https://genepi.food.dtu.dk/resfinder). Two additional tools can be run for their resistance gene content (see below). However, for overlapping hits, ResFinder will be preferentially used.

## AMRFinderPlus [[6](https://www.nature.com/articles/s41598-021-91456-0)]
[AMRFinderPlus](https://github.com/ncbi/amr) is run for their extended databases:
- AMR
- Virulence
- Stress

The AMR genes found by AMRFinderPlus can be ignored in the final output by using the --skip_amrfinder flag

## CARD [[7](https://academic.oup.com/nar/article/51/D1/D690/6764414)]
The [Comprehensive Antibiotic Resistance Database (CARD)](https://github.com/arpcard) can be additionally run by using the --card flag, as the AMR databases are not fully redundant, for a comprehensive AMR overview. When overlapping hits are found, ResFinder hits are preferentially used. So, only additional genes detected are used in the final output.

# Citations
[[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4068535/)] PlasmidFinder and pMLST: in silico detection and typing of plasmids. Carattoli A, Zankari E, Garcia-Fernandez A, Volby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H. Antimicrob. Agents Chemother. 2014. April 28th. [Epub ahead of print]

[[2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6159552/)] Robertson, James, and John H E Nash. “MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies.” Microbial genomics vol. 4,8 (2018): e000206. doi:10.1099/mgen.0.000206

[[3](https://pubmed.ncbi.nlm.nih.gov/32969786/)] Robertson, James et al. “Universal whole-sequence-based plasmid typing and its utility to prediction of host range and epidemiological surveillance.” Microbial genomics vol. 6,10 (2020): mgen000435. doi:10.1099/mgen.0.000435

[[4](https://academic.oup.com/nargab/article/5/3/lqad066/7222077)] Arredondo-Alonso S, Gladstone RA, Pöntinen AK, Gama JA, Schürch AC, Lanza VF, et al. Mge-cluster: a reference-free approach for typing bacterial plasmids. NAR Genom Bioinform. 2023 Sep;5(3):lqad066.

[[5](https://pmc.ncbi.nlm.nih.gov/articles/PMC8914360/)] Bortolaia V, Kaas RF, Ruppe E, Roberts MC, Schwarz S, Cattoir V, et al. ResFinder 4.0 for predictions of phenotypes from genotypes. Journal of Antimicrobial Chemotherapy. 2020 Aug 11. PMID: 32780112 doi: 10.1093/jac/dkaa345 [Epub ahead of print]

[[6](https://www.nature.com/articles/s41598-021-91456-0)] Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984.

[[7](https://academic.oup.com/nar/article/51/D1/D690/6764414)] Alcock et al. 2023. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. Nucleic Acids Research, 51, D690-D699 [PMID 36263822]

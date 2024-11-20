# OpenPedCan-analysis
[![DOI](https://zenodo.org/badge/358689512.svg)](https://zenodo.org/badge/latestdoi/358689512)

The Open Pediatric Cancer (OpenPedCan) project is an open analysis effort that harmonizes pediatric cancer data from multiple sources, performs downstream cancer analyses on these data and provides them on PedcBioPortal and the [NCI's Molecular Targets Platform (MTP)](https://moleculartargets.ccdi.cancer.gov/.
For detailed methods, please see our [methods repository](https://github.com/rokitalab/OpenPedCan-manuscript).

To cite this work, please note the data release used in your work and cite the following:
1. OpenPBTA:
Shapiro, J.A., Gaonkar, K.S., Spielman, S.J., Savonen, C.L., Bethell, C.J., Jin, R., Rathi, K.S., Zhu, Y., Egolf, L.E., Farrow, B.K., et al. (2023). OpenPBTA: The Open Pediatric Brain Tumor Atlas. Cell Genom., 100340. [10.1016/j.xgen.2023.100340](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00115-5.pdf).
2. OpenPedCan:
Geng, Z. et al. The Open Pediatric Cancer Project. bioRxiv (2024) [10.1101/2024.07.09.599086](https://www.biorxiv.org/content/10.1101/2024.07.09.599086v1)
DOI for [all releases](https://zenodo.org/search?q=parent.id%3A6473912&f=allversions%3Atrue&l=list&p=1&s=10&sort=version): `10.5281/zenodo.6473912`. 

The OpenPedCan analyses currently include the following datasets, described more fully below:

- OpenPBTA
- TARGET
- Kids First Neuroblastoma
- Additional PBTA samples from CBTN, PNOC, Mioncoseq, Chordoma Foundation
- DGD (CHOP P30 Panel) - fusion and tumor/normal panels
- GTEx
- TCGA

Open Pediatric Brain Tumor Atlas (OpenPBTA)
In September of 2018, the [Children's Brain Tumor Network (CBTN)](https://cbtn.org/) released the [Pediatric Brain Tumor Atlas (PBTA)](https://cbtn.org/pediatric-brain-tumor-atlas/), a genomic dataset (whole genome sequencing, whole exome sequencing, RNA sequencing, proteomic, and clinical data) for nearly 1,000 tumors, available from the [Gabriella Miller Kids First Portal](https://kidsfirstdrc.org/).
In September of 2019, the Open Pediatric Brain Tumor Atlas (OpenPBTA) Project was launched.
OpenPBTA was a global open science initiative to comprehensively define the molecular landscape of tumors of 943 patients from the CBTN and the PNOC003 DIPG clinical trial from the [Pediatric Pacific Neuro-oncology Consortium](http://www.pnoc.us/) through real-time, collaborative analyses and [collaborative manuscript writing](https://github.com/AlexsLemonade/OpenPBTA-manuscript/) on GitHub, now published in [_Cell Genomics_](10.1016/j.xgen.2023.100340).
Additional PBTA data has been, and will be continually added to OpenPedCan.

Therapeutically Applicable Research to Generate Effective Treatments [(TARGET)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000218.v23.p8)
The Therapeutically Applicable Research to Generate Effective Treatments (TARGET) Initiative is an NCI-funded collection of disease-specific projects that seeks to identify the genomic changes of pediatric cancers. 'The overall goal is to collect genomic data to accelerate the development of more effective therapies.
OpenPedCan analyses include the seven diseases present in the TARGET dataset: Acute Lymphoblastic Leukemia (ALL), Acute Myeloid Leukemia (AML), Clear cell sarcoma of the kidney, Neuroblastoma, Osteosarcoma, Rhabdoid tumor, and Wilm’s Tumor.

Gabriella Miller Kids First Neuroblastoma [(Kids First)](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001436.v1.p1)
The Gabriella Miller Kids First Pediatric Research Program (Kids First) is a large-scale effort to accelerate research and gene discovery in pediatric cancers and structural birth defects.
The program includes whole genome sequencing (WGS) from patients with pediatric cancers and structural birth defects and their families.
OpenPedCan analyses include Neuroblastoma data from the Kids First project.

The Genotype-Tissue Expression [(GTEx)](https://gtexportal.org/home/)
GTEx project is an ongoing effort to build a comprehensive public data resource and tissue bank to study tissue-specific gene expression, regulation and their relationship with genetic variants.
Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq.
OpenPedCan project includes 17,382 GTEx RNA-Seq samples from GTEx v8 release, which span across 31 GTEx groups in the v11 release.

The Cancer Genome Atlas Program [(TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)
TCGA is a landmark cancer genomics program that molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types.
It is a joint effort between NCI and the National Human Genome Research Institute.
OpenPedCan project includes 10,414 TCGA RNA-Seq samples (716 normal and 9698 tumor) from [33 cancer types](https://github.com/d3b-center/OpenPedCan-analysis/blob/0a5c14705a385c99a6a16e34e932e94009b7a11c/analyses/molecular-subtyping-integrate/results/tcga_cancer_groups.tsv) in the release.

DGD [(CHOP P30 Panel)](https://www.chop.edu/cancer-panels)
CHOP's [Division of Genome Diagnostics](https://www.chop.edu/centers-programs/division-genomic-diagnostics) has partnered with CCDI to add somatic panel sequencing data to OpenPedCan and the Molecular Targets Platform.

The OpenPedCan operates on a pull request model to accept contributions from community participants.
The maintainers have set up continuous integration software via GitHub Actions to confirm the reproducibility of analyses within the project’s Docker container.

**New to the project? Please be sure to read the following documentation before contributing:**

1. Learn about the fundamental data used for this project in [**`doc/data-formats.md`**](./doc/data-formats.md) and [**`doc/data-files-description.md`**](./doc/data-files-description.md)
   - A history of data releases can be found in [**`doc/release-notes.md`**](./doc/release-notes.md)
2. See what analyses are being performed in [**`analyses/README.md`**](./analyses/README.md)
3. Read the remainder of this README document in full.
4. Read our contributing guidelines in [**`CONTRIBUTING.md`**](./CONTRIBUTING.md) in full.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Data Description](#data-description)
- [Data Release Process](#data-release-process)
- [How to Obtain OpenPedCan Data](#how-to-obtain-openpedcan-data)
  - [Data Access via Download Script](#data-access-via-download-script)
  - [Data Access via CAVATICA](#data-access-via-cavatica)
- [How to Participate](#how-to-participate)
  - [Planned Analyses](#planned-analyses)
  - [Proposing a New Analysis](#proposing-a-new-analysis)
  - [Implementing an Analysis](#implementing-an-analysis)
    - [Analytical Code and Output](#analytical-code-and-output)
    - [Software Dependencies](#software-dependencies)
    - [Pull Request Model](#pull-request-model)
- [How to Add an Analysis](#how-to-add-an-analysis)
  - [Folder Structure](#folder-structure)
  - [Documenting Your Analysis](#documenting-your-analysis)
  - [Analysis Script Numbering](#analysis-script-numbering)
  - [Output Expectations](#output-expectations)
  - [Docker Image](#docker-image)
    - [Development in the Project Docker Container](#development-in-the-project-docker-container)
      - [RStudio - Local](#rstudio---local)
  - [Development using Amazon EC2 (CHOP)](#development-using-amazon-ec2-chop)
    - [RStudio - EC2](#rstudio---ec2)
  - [Local Development](#local-development)
  - [GitHub Actions (GA)](#github-actions-ga)
    - [Working with the subset files used in GA locally](#working-with-the-subset-files-used-in-ga-locally)
    - [Adding Analyses to Github Actions workflow](#adding-analyses-to-github-actions-workflow)
    - [Adding Analyses with Multiple Steps](#adding-analyses-with-multiple-steps)
      - [1. File and merge a pull request for adding `01-filter-samples.R` to the repository.](#1-file-and-merge-a-pull-request-for-adding-01-filter-samplesr-to-the-repository)
      - [2. File and merge a pull request for adding `02-cluster-heatmap.R` to the repository.](#2-file-and-merge-a-pull-request-for-adding-02-cluster-heatmapr-to-the-repository)
      - [3. File and merge a pull request for the shell script that runs the entirety of `gene-expression-clustering`.](#3-file-and-merge-a-pull-request-for-the-shell-script-that-runs-the-entirety-of-gene-expression-clustering)
    - [Passing variables only in GA](#passing-variables-only-in-ga)
- [Funding](#funding)
- [Contact](#contact)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Data Description

The OpenPedCan dataset includes methylation array, mRNA gene expression, miRNA gene expression, fusion, whole cell proteomics, phospho-proteomics, as well as somatic mutation, copy number, structural and variant results in combined tsv or matrix format.

Below is a summary of biospecimens by sequencing strategy in v14 release:

| Experimental Strategy | Normal | Tumor |
| --------------------- | ------ | ----- |
| Methylation           | 176    | 2606  |
| miRNA-Seq             | 0      | 262   |
| Targeted Sequencing   | 823    | 2261  |
| Phospho-Proteomics    | 0      | 407   |
| Whole Cell Proteomics | 0      | 407   |
| RNA-Seq               | 18115  | 13803 |
| WGS                   | 3491   | 2936  |
| WXS                   | 1288   | 1320  |

[Here](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/analyses/molecular-subtyping-integrate/results/pediatric_cancer_groups.tsv) is a detailed table of pediatric cancer groups in the current release.

## Data Release Process

1. The initial task in the OpenPedcan data release process is to create a `histologies-base.tsv` file via the [D3b Center's histologies-qc](https://github.com/d3b-center/histologies-qc) repository. 
2. A core set of merged matrices from raw genomic algorithm output is created.
The `histologies-base.tsv` file is utilized together with these matrices to run the [pre-release QC](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/data-pre-release-qc) to ensure all samples from the histologies file are represented in the data files, if there is output.
3. The [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) script is run next; this module creates independent specimen lists and analysis output required for molecular subtyping. 
4. The [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh) script is run to subtype all samples and add this information to a final `histologies.tsv` file used in the data release.
5. Once all release files are generated, the [create-subset-files](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/create-subset-files) module is run to create smaller subsets of the data release files for use in continuous integration.

|      Analysis Type     |                                                            Workflow Script                                                            |                                                                  Analysis Module                                                                 | Runtime on EC2 m6i.4xlarge (minutes) |                                                                 Comment                                                                 |
|:----------------------:|:-------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------:|
| Data pre-release QC    |                                                                                                                                       | [data-pre-release-qc](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/data-pre-release-qc)                         | 21.013                               | Performs QC on harmonized pre-release datasets, including primary matrices, before being utilized in all data release analysis modules. |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [copy_number_consensus_call](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/copy_number_consensus_call)           | 67.683                               | Run independently to create the data release copy number consensus file                                                                 |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [independent-samples](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/independent-samples)                         | 00.200                               | Generated using the `histologies-base` file for analysis pre-release and molecular subtyping modules.                                   |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [fusion_filtering](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/fusion_filtering)                               | 81.433                               |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [run-gistic](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/run-gistic)                                           | 56.950                               |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [focal-cn-file-preparation](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/focal-cn-file-preparation)             | 53.217                               |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [fusion-summary](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/fusion-summary)                                   | 00.233                               |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [tmb-calculation](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tmb-calculation)                                 | 13.600                               |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [gene-set-enrichment-analysis](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/gene-set-enrichment-analysis)       | 8.533                                |                                                                                                                                         |
| Analysis pre-release   | [generate-analysis-files.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/generate-analysis-files.sh) | [tp53_nf1_score](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/tp53_nf1_score)                                   | 17.001                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-MB](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-MB)                   | 09.150                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-CRANIO](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-CRANIO)           | 06.083                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-EPN](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-EPN)                 | 07.117                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-embryonal](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-embryonal)     | 05.650                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-chordoma](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-chordoma)       | 01.717                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-EWS](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-EWS)                 | 00.083                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-neurocytoma](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-neurocytoma) | 00.050                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-HGG](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-HGG)                 | 08.383                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-LGAT](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-LGAT)               | 12.517                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-ATRT](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-ATRT)               | 00.067                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-PB](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-PB)               |                                |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-NBL](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-NBL)                 | 03.467                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-pathology](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-pathology)     | 00.567                               |                                                                                                                                         |
| Molecular subtyping    | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [molecular-subtyping-integrate](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/molecular-subtyping-integrate)     | 00.150                               |                                                                                                                                         |
| Analysis pre-release   | [run-for-subtyping.sh](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/scripts/run-for-subtyping.sh)             | [independent-samples](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/independent-samples)                         | 00.412                               | Generated using the `final histologies` file for each cohort and each cohort-cancer_group for downstream analysis modules.              |
| continuous integration |                                                                                                                                       | [create-subset-files](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/create-subset-files)                         | 41.500                               |                                                                                                                                         |

## How to Obtain OpenPedCan Data

We are releasing this dataset on both [CAVATICA](https://cavatica.sbgenomics.com/u/cavatica/opentarget) and AWS S3.
Users performing analyses should always refer to the symlinks in the `data/` directory and not files within the release folder.

**The data formats and caveats are described in more detail in [`doc/data-formats.md`](doc/data-formats.md).
For brief descriptions of the data files, see the [`data-files-description.md`](doc/data-files-description.md) file included in the download.**

Use the [data issue template](https://github.com/d3b-center/ticket-tracker/issues/new?assignees=&labels=data&template=data-question.md&title=) to file issues if you have questions about or identify issues with OpenPedCan data.

### Data Access via Download Script

We have created a shell script that will download the latest release from AWS S3.
macOS users must install `md5sum` before running the download script the first time.
This can be installed with [homebrew](https://brew.sh/) via the command `brew install coreutils` or [conda/miniconda](https://docs.conda.io/projects/conda/en/latest/) via the command `conda install -c conda-forge coreutils`.
_Note: the `download-data.sh` script now has the ability to skip downloads of unchanged files, but if you previously installed md5sum via brew you'll need to run `brew unlink md5sha1sum && brew install coreutils` first to take advantage of this new feature._

Once this has been done, run `bash download-data.sh` to acquire the latest release.
This will create symlinks in `data/` to the latest files.
It's safe to re-run `bash download-data.sh` to check that you have the most recent release of the data.
For downloading DNA methylation-related data, run `bash download-methyl.sh`  in `scripts` folder. 
We will update the default release number whenever we produce a new release.

### Data Access via CAVATICA

For any user registered on CAVATICA, the OpenPBTA and OpenPedcan data can be accessed from the CAVATICA public project below:

- [OpenPBTA Open Access](https://cavatica.sbgenomics.com/u/cavatica/openpbta/)
- [OpenPedcan Open Access](https://cavatica.sbgenomics.com/u/cavatica/opentarget)

The release folder structure in CAVATICA mirrors that on AWS.
Users downloading via CAVATICA should place the data files within the `data/release*` folder and then create symlinks to those files within `/data`.

## How to Participate

### Planned Analyses

There are certain analyses that we have planned or that others have proposed, but which nobody is currently in charge of completing.
Check the existing [issues](https://github.com/d3b-center/OpenPedCan-analysis/issues) to identify these.
If you would like to take on a planned analysis, please comment on the issue noting your interest in tackling the issue in question.
Ask clarifying questions to understand the current scope and goals.
Then propose a potential solution.
If the solution aligns with the goals, we will ask you to go ahead and start to implement the solution.
You should provide updates to your progress in the issue.
When you file a pull request with your solution, you should note that it closes the issue in question.

### Proposing a New Analysis

In addition to the planned analyses, we welcome contributors who wish to propose their own analyses of this dataset as part of the OpenPedCan project.
Check the existing [issues](https://github.com/d3b-center/OpenPedCan-analysis/issues) before proposing an analysis to see if something similar is already planned.
If there is not a similar planned analysis, create a new issue.
The ideal issue will describe the scientific goals of the analysis, the planned methods to address the scientific goals, the input data that is required for the planned methods, and a proposed timeline for the analysis.
Project maintainers will interact on the issue to clarify any questions or raise any potential concerns.

### Implementing an Analysis

This section describes the general workflow for implementing analytical code, and more details are [described below](#how-to-add-an-analysis).
The first step is to identify an existing analysis or propose a new analysis, engage with the project maintainers to clarify the goals of the analysis, and then get the go ahead to move forward with the analysis.

#### Analytical Code and Output

You can perform your analyses via a script (R or Python) or via a notebook (R Markdown or Jupyter).
Your analyses should produce one or more _artifacts_.
Artifacts include both vector or high-resolution figures sufficient for inclusion in a manuscript as well as new summarizations of the data (tables, etc) that are intended for either use in subsequent analyses or distribution with the manuscript.

#### Software Dependencies

Analyses should be performed within the project's [Docker container](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/Dockerfile).
We use a single monolithic container in these analyses for ease of use.
If you need software that is not included, please edit the Dockerfile to install the relevant software or file a [new issue on this repository](https://github.com/d3b-center/OpenPedCan-analysis/issues/new/choose) requesting assistance.

#### Pull Request Model

Analyses are added to this repository via [Pull Requests](https://github.com/d3b-center/OpenPedCan-analysis/compare.
**Please read the [Pull Request section of the contribution guidelines](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/CONTRIBUTING.md#pull-requests) carefully.**
We are using continuous integration software applied to the supplied test datasets to confirm that the analysis can be carried out successfully within the Docker container.

## How to Add an Analysis

Users performing analyses, should **always** refer to the symlinks in the `data/` directory and not files within the release folder, as an updated release may be produced before a publication is prepared.

### Folder Structure

Our folder structure is designed to separate each analysis into its own set of notebooks that are independent of other analyses.
Within the `analyses` directory, create a folder for your analysis.
Choose a name that is unique from other analyses and somewhat detailed.
For example, instead of `gene-expression`, choose `gene-expression-clustering` if you are clustering samples by their gene expression values.
You should assume that any data files are in the `../../data` directory and that their file names match what the `download-data.sh` script produces.
These files should be read in at their relative path, so that we can re-run analyses if the underlying data change.
Files that are primarily graphic should be placed in a `plots` subdirectory and should adhere to the [color palette guide](./figures/README.md#color-palette-usage).
Files that are primarily tabular results files should be placed in a `results` subdirectory.
Intermediate files that are useful within the processing steps but that do not represent final results should be placed in `../../scratch/`.
It is safe to assume that files placed in `../../scratch` will be available to all analyses within the same folder.
It is not safe to assume that files placed in `../../scratch` will be available from analyses in a different folder.

An example highlighting a `new-analysis` directory is shown below.
The directory is placed alongside existing analyses within the `analyses` directory.
In this case, the author of the analysis has run their workflows in R Markdown notebooks.
This is denoted with the `.Rmd` suffix.
However, the author could have used Jupyter notebooks, R scripts, or another scriptable solution.
The author has produced their output figures as `.pdf` files.
We have a preference for vector graphics as PDF files, though other forms of vector graphics are also appropriate.
The results folder contains a tabular summary as a comma separated values file.
We expect that the file suffix (`.csv`, `.tsv`) accurately denotes the format of the added files.
The author has also included a `README.md` ([see Documenting Your Analysis](#documenting-your-analysis)).

```
OpenPedCan-analysis
├── CONTRIBUTING.md
├── README.md
├── analyses
│   ├── existing-analysis-1
│   └── new-analysis
│       ├── 01-preprocess-data.Rmd
│       ├── 02-run-analyses.Rmd
│       ├── 03-make-figures.Rmd
│       ├── README.md
│       ├── plots
│       │   ├── figure1.pdf
│       │   └── figure2.pdf
│       ├── results
│       │   └── tabular_summary.csv
│       └── run-new-analysis.sh
├── data
└── scratch
```

### Documenting Your Analysis

A goal of the OpenPedCan project is to create a collection of workflows that are commonly used for atlas papers.
As such, documenting your analytical code via comments and including information summarizing the purpose of your analysis is important.

When you file the first pull request creating a new analysis module, add your module to the [Modules At A Glance table](analyses#modules-at-a-glance).
This table contains fields for the directory name, what input files are required, a short description, and any files that you expect other analyses will rely on.
As your analysis develops and input or output files change, please check this table remains up to date.
This step is included in the pull request reproducibility checklist.

When an analysis module contains multiple steps or is nearing completion, add a `README.md` file that summarizes the purpose of the module, any known limitations or required updates, and includes examples for how to run the analyses to the folder.

### Analysis Script Numbering

As shown above, analysis scripts within a folder should be numbered from `01` and are intended be run in order.
If the script produces any intermediate files, these files should be placed in `../../scratch`, which is used as described above.
A shell script that runs all analytical code in the intended order should be added to the analysis directory (e.g. `run-new-analysis.sh` above).
See the [continuous integration instructions for adding analyses with multiple steps](#adding-analyses-with-multiple-steps) for more information.

### Output Expectations

The GA system that we use will generate, as artifacts, the contents of the `analyses` directory applied over a small test dataset.
Our goal is to capture all of the outputs that produce final results as artifacts.
Files that are primarily graphic should be placed in a `plots` subdirectory of the analysis's folder.
Plots should use the specified color palettes for this project.
See more [specific instructions on how to use the color palette here](./figures/README.md#color-palette-usage).
Files that are primarily tabular results files should be placed in a `results` subdirectory of the analysis's folder.
Files that are intermediate, which means that they are useful within an analysis but do not provide final outputs should be placed in `../../scratch`.

### Docker Image

We build our project Docker image from a versioned [`tidyverse`](https://hub.docker.com/r/rocker/tidyverse) image from the [Rocker Project](https://www.rocker-project.org/) (v4.2.3).

To add dependencies that are required for your analysis to the project Docker image, you must alter the project [`Dockerfile`](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/Dockerfile).
The `Dockerfile` can be directly edited to install dependencies, if you are developing using a branch on the [d3b-center/OpenPedCan-analysis](https://github.com/d3b-center/OpenPedCan-analysis) repository.
If you are developing using a branch on your fork of the d3b-center/OpenPedCan-analysis repository, create a branch on the d3b-center/OpenPedCan-analysis repository to edit the `Dockerfile` to install dependencies, e.g. <https://github.com/d3b-center/OpenPedCan-analysis/pull/36>, so [the GitHub action for checking docker image build](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/.github/workflows/build-docker.yml) can run with the Docker Hub credentials saved in the d3b-center/OpenPedCan-analysis repository.

- R packages installed on this image will be installed from an [MRAN snapshot](https://mran.microsoft.com/documents/rro/reproducibility#reproducibility) corresponding to the last day that R 3.6.0 was the most recent release ([ref](https://hub.docker.com/r/rocker/tidyverse)).
  - Installing most packages, from CRAN or Bioconductor, should be done with our `install_bioc.R` script, which will ensure that the proper MRAN snapshot is used. `BiocManager::install()` should _not_ be used, as it will not install from MRAN.
  - R packages that are not available in the MRAN snapshot can be installed via github with the `remotes::install_github()` function, with the commit specified by the `ref` argument.
- Python packages should be installed with `pip3 install` with version numbers for all packages and dependencies specified.
  - As a secondary check, we maintain a `requirements.txt` file to check versions of all python packages and dependencies.
  - When adding a new package, make sure that all dependencies are also added; every package should appear with a specified version **both** in the `Dockerfile` and `requirements.txt`.
- Other software can be installed with `apt-get`, but this should _never_ be used for R packages.

If you need assistance adding a dependency to the Dockerfile, [file a new issue on this repository](https://github.com/d3b-center/ticket-tracker/issues/new) to request help.

#### Development in the Project Docker Container

If you are new user download Docker from [here](https://docs.docker.com/get-docker/)
You may need to run each of the docker commands with the platform (eg: `--platform linux/amd64`).
The most recent version of the project Docker image, which is pushed to Docker Hub after a pull request gets merged into the dev branch, can be obtained via the command line with:

```bash
docker pull pgc-images.sbgenomics.com/d3b-bixu/openpedcanverse:latest
```

Development should utilize the project Docker image.
An analysis that is developed using the project Docker image can be efficiently rerun by another developer or the original developer (after a long time since it is developed), without dependency or numerical issues.
This will significantly facilitate the following tasks that are constantly performed by all developers of the OpenPedCan-analysis project.

- Review another developer's pull request, including code and results. For more information about pull request and review, see [the guideline for how to contribute to the OpenPedCan-analysis repository](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/CONTRIBUTING.md#contribution-guidelines-for-the-openpbta-analysis).
- Update the results of an analysis module that is developed by another developer. For example, rerun the same analysis module with new data.
- Update the code of an analysis module that is developed by another developer. For example, add a new feature to a module, or refactor a module.

**If you are a Mac or Windows user, the default limit for memory available to Docker is 2 GB.
You will likely need to increase this limit for local development.**
[[Mac documentation](https://docs.docker.com/docker-for-mac/#resources), [Windows documentation](https://docs.docker.com/docker-for-windows/#advanced)]

##### RStudio - Local

Using `rocker/tidyverse:3.6.0` as our base image allows for development via RStudio in the project Docker container.
If you'd like to develop in this manner, you may do so by running the following and changing `<password>` to a password of you choosing at the command line:

```bash
docker run -e PASSWORD=<password> -p 8787:8787 pgc-images.sbgenomics.com/d3b-bixu/openpedcanverse:latest
```

You can change the volume that the Docker container points to either via the [Kitematic GUI](https://docs.docker.com/kitematic/userguide/) or the [`--volume` flag](https://docs.docker.com/storage/volumes/) to `docker run`.
For example, from your cloned `OpenPedCan-analysis` folder, run the command below:

```bash
docker run --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/OpenPedCan-analysis pgc-images.sbgenomics.com/d3b-bixu/openpedcanverse:latest
```

Once you've set the volume, you can navigate to `localhost:8787` in your browser if you are a Linux or Mac OS X user.
The username will for login will be `rstudio` and the password will be whatever password you set with the `docker run` command above.

If you are a new user, you may find [these instructions](https://github.com/AlexsLemonade/RNA-Seq-Exercises/blob/master/docker-pull.md) for a setting up a different Docker [container](https://www.andrewheiss.com/blog/2017/04/27/super-basic-practical-guide-to-docker-and-rstudio/) or [this guide](https://www.andrewheiss.com/blog/2017/04/27/super-basic-practical-guide-to-docker-and-rstudio/) from Andrew Heiss helpful.

You can also run the analysis on the terminal once you have a docker container running locally by running `docker exec` helpful information [here](https://buildvirtual.net/docker-exec-what-does-it-do/)

```
docker exec -ti <CONTAINER_NAME> bash
```

If you set the `PWD:/home/rstudio/OpenPedCan-analysis` above, then you can navigate to `home/rstudio/OpenPedCan-analysis` and begin.

### Development using Amazon EC2 (CHOP)

Many analyses will require Amazon EC2 for development.
For this, we have created a template image in `Mgmt-Console-Dev-chopd3bprod`.
Navigate to the Service Catalog and select `d3b-research-instance`.
The standard mount comes with a default 100 GB root volume.
Below are the instance names, hourly rates, vCPUs, and memory.

| Instance name | Hourly rate | vCPU | Memory |
| ------------- | ----------- | ---- | ------ |
| m6i.large     | $0.096      | 2    | 8 GB   |
| m6i.xlarge    | $0.192      | 4    | 16 GB  |
| m6i.2xlarge   | $0.384      | 8    | 32 GB  |
| m6i.4xlarge   | $0.768      | 16   | 64 GB  |
| m6i.8xlarge   | $1.536      | 32   | 128 GB  |

#### RStudio - EC2

To use RStudio from EC2, run docker using the following command:

```bash
docker run --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 80:8787 -v $PWD:/home/rstudio/OpenPedCan-analysis pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:latest
```

Then, paste the instance IP address into your browser to start RStudio.

### Local Development

While we encourage development within the Docker container, it is also possible to conduct analysis without Docker if that is desired.
In this case, it is important to ensure that local or personal settings such as file paths or installed packages and libraries are not assumed in the analysis.

### GitHub Actions (GA)

We use GitHub Actions (GA) to ensure that the project Docker image will build if there are any changes introduced to the [`Dockerfile`](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/Dockerfile) and that all analysis code will execute.

We have put together data files specifically for the purpose of GA that contain all of the features of the full data files for only a small subset of samples.
You can see how this was done by viewing [this module](https://github.com/d3b-center/OpenPedCan-analysis/tree/dev/analyses/create-subset-files).
We use the subset files to cut down on the computational resources and time required for testing.

Provided that your analytical code points to the symlinks in the `data/` directory per [the instructions above](#how-to-add-an-analysis), adding the analysis to the GA (see below) will run your analysis on this subset of the data.
Do not hardcode sample names in your analytical code: there is no guarantee that those samples will be present in the subset files.

#### Working with the subset files used in GA locally

If you would like to work with the files used in GA locally, e.g., for debugging, you can obtain them from AWS by running the following in the root directory of the project:

```
bash scripts/download-testing-files.sh
```

Running this will change the symlinks in `data` to point to the files in `data/testing`.

#### Adding Analyses to Github Actions workflow

For an analysis to be run in a Github Actions workflow, it must be added to [`.github/continuous_integration.yml`](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/.github/workflows/continuous_integration.yml). Here is an example of a step in that workflow:

```yaml
# Molecular subtyping modules
- name: Molecular Subtyping - MB
  entrypoint: molecular-subtyping-MB/run-molecular-subtyping-mb.sh
  openpbta_subset: 0
```

Each analysis entrypoint will be run with the container image built in the earlier stage in the workflow, and must be present in the `analyses/` directory in the repository. The section below shows how each entrypoint will be invoked in the actions workflow:

```yaml
- name: Run Analysis
  uses: docker://pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:analysisjob
  with:
    entrypoint: analyses/${{ matrix.entrypoint }}
  env:
    OPENPBTA_SUBSET: ${{ matrix.openpbta_subset }}
    OPENPBTA_TESTING: ${{ matrix.openpbta_testing }}
    RUN_FOR_SUBTYPING: ${{ matrix.run_for_subtyping }}
    OPENPEDCAN_POLYA_STRAND: ${{ matrix.openpedcan_polya_strand }}
```

In this workflow, a new container image is built and pushed to the SBG DockerHub under tag `pgc-images.sbgenomics.com/d3b-bixu/open-pedcan:analysisjob`, then pulled and each analysis entrypoint is run in parallel with a Github Actions matrix build.

Because of the dependency on the image in DockerHub, the following changes will need to be made on any fork of this repository before running this job:

1. The registry under the `registry:` parameter in the step below will need to be changed to a registry that you control.

```yaml
- name: Login to DockerHub
  uses: docker/login-action@v2
  with:
    registry: pgc-images.sbgenomics.com
    username: ${{ secrets.DOCKER_HUB_USERNAME }}
    password: ${{ secrets.DOCKER_HUB_ACCESS_TOKEN }}
```

2. The image name under `images:` will need to be changes to fit the registry from step 1.

```yaml
- name: Docker meta
  id: meta
  uses: docker/metadata-action@v4
  with:
    images: pgc-images.sbgenomics.com/d3b-bixu/open-pedcan
    tags: |
      type=raw,value=analysisjob
      # Only tag the image with latest if we're building on the default
      # branch (e.g., dev).
      type=raw,value=latest,enable={{is_default_branch}}
```

3. A Username and Access Token for the new Docker Hub registry will need to be added in [secrets](https://docs.github.com/en/actions/security-guides/encrypted-secrets#creating-encrypted-secrets-for-a-repository) under the forked repository will need to be added.

To add a new analysis job, take the template below and value each missing prompt, then add it to [`.github/continuous_integration.yml`](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/.github/continuous_integration.yml).

```yaml
- name: <Name of the analysis to run>
  entrypoint: <Path to the analysis script from the analyses/ directory>
```

Optionally, environment variables for `OPENPBTA_SUBSET`, `OPENPBTA_TESTING`, `RUN_FOR_SUBTYPING`, and `OPENPEDCAN_POLYA_STRAND` can be passed in using the syntax below.

```yaml
- name: <Name of the analysis to run>
  entrypoint: <Path to the analysis script from the analyses/ directory>
  openpbta_subset: <Value for OPENPBTA_SUBSET>
  openpbta_testing: <Value for OPENPBTA_TESTING>
  run_for_subtyping: <Value for RUN_FOR_SUBTYPING>
  openpedcan_polya_strand: <Value for OPENPEDCAN_POLYA_STRAND>
```

#### Adding Analyses with Multiple Steps

There is a different procedure for adding an analysis comprised of multiple scripts or notebooks to GA.
Per [the contribution guidelines](https://github.com/d3b-center/OpenPedCan-analysis/blob/dev/CONTRIBUTING.md#size-and-composition-of-pull-requests), each script or notebook should be added via a separate pull request.
For each of these pull requests, the individual script or notebook should be added as its own run in the `.github/continuous_integration.yml` file.
This validates that the code being added can be executed at the time of review.

Once all code for an analysis has been reviewed and merged, a final pull request for the analysis that is comprised of the following changes should be filed:

- A shell script that will run all script and/or notebooks in the analysis module.
- The multiple runs from the module that are in the `continuous_integration.yml` file are replaced with a single run that runs the shell script.

If the `gene-expression-clustering` analysis instead required two scripts run sequentially (`01-filter-samples.R` and `02-cluster-heatmap.R`), we would follow the procedure below.

##### 1. File and merge a pull request for adding `01-filter-samples.R` to the repository.

In this pull request, we would add the following change to `.github/continuous_integration.yml`.

```yaml
- name: Run Filter Samples
  entrypoint: gene-expression-clustering/01-filter-samples.sh
```

##### 2. File and merge a pull request for adding `02-cluster-heatmap.R` to the repository.

In this pull request, we would add the following change to `.github/continuous_integration.yml`.
This would be added _below_ the `Filter Samples` run.

```yaml
- name: Cluster samples and plot heatmap
  entrypoint: gene-expression-clustering/02-cluster-heatmap.sh
```

##### 3. File and merge a pull request for the shell script that runs the entirety of `gene-expression-clustering`.

In this pull request, we would add a shell script that runs `01-filter-samples.R` and `02-cluster-heatmap.R`.
Let's call this shell script `run-gene-expression-clustering.sh` and place it in the analysis directory `analyses/gene-expression-clustering`.

The contents of `analyses/gene-expression-clustering/run-gene-expression-clustering.sh` may look like:

```
#!/bin/bash
# This script runs the gene-expression-clustering analysis
# Author's Name 2019

set -e
set -o pipefail

Rscript --vanilla analyses/gene-expression-clustering/01-filter-samples.R
Rscript --vanilla analyses/gene-expression-clustering/02-cluster-heatmap.R
```

We would remove the runs `Filter Samples` and `Cluster Samples and Plot Heatmap` from `.github/continuous_integration.yml` and instead replace them with a single run:

```yaml
- name: Cluster samples and plot heatmap
  entrypoint: gene-expression-clustering/run-gene-expression-clustering.sh
```

#### Passing variables only in GA

The analyses run in GA use only a small portion of the data so that tests can be run quickly.
For some analyses, there will not be enough samples to fully test code without altering certain parameters passed to methods.
The preferred way to handle these is to run these analyses through a shell script that specifies default parameters using environment variables.
The default parameters should be the ones that are most appropriate for the full set of data.
In GA, these will be replaced.

We might decide that it makes the most sense to run an analysis using a more permissive statistical significance threshold in GA so that some "significant" pathways still appear and subsequent code that examines them can be tested.
We'd first write code capable of taking command line parameters.
In R, we could use `optparse` to specify these in a script - imagine it's called `pathway_sig.R` and it contains an option list:

```R
option_list <- list(
  optparse::make_option(
    c("-a", "--alpha"),
    type = "double",
    help = "pathway significance threshold",
  )
)
```

Then we would create a shell script (perhaps `run_pathway_sig.sh`) that uses a default environment variable. If `OPENPBTA_PATHSIG` is defined, it will be used. Otherwise, a value of 0.05 is used.
Note: the `-` before the `0.05` below is necessary notation for a default parameter and _not_ designating a negative 0.05.

```bash
PATHSIG=${OPENPBTA_PATHSIG:-0.05}

Rscript analyses/my-path/pathway_sig.R --alpha $PATHSIG
```

We can override this by passing environment variables in `.github/continuous_integration.yml`.
For testing, we might want to use an alpha level of 0.75 so that at least some "significant" pathways appear, which allows testing subsequent code that depends on them.
The name command in the `.github/continuous_integration.yml` is used to specify these parameters.

```yaml
- name: Cluster samples and plot heatmap
  entrypoint: OPENPBTA_PATHSIG=0.75 my-path/run_pathway_sig.sh
```

In this example `OPENPBTA_PATHSIG=0.75` species an environment variable `OPENPBTA_PATHSIG` that is set to 0.75.
Any environment variables prefixed with `OPENPBTA_` are passed to the specified shell script.
Environment variables without this prefix are not passed.


## Funding
[**OpenPBTA**](https://github.com/AlexsLemonade/OpenPBTA-analysis) was funded through the [Children's Brain Tumor Network (CBTN)](https://cbtn.org/) by the following donors who provided leadership level support: CBTN Executive Council members, Brain Tumor Board of Visitors, Children's Brain Tumor Foundation, Easie Family Foundation, Kortney Rose Foundation, Lilabean Foundation, Minnick Family Charitable Fund, Perricelli Family, Psalm 103 Foundation, and Swifty Foundation. 
Additional funding was provided by Alex’s Lemonade Stand Foundation (ALSF) Childhood Cancer Data Lab, ALSF Young Investigator Award, ALSF Catalyst Award, ALSF Catalyst Award, ALSF CCDL Postdoctoral Training Grant, Children’s Hospital of Philadelphia Division of Neurosurgery, Australian Government, Department of Education, St. Anna Kinderkrebsforschung, Austria, the Mildred Scheel Early Career Center Dresden P2, funded by the German Cancer Aid, NIH Grants 3P30 CA016520-44S5, U2C HL138346-03, U24 CA220457-03, K12GM081259, R03-CA23036, NIH Contract Nos. HHSN261200800001E and 75N91019D00024, Task Order No. 75N91020F00003, and the Intramural Research Program of the Division of Cancer Epidemiology and Genetics of the National Cancer Institute.

Inaugural funding for **OpenPedCan** was provided in part by [NCI's Childhood Cancer Data Initiative](https://www.cancer.gov/research/areas/childhood/childhood-cancer-data-initiative) through NIH Task Order No. 75N91020F00003, the CBTN, and the Children’s Hospital of Philadelphia Division of Neurosurgery.


## Contact
Please submit an issue with any questions, bugs, or feature requests.





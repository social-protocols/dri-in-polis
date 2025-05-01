# Deliberative Reason Index for Pol.is

This repository contains a Julia reproduction of the Deliberative Reason Index (DRI) calculation, as well as some experiments I have run related to the DRI calculation and its possible integration with [pol.is](https://pol.is/)

The DRI calculation is based on the R implementation from the supplementary materials to the paper [*How Deliberation Happens: Enabling Deliberative Reason*](https://www.cambridge.org/core/journals/american-political-science-review/article/how-deliberation-happens-enabling-deliberative-reason/6558F69855ADA8B15BF2EC2E5D403E71).

The output of the experiments and a brief writeup of the methodology of each experiment are below.

## Experiments

- [**DRI Reference Implementation**](docs/reference-implementation.md): Re-produces the calculations and some charts from the above paper.
- [**Proof-of-Concept DRI in pol.is**](docs/poc-polis.md): Experiment with a variant of the DRI that uses pol.is data.
- [**Random-Tagging Analysis**](docs/random-tagging.md): Experiment with randomly classifying statements as preferences/considerations to see if it makes a difference.
- [**Statement Subset Selection Analysis**](docs/statement-subset.md): Experiment with using only a subset of statements -- the most or least controversial. 
- [**Dual PCA**](docs/dual-PCA.md): Experiment on DRI using reduced-dimension data. 
- [**Validity Test**](docs/validity-test.md): A validity test showing positive DRI can be produced from random data.

## Background and Attribution

The original DRI calculation is described in:

**NIEMEYER, SIMON, FRANCESCO VERI, JOHN S. DRYZEK, and ANDRÉ BÄCHTIGER. 2024. “How Deliberation Happens: Enabling Deliberative Reason.” American Political Science Review 118(1): 345–62. doi: 10.1017/S0003055423000023.**

**Important Note:**  
The supplementary materials (including the original R code and input data files) are hosted on [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY) and are governed by terms that restrict their redistribution. This repository solely provides an independent Julia port for calculating the DRI.

## Installation

1. **Install Julia:**  
   Download and install Julia from [julialang.org](https://julialang.org/downloads/).

2. **Clone the Repository:**  
   ```bash
   git clone https://github.com/social-protocols/dri-in-polis.git
   cd dri-julia
   ```

3. **Install Dependencies:**

    julia --project -e 'using Pkg; Pkg.instantiate()'

## Downloading Data

### Downloading Data for Reference Implementation

- Go to https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY
- Click the download icon next to Data1_Raw_Input.tab
- Choose "Comma Separated Values (Original Format)"
- Save file under Input/ 


### Downloading Data for Pol.is Implementation

Run:

    git submodule update --init --recursive


## Running

### Running Reference Implementation

Run:

    julia --project reference-implementation.jl

This will produce data and plots in the directory local-output/reference-implementation.

### Running Pol.is Implementation

Run:

    julia --project poc-polis.jl

This will produce data and plots in the directory local-output/poc-polis.

By default, this runs the analysis for the vtaiwan.uberx case. You can optionally specify any case under polis/openData/, which is a clone of https://github.com/compdemocracy/openData, as the first argument.

The second argument can be the name of a correlation method, currently "phi", "pearson", or "pearson_binary"

    julia --project poc-polis.jl american-assembly.bowling-green pearson_binary

This will produce data and plots in the directory local-output/$CASE.

### Running Random Tagging Analysis

Run:

    julia --project random-tagging.jl

or to specific a number of samples and a method: 

    julia --project random-tagging.jl 10000 method2
 
This will produce a .png for each case with a histogram of the random-tag DRI values (pre and post) along with a p-value under local-output/. P-values are also output to a .CSV file in the same directory.

### Running Statement Subset Analysis

Run:

    julia --project statement-subset.jl

This will place output into local-output/statement-subset

### Running DRI Validity Check

Run:

    julia --project validity-check.jl

This will place output into local-output/statement-subset


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


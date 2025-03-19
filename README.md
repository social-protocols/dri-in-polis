# Deliberative Reason Index for Pol.is

This repository contains a Julia implementation of the Deliberative Reason Index (DRI) modified to work with [pol.is](https://pol.is/) data.

The DRI calculation is based on the R implementation from the supplementary materials to the paper *How Deliberation Happens: Enabling Deliberative Reason*.

## Background and Attribution

The original DRI calculation is described in:

**NIEMEYER, SIMON, FRANCESCO VERI, JOHN S. DRYZEK, and ANDRÉ BÄCHTIGER. 2024. “How Deliberation Happens: Enabling Deliberative Reason.” American Political Science Review 118(1): 345–62. doi: 10.1017/S0003055423000023.**

**Important Note:**  
The supplementary materials (including the original R code and input data files) are hosted on [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY) and are governed by terms that restrict their redistribution. This repository solely provides an independent Julia port for calculating the DRI.

## Reference Implementation

This will produce similar data and plots to those from the research paper. The DRI calculations are identical. Chart styles are different, and exact data point positions in the charts differ slightly because of random jitter (added in reference implementation).

## Pol.is implementation

Here we use pol.is data as the input. Because user votes in pol.is are binary agree/disagree (or pass) values, we use a different methods for calculating the correlation between two users values. By default, we use the Phi coefficient, but Pearson's or other correlation coefficients can be plugged in. 

NOTE: as of this initial implementation, statements are classified RANDOMLY as considerations/preferences.

## Installation

1. **Install Julia:**  
   Download and install Julia from [julialang.org](https://julialang.org/downloads/).

2. **Clone the Repository:**  
   ```bash
   git clone https://github.com/johnwarden/dri-in-polis.git
   cd dri-julia
   ```

## Downloading Data for Reference Implementation

- Go to https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY
- Click the download icon next to Data1_Raw_Input.tab
- Choose "Comma Separated Values (Original Format)"
- Save file under Input/ 

## Running Reference Implementation

   julia --project -e 'using Pkg; Pkg.instantiate()'
   julia --project reference_dri.jl

This will produce data and plots in the directory Output/reference.

## Running Pol.is Implementation

   julia --project -e 'using Pkg; Pkg.instantiate()'
   julia --project polis_dri.jl

### Arguments

By default, this runs the analysis for the vtaiwan.uberx case. You can optionally specify any case under polis/openData/, which is a clone of https://github.com/compdemocracy/openData, as the first argument.

The second argument can be the name of a correlation method, currently "phi", "pearson", or "pearson_binary"


   julia --project polis_dri.jl american-assembly.bowling-green pearson_binary

This will produce data and plots in the directory Output/$CASE.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.





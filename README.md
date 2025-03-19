# Deliberative Reason Index for Pol.is

This repository contains a Julia implementation of the Deliberative Reason Index (DRI) modified to work with [pol.is](https://pol.is/) data.

The DRI calculation is based on the R implementation from the supplementary materials to the paper *How Deliberation Happens: Enabling Deliberative Reason*.

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


## Downloading Data for Reference Implementation

- Go to https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY
- Click the download icon next to Data1_Raw_Input.tab
- Choose "Comma Separated Values (Original Format)"
- Save file under Input/ 

## Reference Implementation

   julia --project reference_dri.jl

This will produce data and plots in the directory Output/reference.

These are very similar to those from the research paper. The DRI calculations are identical. Chart styles are different, and exact data point positions in the charts differ slightly because of random jitter (added in reference implementation).

## Pol.is Implementation

   julia --project polis_dri.jl

Here we use pol.is data as the input. Because user votes in pol.is are binary agree/disagree (or pass) values, we use a different methods for calculating the correlation between two users values. By default, we use the Phi coefficient, but Pearson's or other correlation coefficients can be plugged in. 

NOTE: as of this initial implementation, statements are classified RANDOMLY as considerations/preferences.


### Arguments

By default, this runs the analysis for the vtaiwan.uberx case. You can optionally specify any case under polis/openData/, which is a clone of https://github.com/compdemocracy/openData, as the first argument.

The second argument can be the name of a correlation method, currently "phi", "pearson", or "pearson_binary"

   julia --project polis_dri.jl american-assembly.bowling-green pearson_binary

This will produce data and plots in the directory Output/$CASE.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Sample Output Charts

Below are sample output charts generated from the reference implementation and the pol.is implementations. These files are hosted on GitHub and can be viewed directly in your browser.

### Reference Implementation

Some sample figures produced by the reference implementation:

- **Fig2_a_Pre:**  
  ![Fig2_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig2_a_Pre.png)
- **Fig2_b_Post:**  
  ![Fig2_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig2_b_Post.png)
- **Fig3_1Control_a_Pre:**  
  ![Fig3_1Control_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_1Control_a_Pre.png)
- **Fig3_1Control_b_Post:**  
  ![Fig3_1Control_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_1Control_b_Post.png)
- **Fig3_2Brief_a_Pre:**  
  ![Fig3_2Brief_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_2Brief_a_Pre.png)
- **Fig3_2Brief_b_Post:**  
  ![Fig3_2Brief_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_2Brief_b_Post.png)
- **Fig3_3Building_a_Pre:**  
  ![Fig3_3Building_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_3Building_a_Pre.png)
- **Fig3_3Building_b_Post:**  
  ![Fig3_3Building_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_3Building_b_Post.png)

### Pol.is Implementation

The pol.is implementation produces PNG charts. Below are sample outputs for two cases:

#### vtaiwan.uberx Case

- **Correlation (Phi):**  
  ![vtaiwan.uberx - Phi](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/vtaiwan.uberx/Figures/polis_dri_plot-phi.png)

- **Correlation (Pearson Binary):**  
  ![vtaiwan.uberx - Pearson Binary](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/vtaiwan.uberx/Figures/polis_dri_plot-pearson_binary.png)

#### london.youth.policing Case

- **Correlation (Phi):**  
  ![london.youth.policing - Phi](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/london.youth.policing/Figures/polis_dri_plot-phi.png)

- **Correlation (Pearson Binary):**  
  ![london.youth.policing - Pearson Binary](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/london.youth.policing/Figures/polis_dri_plot-pearson_binary.png)

  ![london.youth.policing - Phi](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/london.youth.policing/Figures/polis_dri_plot-phi




# Deliberative Reason Index for Pol.is

This repository contains a Julia implementation of the Deliberative Reason Index (DRI) modified to work with [pol.is](https://pol.is/) data.

The DRI calculation is based on the R implementation from the supplementary materials to the paper [*How Deliberation Happens: Enabling Deliberative Reason*](https://www.cambridge.org/core/journals/american-political-science-review/article/how-deliberation-happens-enabling-deliberative-reason/6558F69855ADA8B15BF2EC2E5D403E71).

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



## Reference Implementation

The reference implementation re-produces the calculations and some of the charts from the research paper. The DRI calculations are identical. Chart styles are different, and exact data point positions in the charts differ slightly because of random jitter (added in reference implementation). See [charts](#reference-implementation-charts) below.

### Downloading Data for Reference Implementation

- Go to https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6OKNOY
- Click the download icon next to Data1_Raw_Input.tab
- Choose "Comma Separated Values (Original Format)"
- Save file under Input/ 

### Running Reference Implementation

Run:

    julia --project reference_dri.jl

This will produce data and plots in the directory Output/reference.

### Output

And CSV files are [here](https://github.com/social-protocols/dri-in-polis/tree/master/SampleOutput/reference)

**Charts:**

![Fig2_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig2_a_Pre.png)
![Fig2_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig2_b_Post.png)
![Fig3_1Control_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_1Control_a_Pre.png)
![Fig3_1Control_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_1Control_b_Post.png)
![Fig3_2Brief_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_2Brief_a_Pre.png)
![Fig3_2Brief_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_2Brief_b_Post.png)
![Fig3_3Building_a_Pre](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_3Building_a_Pre.png)
![Fig3_3Building_b_Post](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/reference/Figures/Fig3_3Building_b_Post.png)


## Pol.is Implementation

Here we use pol.is data as the input. Because user votes in pol.is are binary agree/disagree (or pass) values, we use a different methods for calculating the correlation between two users values. By default, we use the Phi coefficient, but Pearson's or other correlation coefficients can be plugged in. 

**NOTE**: currently, statements are classified RANDOMLY as considerations/preferences.

### Downloading Data for Pol.is Implementation

Run:

    git submodule update --init --recursive

### Running Pol.is Implementation

Run:

    julia --project polis_dri.jl

This will produce data and plots in the directory Output/polis.

By default, this runs the analysis for the vtaiwan.uberx case. You can optionally specify any case under polis/openData/, which is a clone of https://github.com/compdemocracy/openData, as the first argument.

The second argument can be the name of a correlation method, currently "phi", "pearson", or "pearson_binary"

    julia --project polis_dri.jl american-assembly.bowling-green pearson_binary

This will produce data and plots in the directory Output/$CASE.

### Output Charts

Below are sample outputs for two pol.is cases:

#### vtaiwan.uberx Case

- **Correlation (Phi):**  
  ![vtaiwan.uberx - Phi](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/vtaiwan.uberx/Figures/polis_dri_plot-phi.png)


#### london.youth.policing Case

- **Correlation (Phi):**  
  ![london.youth.policing - Phi](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/polis/london.youth.policing/Figures/polis_dri_plot-phi.png)



## DRI Resampling Analysis 

The pol.is implementation shows us that that even designating statements as considerations or preferences *randomly* results in a positive DRI. This is because, generally, the intersubjective-agreement (pairwise Spearman correlation) on a random group of statements will correlate with intersubjective agreement on another random group of statements.

This naturally raises the question of whether inter-subjective agreement on considerations correlates with inter-subjective agreement on preferences to a *greater degree* than it does to a random set of statements.


A simple way to test this hypothesis is to ignore the way the statements are tagged in the dataset, and instead tag them randomly, and see if we get a similar result. If we repeat this process n times, we can count the number of times the resulting DRI is higher than the actual DRI. This gives us an estimate of the probability that the actual DRI would be as high as it is if tags were random.

resampling_dri.jl tests this hypothesis using the data from all the cases in the paper. It randomly tags the statements and calculates the DRI, repeating this process n times for each case, then produced a histogram and a p-value: which is the probability that random tagging would result in a DRI as extreme.

### Running Resampling Analysis

Run:

    julia --project resampling_dri.jl

This will produce a .png for each case with a histogram of the random-tag DRI values (pre and post) along with a p-value under Output/. P-values are also output to a .CSV file in the same directory.

### Results

The aggregate p-value (using Fisher's method) is 0.6616938467209219

### Alternative Methodology

One subtle problem with this methodology: the dataset contains users *ranks* of each consideration/preferences. If we change the set of items counted as a consideration or preference, the rankings lose their meaning.

To account for this, we implemented a slightly modified methodology. First, calculate the DRI based on a *random subset of considerations* against *all preferences*. This should give us a value that *approximates* the actual DRI (I have confirmed that the average of this value over many samples is approximately equal to the actual DRI). This way, we are no longer calculating Spearman correlations based on a mixed set of considerations/preferences.

Next, we calculate a hypothetical DRI using the same subset of considerations against all other considerations.

Repeat this n times and calculate the number of times the resampled DRI is greater than the control DRI to obtain a p-value.

### Running Resampling Analysis with Alternative Methodology

Run:

   julia --project resampling_dri.jl 10000 method2


### Results

The aggregate p-value (using Fisher's method) is 0.3043455661723067.


## Output

Charts for all cases can be found [here](https://github.com/social-protocols/dri-in-polis/tree/master/SampleOutput/resampling)

Below are some sample output charts:

![Case 18: GBR Futures](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/resampling/dri_distribution_case18.0_method2.png)

![Case 1: Biobanking WA](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/SampleOutput/resampling/dri_distribution_case11.0_method1.png)




## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


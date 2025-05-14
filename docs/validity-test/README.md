# **Technical Note: A Validity Test of the Deliberative Reason Index Using Randomized Input**  
_Jonathan Warden_  

---

## 1. Purpose & Scope  

The **Deliberative Reason Index (DRI)**, as proposed by Niemeyer, Veri, Dryzek, and Bächtiger (2024), aims to quantify how closely participants’ agreement on policy preferences aligns with their agreement on underlying beliefs and values. High DRI scores (near +1) are taken to signal effective, mutually understood reasoning or **intersubjective consistency**. In this technical note we show—via three brief simulation exercises—that **high DRI values can arise even when there is no intersubjective consistency**. Our goal is to provide a concise, reproducible demonstration (with [code posted here](https://github.com/social-protocols/dri-in-polis)) that invites further discussion and refinement of the index.

## 2. What DRI Measures  

- **Conceptual aim:** DRI is designed to capture intersubjective consistency as the degree to which participants who similarly rank a set of “consideration” statements also similarly rank preferences, and vice versa.  
- **Why it matters:** Unlike measures of mere preference convergence or procedural quality, DRI focuses on the coherence between *reasons* and *conclusions*, a core ideal in deliberative‐democratic theory.

## 3. Brief Method: Calculating DRI  

1. **Survey data**: Each participant provides (a) an ordinal ranking of ~20–40 consideration statements and (b) a ranking of ≤10 policy options.  
2. **Pairwise correlations**: For every pair (i, j), compute Spearman’s ρ for their consideration-rankings (ρ<sub>C</sub>) and for their preference-rankings (ρ<sub>P</sub>).  
3. **Consistency distance**: Plot (ρ<sub>C</sub>, ρ<sub>P</sub>) on an x–y scatter; measure each point’s orthogonal distance to the 45° line (y = x). Small distances indicate high deliberative consistency.  
4. **Aggregate**: Average these distances for each individual (yielding DRI<sub>Ind</sub> scaled to [–1, +1]) and then across the group to obtain a single **group DRI**.  

## 4. Overview of Simulations  

To test whether high DRI is unique to high intersubjective consistency, we apply the same calculation to three contrived data scenarios:  

- **A. Purely random data**.  
- **B. Shuffled real data** (considerations and preferences randomly shuffled across participants).  
- **C. Mismatched topics** (considerations on Issue X paired with preferences on Issue Y).  

Under all three scenarios DRI can reach levels exceeding 0.5 even absent any logical link between considerations and preferences. These results are fully [reproducible](https://github.com/social-protocols/dri-in-polis).

### 5 Why Random Data Can Yield High DRI

#### 5.1 Pure noise still hugs the diagonal.  

Imagine 100 people each randomly ranking 10 reasons and 10 preferences -- no logic whatsoever.  

With totally random ranking data, spearman ρs will be normally distributed around 0. This means that when you plot Spearman ρs for considerations (x-axis) against preferences (y-axis), the cloud of points clusters around (0,0), which lies on the 45° diagonal. Hence the distance to the diagonal is small on average and the DRI is positive, as illustrated below.
![Figure 1. IC plot, 10×10 random data (DRI = 0.23)](./random-diffuse.png)


Increasing the number of considerations and preference statements tightens the cluster and thus increases the DRI.

![Figure 1. IC plot, 100×100 random data (DRI = 0.23)](./random-concentrated.png)

<!--
Adding statements that everybody agrees on shifts the cluster upwards and to the right.

![Figure 1. IC plot, mix random and complete agreement](./random-corresponding-agreement.png)
-->

**But you might object:** "Real data is not random!"

---

#### 5.2: Keep the shape, break the logic.  

Shuffle everybody’s preferences (each individual gets the preference of a different randomly selected individual) while leaving their consideration rankings intact.

The marginal distributions (the “shape”) stay the same, but any genuine alignment between reasons and preferences is destroyed. 

Yet in the case below (the GBR future study) the DRI only falls from **0.65** to about **0.63**.

![Figure 2. Actual vs. shuffled GBR Futures (DRI 0.49 → 0.34)](./shuffled-against-standard-18.0.png)

**Objection again:** “Okay, but the DRI is high because there is consistency in overall agreement levels.”

---

### 5.3: Franken-Study: Cross-Topic Pairings

Pair participants' rankings of considerations from one survey (the GBR Futures study) with a different set of users' rankings of preferences from another wholly unrelated survey (the FNQCJ Study).

Pre-deliberation DRI of 0.09 jumps to 0.55 post deliberation.

![Figure 3. Franken-study (GBR vs. FNQCJ) DRI ↦ 0.55](./frankenstudy-18.0-3.0.png)


## 6. Summary of Experiments

| Simulation        | N users | Considerations × Preferences | DRI (example) | Note                   |
|-------------------|---------|--------------------|---------------|------------------------|
| Random (10×10)    | 100     | 10 × 10            | 0.23          | Pure noise             |
| Random (100×100)  | 100     | 100 × 100          | 0.77          | Tighter clustering     |
| Shuffled GBR Futures   | 7     | 35 × 7             | 0.63          | Permuted preferences   |
| Franken GBR→FNQCJ | 7      | 35 × 5             | 0.55          | Cross-topic pairing    |

In all of these experiments, preferences are drawn **randomly** from some distribution -- meaning there is no possible logic linking rankings of considerations to preferences.

The high DRI occurs in all these cases because the points on the DRI plot fall in a:

1. **cluster** that is 
2. **near the diagonal**.

The cluster will be near the diagonal whenever the horizontal mean is approximately the same as the vertical mean.

This can happen spuriously with:

1. **totally random data**, because random Spearman correlations concentrate around (0,0), and  
2. corresponding levels of **overall agreement** on considerations and preferences (e.g. shuffled or cross-topic data), even when no logic ties them together.  

Although genuine intersubjective consistency can also contribute to a high DRI, DRI cannot on its own distinguish between intersubjective consistency and corresponding overall agreement.


### 7. Genuine Intersubjective Consistency

However, not all high DRI values are spurious. For example in the Fremantle Bridge case the high DRI can't be explained away so easily. Unlike in the shuffled GBR Futures experiment or the Frankenstudy, the DRI **falls** significantly to (0.34 ↦ 0.06) when the logical link between consideration and preferences is broken via shuffling. A high Pearson's ρ of .69 adds further evidence of intersubjective consistency. 

![Figure 2. Actual vs. shuffled Fremantle Bridge](./shuffled-vs-standard-pre-post-12.0.png)


## 8. Conclusions & Further Research

1. **Permutation tests for statistical significance.**  
    - For each of the 19 cases in Niemeyer et al. (2024), shuffle profiles 1,000–10,000 times to build a null distribution of DRI and the **DRI delta** (post – pre).  
    - Compute p-values: *Is the observed increase in DRI larger than expected if there was no intersubjective consistency*
1. **Domain-restriction simulations.**  
    - Simulate groups with a known linear logic linking considerations to preferences. Vary the range of generated ρ<sub>C</sub> to test how DRI and Pearson *r* respond.
1. **Cross-case meta-analysis.**  
    - Apply these tests across all 19 forums; examine whether design features (group-building, complexity, decision-impact) predict *statistically significant* gains in intersubjective consistency.
1. **Refining the index.**  
    - Explore a “corrected DRI” subtracting expected diagonal proximity under a null.
    - Compute the Pearson correlation or regression coefficient between ρ<sub>C</sub> and ρ<sub>P</sub> across pairs; compare its sensitivity and robustness to DRI.

These tests may help discern the extent to which a high DRI truly reflects **deliberative reason** rather than a statistical mirage, and may point the way to a refinement of DRI that better measures intersubjective consistency.

All code for reproducing these results is in [this repository](https://github.com/social-protocols/dri-in-polis).

### Invitation to Collaborators

We welcome co-authors interested in advancing a more robust measure of deliberative integration. Please contact me:

- https://x.com/johnwarden
- https://bsky.app/profile/jonathanwarden.bsky.social


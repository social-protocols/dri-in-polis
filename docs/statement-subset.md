# Statement Subset Selection

## The FNQCJ Case: Agreement on Considerations, Disagreement on Preferences

The plots for the FNQCJ case in the How Deliberation Paper stand out: pre-deliberation, the intersubjective agreement values are in a tight vertical band in the right-half of the plot, indicating a high level of intersubjective **agreement on considerations, but disagreement on preferences**.

Below are the plots from the paper for this case (reproduced using the Julia [reference implementation](reference-implementation.md)).

![FNQCJ Case All Considerations](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/published-output/statement-subset/case-3.0-all-considerations-considerations.png)


## Focus on Controversial Considerations

This naturally raises the question of whether there is some subset of consideration that is more controversial, and what the DRI analysis would look like if the survey only included the more controversial considerations? 

More generally, how sensitive are the results to choices made by the survey designers of what questions to include?

This experiment recalculates the DRI using only subset of considerations statements -- simulating what might have happened if the less controversial statements had been left off the original survey. Here's the results for the FNQCJ case using the 25% of considerations with the highest variance in participant ranking.

![FNQCJ Case Top Quartile](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/published-output/statement-subset/case-case-3.0-top-quartile-considerations.png)

Both pre- and post-deliberation, intersubjective-agreement values are now much more spread out horizontally. Otherwise *the charts tell much the same story* as the original chart: in both cases, the main change pre- and post- deliberation is that *the points have moved from the bottom to the top half*. That is, in both cases, there has been an increase in agreement on preferences.

But in the original chart, there is a very large change in DRI, but in the high-controversy statement chart, there is almost no change. 

**So the results are very sensitive to the choice of considerations to include in the survey. Specifically, the inclusion of a large number of non-controversial considerations.**


### Diluting Surveys

Looking from another direction, consider a case where there is a great deal of intersubjective disagreement on considerations. What would happen if the survey were *diluted* with a large number of considerations with broad agreement? Adding enough of these statements would shift values to a narrow band on the right, just like we see in the original FNQCJ plot.

Here's what happens if we only include the **least** controversial statements for the FNQCJ case. The plot becomes even more concentrated in a narrow band on the right-hand side of the plot.

![FNQCJ Case Bottom Quartile](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/published-output/statement-subset/case-case-3.0-bottom-quartile-considerations.png)

Further, the *delta* DRI between pre- and post- is greatly exaggerated

### All Results

The following chart shows the DRI delta (between pre- and post- deliberation), depending on whether the results include the top quartile (most controversial) statements, all statements, or bottom quartile (least controversial) statements.

![DRI Delta Comparison All Cases](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/published-output/statement-subset/delta-dri-comparison.png)


Charts for all cases can be found [here](https://raw.githubusercontent.com/social-protocols/dri-in-polis/master/published-output/statement-subset)





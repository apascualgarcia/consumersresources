### Files

* Networks, where they are: G, and for A both optimized and random **-> "matrices" folder**
* Create a single file with:
   * Network name (G and A used)
   * Topology values (connectances and nestedness)
   * Feasibility volume.
   * Dynamical stability volume.
   * Rate of return.
   **-> This is being computed, you can find an example of it in the results folder ("low-res" simulations)**


### Figures

(Some comments apply for the same figure in both feasibility and dynamical stability, although for some the changes are needed only in one of them)

#### Figure 1 and SM

* We need a couple of examples (two networks) of the surface for feasibility and dynamical stability vs. alpha (as in the thesis). Did the new procedure solve the patchiness of dynamical stability? Particularly relevant the one in which there is a fit for the curve separating the feasible vs unfeasible region (function of S and gamma) for a single alpha. I would select these two networks, eta_G = 0.45 and connectance the highest and the smallest.

* We need to show for a couple of the networks selected the decay of the volume vs. alpha (with the fit) for both feasibility and dynamical stability (two figures, the same two networks I mentioned for example).

* Figures MonteCarlo. We need a figure showing the value of the optimization function vs the step. One in which E becomes negative.

#### Figure 2A and 2B and 1 SM

* Figure `feasible_volume_Nr25_Nc25.pdf` (same for dynamical stability if applies).
    * Remove 10^3  in the axis label and add 10^-3 in the ticks. **-> For some reason (I think because the tick labels are custom) I cannot manage to make it work. If you find a way to make it work, please let me know.**
    * Change the order in the boxes so that optimized and random appear together for each value of alpha. It is misleading now because one tends to compare those boxes that are next to each other but, for optimized and random, contiguous boxes belong to different alpha values **-> I am not sure I see your point, the background alternates between white and grey to show a different value of alpha**. I think the most natural order should be FC-RM-OM because we show in this way the effect of our choices (FC, no choices), RM (a given connectance) and then OM (completely manipulated) **-> That makes sense, I switched RM and OM and will add FC at the start for the final figure**
    * I think the line connecting the boxes is distracting. It connects the means that are more similar than the medians that are clearly distinct. So either connect the medians or simply remove it. **-> I removed them**
   * P-values figures.
       * I think the test we need is one in which we test if both samples come or not from the same distribution. The test you chose can't reject the hypothesis of both distributions being equal (at least as you enunciated it in the title of the figure, it just rejects that one is larger but they still could be equal). We would like to reject the hypothesis that both distros are equal (and, since one has a median larger than the other one, we can conclude that it is signficantly higher). I think we could use a [Mann-Whitney test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Assumptions_and_formal_statement_of_hypotheses).
       * The plot could be one in which the p-value is shown for the comparison OM-RM only (since are those we have doubts on whether they are equal or not), and we indicate whether OM > RM (or RM > OM) with a color or a symbol (e.g. red square means median OM > RM and green circle means median RM > OM). I think we don't need to compare OM and RM with FC because it is clearly different, just OM and RM.   
       * Substitute description of the test from the title by "feasible volume decay". Remove "feasibility" from the y axis (same for dynamical stability)
   * Silly detail but some reviewers are picky, the subscript G and A in eta and kappa should have mathrm format, since it is not an index nor a variable (i.e. not italic). **-> Updated in my code but I haven't generated the relevant figures again yet**

* Figure dominant eigenvalue.
    * y axis - change to "rate of return". Actually, if I remember correctly, you took for each dynamically stable point in the spaceof parameters the eigenvalue, and then you averaged, it would be the "mean rate of return". **-> Yes that's correct, it's done**
    * Are the differences significant? We need the test for the dominant eigenvalue as well. (just OM vs RM)
    * Previous comments reg. the order of the boxes and the x axis/ticks also apply for this figure. **->Same comments**
    * We need to plot dominant eigenvalue for connectance vs. nestedness (same six plots than for dynamical volume).

#### Figure 3 (or 4)

* Effective competition.
    * Generate  again the figures with the correct labels and understand well what are the maths behind. (e.g. dominant eigenvalue of C or A or what). **->Can you please show me the final effective competition definition we agreed on?**

* Functional groups.
    * I will finish it.


.

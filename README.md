# R codes for the paper: *"Neurotrophin-3 produced by motor neurons non-cell autonomously regulate the development of pre-motor interneurons in the developing spinal cord"*

## Neuron quantifications (Figure 1)

**Script**: `S02-mixed_models.R` 

**Analysis summary**: Compare the cell count between two groups: controls versus constitutive or conditional mutants either with a mixed model and involving a random embryo effect if the random effect was
significant or a linear regression if it was deemed non-significant. A correction for multiple
testing with the False Discovery Rate procedure was applied to the obtained p-values.

## Spatial distribution of the motor neurons and ventral IN populations (Figures 2, 4, S1 and S3)

**Scripts**: `R01-2D_densityPlot.Rmd` (2D densityPlot) ; `S01-statistics_2DDensity.R` (2D densityPlot statistics)

**Analysis summary**: DV and ML positions were respectively defined as (dIN ∗ sin αIN)/H and (dIN ∗ cos αIN)/W, where the distance (dIN) and angle (αIN) were measured from the ventral limit of the central canal to the cell soma using the ruler analysis tool in Image J.
ML versus DV values were plotted using R with a 2D kernel density with a contour line made of splines from Bezier curves. For each section and level, a non parametric multivariate ANOVA-type statistic was performed to assess if the locations of DV and ML in the two conditions are different.


## Software versions:

R 4.4.2
lme4 1.1-37
lmerTest 3.1-3
npmv 2.4.1
ggplot2 3.5.2
bezier 1.1.2
pracma 2.4.4


# Ideas and TODO

- Synthetic data - continous mixture of two different populations

- Language - avoid "sample" whenever possible - "draw" (from bootstrap), "observation" (a row in the data)

- Parallel coordinates plots to show more dimensions

- Find a threshold (number of samples kept) where the orginal MDS structure ceases to be recognizable
  - How to evaluate this??? Maybe residual error after procrustes? But what about qualitative structure?

- Use Enterotype example data from phyloseq - they are know to have had issues
  - Generally requires getting the handle on Unifrac





Manuscript outline:

Name: metagenboot: Visual Diagnostics for Metagenomic Datasets via Dirichlet-Multinomial Bootstrapping

# Introduction

- Few methods to visualize confidence in MDS
- A general method - any analysis can be put through resampling and we provide tools to make this easy
- Resampling already suggested in "Ten quick tips for effective dimensionality reduction" https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006907
- Related to DISTATIS (https://personal.utdallas.edu/~herve/abdi-distatis2005-old.pdf) but DISTATIS makes little sense for NMDS, also we are not trying to find a consensus.


jackknife beta diversity (in QIIME jackknifed_beta_diversity.py) - works by simply stacking the PCoA plots over one another (maybe works because the jackknife = removing one sample at a time so the changes are unlikely to be large) + confidence ellipses

- QIIME2 has rarefaction resampling for alpha diversity: https://docs.qiime2.org/2019.1/tutorials/moving-pictures/



# Resampling models

- Jackknife
- Dirichlet-multinomial 
  - Models technical noise
  - Prior choice (max. likelihood / Bayes / fixed)
- Rarefying 
  - Including mostly for completeness and comparisons
  
# Evaluating reliability of dimensionality reduction techniques

Dunthorn 17: We see that the two clusters for Panama are probably not real

Confidence ellipses (as used in QIIME) likely misleading as the model-based resamples (or even rarefied samples) are nowhere resembling an ellipse.

Why not do a big MDS over all samples? This has two problems:
One practical, because NMDS (at least the vegan implementation) tends
to have convergence issues for large datasets and takes looooong to
compute.
Second is theoretical: I've just artificially added similar points to
the dataset. The more samples I take, the more the optimization
function of NMDS will be rewarded for keeping the bootstrap samples
close to the original, even when this means representing the distances
between the original points poorly. And this is actually what happens,
at least to an extent - for example when using the DESeq2 model to
resample Tijana's data, NMDS ran on the big matrix resulted in the
original observations neatly grouped by the DESeq2 predictors, making
the plot hugely different from what you see when you run NMDS on the
originals alone.

Falsifying approach: Reliable clustering survives high-noise resampling. Homogenity must be visible after low-noise resampling.



# Diagnosing Coverage of Biological Variability



# Appendix

- Can be used with other methods: Using resampling for species richness/Shannon diversity


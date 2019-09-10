# Ideas and TODO



- Find a threshold (number of samples kept) where the orginal MDS structure ceases to be recognizable
  - How to evaluate this??? Maybe residual error after procrustes? But what about qualitative structure?

- Bootstrap samples instead of OTUs - run MDS on a subset of samples, use procrustes to match the overlap.

- Try to find a way to plot the are occupied by points and samples instead of the actual points (probably tricky)

- Use Enterotype example data from phyloseq - they are know to have had issues

- Add DISTATIS over resampled as an option

Manuscript outline:

Name: metagenboot: Diagnostics for Metagenomic Datasets via Model-Based Bootstrapping

# Introduction

- Few methods to visualize confidence in MDS
- No (check) way to choose differential abundance library
- A general method - any analysis can be put through resampling and we provide tools to make this easy
- Resampling already suggested in "Ten quick tips for effective dimensionality reduction" https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006907
- Related to DISTATIS (https://personal.utdallas.edu/~herve/abdi-distatis2005-old.pdf) but DISTATIS makes little sense for NMDS, also we are not trying to find a consensus.

# Resampling models

- Dirichlet-multinomial 
  - Models technical noise
  - Prior choice (max. likelihood / Bayes / fixed)
- Rarefying 
  - Including mostly for completeness and comparisons
- DESeq2
  - Requires groups
- DESeq2 - ungrouped
  - TODO
- ANCOM
  - Probably only posssible to sample on the unobserved scale
  
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

Testing with an experiment?
- Run MDS on simulated smaller experiments and the full experiment
- Sho



# Diagnosing Coverage of Biological Variability


# Choosing Differential Abundance Model



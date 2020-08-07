GLOBAL OVERVIEW

HERMES version: 2.2

The HERMES index is a method of quantifying molecular evolution of mitochondrial genomes
in different species and clusters; this method was originally proposed in Plazzi et al.
(2016). The index relies on maximum likelihood factor analysis to summarize different measures
that are typically found to be linked with evolutionary rates; it is intended to be computed
a posteriori, i.e. after a phylogenetic and genomic analysis. As different empirical measures
are merged together in a single score, it is a “hyper-empirical” index; moreover, it is a
relative measure, because all species are compared with an outgroup: therefore, it was called
Hyper-Empirical Relative Mitochondrial Evolutionary Speed (HERMES) index.

The present Python script performs data collection and then calls a dedicated R script to
complete the factor analysis. The mitogenomic features that are currently implemented in
HERMES-v2.2.py are:

the percentage of Unassigned Regions (URs);

the Amount of Mitochondrial Identical Gene Arrangements (AMIGA) index;

the absolute value of the Strand Usage (SU) skew;

the root-to-tip distance computed over a given phylogenetic tree;

the Maximum Likelihood (ML) pairwise distance from a given outgroup;

the AT content;

the AT skew;

the GC skew;

the number of (annotated) genes;

the length of the molecule;

the Codon Adaptation Index (CAI), as defined in Sharp and Li (1987) and Xia (2007).

For each possible combination of at least two of these variables, a factor analysis is carried
out. Normalization and varimax rotation are used, factor scores are found using correlation
preserving, and correlations are found using the Pearson method; given the possible presence
of missing values, missing data are set to be imputed using the median.

All the variables are pooled together for each species into the value of a single loading: we
define this score as the HERMES score of a given species.

The best-performing variable set and the goodness-of-fit of the analysis is assessed following
the recommendations of Hu and Bentler (1999): Tucker-Lewis Index (TLI) greater than 0.95; root
mean square of the residuals (SRMR) smaller than 0.08; root mean squared error of approximation
(RMSEA) less than 0.06; moreover, the  Kaiser-Meyer-Olkin index (KMO) is taken into account on
this regard.

REFERENCES

Hu L-T, Bentler PM. 1999. Cutoff criteria for fit indexes in covariance structure analysis:
Conventional criteria versus new alternatives. Structural Equation Modeling: A Multidisciplinary
Journal. 6:1-55.

Plazzi F, Puccio G, Passamonti M. 2016. Comparative Large-Scale Mitogenomics Evidences
Clade-Specific Evolutionary Trends in Mitochondrial DNAs of Bivalvia. Genome Biol Evol.
8:2544-2564.

Sharp PM, Li WH. 1987. The Codon Adaptation Index--a Measure of Directional Synonymous Codon
Usage Bias, and Its Potential Applications. Nucleic Acids Res. 15:1281-1295.

Xia X. 2007. An Improved Implementation of Codon Adaptation Index. Evol Bioinform. 3:53-58.

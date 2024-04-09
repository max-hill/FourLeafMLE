# Maximum likelihood estimation for quartets under the CFN model

One of the core motivations evolutionary biology is to recover the tree of life that describes how species evolve over time. One common approach is to use DNA sequence alignment data from present day species to estimate the tree topology. By tree topology, we mean the branching of a binary tree. More than that, we can estimate how closely species are related by assign branch weights to the edges of a tree. By fixing a phylogenetic model, we are assuming a tree topology and structure of the branch lengths. Maximum likelihood estimation takes the phylogenetic data and recovers the branch lengths for a particular model. In the many phylogenetic problems there is a large number of species to compare. In virology, where DNA alignment is replaced by RNA sequence alignment there may be a vast number of observed species. To mitigate the computational explosion of statistical inference on these large scale models, practicioners use quartet methods." A quartet method does statistical inference for many subsets of four species and then puzzle pieces the tree of life from those outputs. The number of quartets is quartic in the number of species, so it is important to be able to be able to do this estimation correctly as quickly as possible. However, a big difficulty with MLE is that the likelihood function is non-convex and the maximum may lie on the boundary of the model. This package resolves these two obstacles using computer algebra and the theory of maximum likelihood degrees. Our main contribution puts algebra into practice by implementing a custom tailored solver for biologists. 

## Quartet method motivation
Quartet methods are a canonical method and hold a special place in phylogenetics.

## four_leaf_mle
This is the main function.
It takes site pattern data as an input and returns a vector of information for the mle.
This information include the maximum value of the log-likelihood function. 
It also includes a tree configuration that achieves that. A todo is to make sure we aren't throwing out any trees that also achieve that. 
It also returns the branch lengths. 

## Showcase: our method is fast and/or our method can handle LBA
Two ways that our function is useful. 
One idea to show the former is to compare with heuristics using EM. 
p-values for model selection 

## Felsenstein zone plot 
Our implementation leads to conjectures on the geometry of phylogenetic data. There should be some algebraic curves that cut out the boundaries. 

## Test functions
We have several test functions to ensure correctness of the method

## compute_RX_MLE functions for case analysis
This problem breaks down to computing MLE's for ten different subcases. (This is in contrast to the binary latent class model with only a single latent variable (represented by the internal node) which has many many more subcases; the boundary stratification is very very different from this CFN case with two internal nodes)
The first two cases are solved using homotopy continuation. 
The other eight can be solved using formulas.







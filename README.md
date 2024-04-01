# Maximum likelihood estimation for quartets under the CFN model

Quartet methods are a canonical method 

## four_leaf_mle
This is the main function.
It takes site pattern data as an input and returns a vector of information for the mle.
This information include the maximum value of the log-likelihood function. 
It also includes a tree configuration that achieves that. A todo is to make sure we aren't throwing out any trees that also achieve that. 
It also returns the branch lengths. 

## Test functions
We have several test functions to ensure correctness of the method

## compute_RX_MLE functions for case analysis
This problem breaks down to computing MLE's for ten different subcases. (This is in contrast to the binary latent class model with only a single latent variable (represented by the internal node) which has many many more subcases; the boundary stratification is very very different from this CFN case with two internal nodes)
The first two cases are solved using homotopy continuation. 
The other eight can be solved using formulas.

## Felsenstein zone plot 
TBD

## Comparing to heuristics?





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

### Speed
Consider a phylogenetic tree with 16 leaves. There are 1820 quartets in this tree. If it takes on average one minute to estimate the parameters for a single quartet, then it would take over 30 hours to estimate all of the quartets. This gives us a scale of the importance of being able to estimate a quartet quickly. 

In the following experiment we demonstrate the speed of our method. The biggest takeaway is that our method is able to compute the MLE of a quartet in a sixth of a second on average. 

For generic data, 

we ran the software 1000 times and found the average time solve, min and max times and/or quartiles:

For sparse data, 

we ran the software 1000 times and found the average time solve, min and max times and/or quartiles:

Open question: characterize the data discriminants and log Voronoi cells

### Global optimization 
For large trees it is infeasible to consider the MLEs of all the quartests at once. Instead, practicioners initialize all of the parameters of the tree. Next, they run the EM algorithm to compute a local maximum by only updating one quartet at a time. If the EM algorithm converges to a suboptimal local max, then this can lead to slow downs in reconstructing the large tree. Instead of relying a local method to compute the MLE, we use algebraic methods to find the global optimum when estimating a quartet. This implementation can be viewed as a black box solver for *any* phylogenetic reconsutruction method using quartets.  

(move to quartet motivation section)

### Long branch attraction
p-values for model selection 

## Felsenstein zone and the geometry of data 
Our implementation leads to conjectures on the geometry of phylogenetic data. There should be some algebraic curves that cut out the boundaries. 

## Model stratification and ML degrees: The function compute_RX_MLE functions for case analysis
This problem breaks down to computing MLE's for ten different subcases. (This is in contrast to the binary latent class model with only a single latent variable (represented by the internal node) which has many many more subcases; the boundary stratification is very very different from this CFN case with two internal nodes)
The first two cases are solved using homotopy continuation. 
The other eight can be solved using formulas.


### Test functions
We have several test functions to ensure correctness of the method 







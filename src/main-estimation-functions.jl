###_____________________________________________________________________________
##
## Main functions for computing the MLE or a list of local maxima
##______________________________________________________________________________

"""
    test_data_genericity(SITE_PATTERN_DATA)

# Description

Test whether a data vector is generic. Here, 'generic' means that for each pair of leaves, i and j, there
exists at least one site such that the base pair at leaf i differs from the base pair at leaf j. The
functions in this package cannot handle non-generic data without errors.

# Examples

```
julia> test_data_genericity([1,0,0,0,0,0,0,0])
false

julia> test_data_genericity([1,3,5,0,0,9,5,4])
true

"""
function test_data_genericity(SITE_PATTERN_DATA)
    all([any(>(0),round.(SITE_PATTERN_DATA[subset], digits = 10))
         for subset in [[5,6,7,8], [3,4,7,8], [2,4,6,8], [3,4,5,6], [2,4,5,7], [2,3,6,7]]])
end


"""
    fourLeafMLE(SITE_PATTERN_DATA)

# Description

Attempt to return a unique global maximum by computing the local maxima over the interior of the model as
well as all submodels, and returning a list of only those which maximize the likelihood.

# Arguments

- `SITE_PATTERN_DATA' : A site frequency vector of length 8, specifying the number of times each site
                        pattern is observed in the data. The patterns are ordered as [xxxx, xxxy, xxyx,
                        xxyy, xyxx, xyxy, xyyx, xyyy] where x and y signify different nucleotides. For
                        example, the vector `SITE_PATTERN_DATA=[14,2,5,10,7,4,3,6]` means that pattern
                        xxxx is observed 14 times, pattern xxxy is observed 2 times, pattern xxyx is
                        observe 5 times, and so forth.

# Output

A list of vectors corresponding to trees (or submodels) which maximize the likelihood. Each vector
corresponds to a located maximum and takes the form

```

[logL, reduced tree class, list of compatible binary quartet topologies,
 branch lengths (excluding those which equal 0 or 1), branch length names, labels, verbal description]

```

Here, the 'reduced tree class' refers to a scheme developed by the authors to categorize submodels by the
type of computations required to find their local maxima. There are 10 nontrivial classes (i.e., which
having nonzero likelihoods for generic data), denoted R1, R2, ..., R10. For example, R1 = binary quartet,
R2 = star tree, and R3 through R10 refer to additional submodel types whose details can be found in the
comments of the file `aux-functions-for-R3-through-R10.jl'.

If two or more local maxima are found which maximize the likelihood, then this function will return more
than one maximum likelihood estimator. There are two possible reasons for this (1) the MLE is indeed not
unique, or (2) the floating point arithmetic is insufficiently precise to distinguish which maxima has
greater likelihood.


# Example usage:

```

random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,1)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
fourLeafMLE(SITE_PATTERN_DATA)
listMaxima(SITE_PATTERN_DATA)

```

"""
function fourLeafMLE(SITE_PATTERN_DATA)
    "Attempt to return the global maxima."
    test_data_genericity(SITE_PATTERN_DATA) || error("Data does not satisfy genericity conditions")
    tmp = [compute_R1_and_R2_MLE(SITE_PATTERN_DATA);
           compute_R3_MLE(SITE_PATTERN_DATA); compute_R4_MLE(SITE_PATTERN_DATA);
           compute_R5_MLE(SITE_PATTERN_DATA); compute_R6_MLE(SITE_PATTERN_DATA);
           compute_R7_MLE(SITE_PATTERN_DATA); compute_R8_MLE(SITE_PATTERN_DATA);
           compute_R9_MLE(SITE_PATTERN_DATA); compute_R10_MLE(SITE_PATTERN_DATA)]
    # remove any MLEs with at least one Hadamard parameter indistinguishable from zero:
    tol=10e-9
    filter!(x -> all(i -> i >= tol, x[4]), tmp)
    # keep only the MLEs with the biggest logL:
    maximum_logL = maximum(x[1] for x in tmp) 
    tmp = [x for x in tmp if x[1] == maximum_logL] 
    return tmp
end

"""
    listMaxima(SITE_PATTERN_DATA)

# Description

Return a list of local maxima or critical points, sorted by log-likelihood. The only difference with the
function `fourLeafMLE()' is that `listMaxima' does not filter out critical points which do not achieve
the global optimum. See the documentation for `fourLeafMLE()'.

# Example usage:

```

random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,1)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
listMaxima(SITE_PATTERN_DATA)

```

"""
function listMaxima(SITE_PATTERN_DATA)
    
    test_data_genericity(SITE_PATTERN_DATA) || error("Data does not satisfy genericity conditions")
    tmp = [compute_R1_and_R2_MLE(SITE_PATTERN_DATA);
           compute_R3_MLE(SITE_PATTERN_DATA); compute_R4_MLE(SITE_PATTERN_DATA);
           compute_R5_MLE(SITE_PATTERN_DATA); compute_R6_MLE(SITE_PATTERN_DATA);
           compute_R7_MLE(SITE_PATTERN_DATA); compute_R8_MLE(SITE_PATTERN_DATA);
           compute_R9_MLE(SITE_PATTERN_DATA); compute_R10_MLE(SITE_PATTERN_DATA)]
    # remove any MLEs with at least one Hadamard parameter indistinguishable from zero:
    tol=10e-9
    filter!(x -> all(i -> i >= tol, x[4]), tmp)
    # Sort the array of arrays based on the comparison function
    compare_func(a, b) = (a[1] > b[1])
    sorted_tmp = sort(tmp, lt=compare_func)
    return sorted_tmp
end

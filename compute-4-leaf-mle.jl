# import Pkg; Pkg.add("Distributions")
using Distributions,HomotopyContinuation,LinearAlgebra
using Random # for the function randperm
# using ColorSchemes # for colorscheme
# using ColorTypes
using Plots
include("constants.jl")
include("general-auxilliary-functions.jl")
include("aux-functions-for-homotopy-step.jl")
include("aux-functions-for-R1-and-R2.jl")
include("aux-functions-for-R3-through-R10.jl")

# in emacs, run (setenv "JULIA_NUM_THREADS" "12")
# Threads.nthreads()

###_____________________________________________________________________________
##
## Compute list of MLEs
##______________________________________________________________________________

function test_data_genericity(SITE_PATTERN_DATA)
    all([any(>(0),round.(SITE_PATTERN_DATA[subset], digits = 10))
         for subset in [[5,6,7,8], [3,4,7,8], [2,4,6,8], [3,4,5,6], [2,4,5,7], [2,3,6,7]]])
end

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
    maximum_logL = maximum(x[1] for x in tmp)    # shouldn't I be rounding the logL to like 10 decimal places first?
    tmp = [x for x in tmp if x[1] == maximum_logL] 
    return tmp
end

function listMaxima(SITE_PATTERN_DATA)
    "Return a list of local maxima or critical points, sorted by log-likelihood."
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

#= Example usage:
random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,1)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
fourLeafMLE(SITE_PATTERN_DATA)
listMaxima(SITE_PATTERN_DATA)
=#


###__________________________________________________
##
## Make package
##___________________________________________________
import Pkg; Pkg.add("PkgTemplates")
using PkgTemplates
t= Template(;
            user="max-hill",
            license="MIT",
            authors=["Max Hill", "Jose Israel Rodriguez"],
            plugins=[TravisCI(), Codecov(), Coveralls(), AppVeyor()],)


generate("fourLeafMLE",t)

# first create github repo. Then goto your template repo (~/.julia/dev/MyExample) and run
# git remote set-url origin https://github.com/max-hill/MyExample.git
# then push

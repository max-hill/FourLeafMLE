__precompile__()

module FourLeafMLE

export
    fourLeafMLE,
    listMaxima,
    make_classification_plot,
    generate_averaged_classification_plot_data,
    make_averaged_classification_plot,
    generate_classification_plot_data,
    compute_probability_vector,
    computeProbabilityVector,
    compute_hadamard_parameters,
    test_data_genericity,
    Multinomial # from the package Distributions


using Distributions
using HomotopyContinuation
using LinearAlgebra # for matrix multiplication
using Random # for permutation function in tests
using Pipe
using Plots, ProgressBars, LaTeXStrings # for plotting

include("constants.jl")
include("general-auxilliary-functions.jl")
include("aux-functions-for-homotopy-step.jl")
include("aux-functions-for-R1-and-R2.jl")      # these functions use homotopy
include("aux-functions-for-R3-through-R10.jl") # these don't use homotopy
include("main-estimation-functions.jl")
include("plots.jl")
include("../test/test_functions.jl")

end

# Note to self: to run Julia with multiple threads in the emacs REPL, first run
# (setenv "JULIA_NUM_THREADS" "12").


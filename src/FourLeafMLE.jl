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
    


using Distributions
using HomotopyContinuation
using LinearAlgebra
using Random
using Pipe
using Plots, ProgressBars, LaTeXStrings # for plotting

include("constants.jl")
include("general-auxilliary-functions.jl")
include("aux-functions-for-homotopy-step.jl")
include("aux-functions-for-R1-and-R2.jl")      # these functions use homotopy
include("aux-functions-for-R3-through-R10.jl") # these don't use homotopy
include("main-estimation-functions.jl")
include("plots.jl")

end

#=INSTRUCTIONS FOR USE:
To load the module, run the above code, and then run

`import .FourLeafMLE'

You can then use the exported functions. See docstrings for instructions and
examples of how to use them.
=#


# NOTE: For optimal performance, run Julia with multiple threads.
# To do this in emacs REPL, first run (setenv "JULIA_NUM_THREADS" "12")
# To check the number of treads, run Threads.nthreads() in Julia.


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
## Generate random data
##______________________________________________________________________________

# SITE_PATTERN_DATA = rand(Multinomial(N,8))
random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,1)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
#SITE_PATTERN_DATA = [125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 126.0]

###_____________________________________________________________________________
##
## Compute list of MLEs
##______________________________________________________________________________

function test_data_genericity(SITE_PATTERN_DATA)
    all([any(>(0),round.(SITE_PATTERN_DATA[subset], digits = 10))
         for subset in [[5,6,7,8], [3,4,7,8], [2,4,6,8], [3,4,5,6], [2,4,5,7], [2,3,6,7]]])
end

function fourLeafMLE(SITE_PATTERN_DATA)
    "Attempt to return the global maximum."
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

fourLeafMLE(SITE_PATTERN_DATA)
listMaxima(SITE_PATTERN_DATA)

###_____________________________________________________________________________
##
## Felsenstein zone plot (distance version)
##______________________________________________________________________________

#=
The coloring here isn't good, because there may be more than one MLE type which achieves maximal likelihood.
=#

##
### Plot with hadamard axes
x_values, y_values, z_values = Float64[], Float64[], Float64[]
x_increments,y_increments= collect(0:0.05:0.95), collect(0:0.05:0.95)
@time begin
    for x in x_increments
        for y in y_increments
            push!(x_values,x)
            push!(y_values,y)
            print("$x,$y\n")
            model=computeProbabilityVector([y,x,y,x,x],1)
            z=0
            for i in 1:5
                SITE_PATTERN_DATA=rand(Multinomial(500,model))
                if 1 in fourLeafMLE(SITE_PATTERN_DATA)[1][3] # test if the MLE is compatible with the true topology τ=1
                    z=z+1.0
                end
            end
            push!(z_values,z)
        end
    end
end

z_values = z_values / maximum(z_values) # normalize z
plot()
scatter(x_values, y_values, zcolor=z_values, c = cgrad([:white, :black]), legend=false, markershape=:square, markersize=8, markerstrokewidth=0, aspect_ratio=:equal, xlims=(-.05, 1.0), ylims=(-.05, 1.0))
xlabel!("x ")
ylabel!("y")
title!("Felsenstein zone plot (Hadamard)")


# # plot with probability axes like in Swofford 2001
# p_a,p_b = (1 .- x_values)/2, (1 .- y_values)/2
# plot()
# scatter(p_a, p_b, zcolor=z_values, c = cgrad([:white, :black]), legend=false, markershape=:square, markersize=7, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 0.51), ylims=(0, 0.51))
# xlabel!("p_a")
# ylabel!("p_b", rotate=true)
# title!("Felsenstein zone plot (probability of difference)")



# Include some tests with data from Chor 2003


###_____________________________________________________________________________
##
## Felsenstein zone plot (axes in evolutionary distance)
##______________________________________________________________________________

@time begin
    m=1 # number of samples per parameter regime (replicates)
    x_values, y_values, τ1_scores, τ3_scores, τ2_scores, avg_number_of_maximizers = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    strict_τ_records = Any[]
    step_size=0.01
    dx_increments, dy_increments = collect(0.01:step_size:1.49), collect(0.01:step_size:1.49) # distance increments
    for x in exp.(-2*dx_increments)
        for y in exp.(-2*dy_increments)
            push!(x_values,x)
            push!(y_values,y)
            model=computeProbabilityVector([y,x,y,x,x],1)
            model=round.(model, digits = 16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            print("x=",x,"; y=",y,"\n",model,"\n")
            τ1_score, τ3_score, τ2_score = zeros(3)
            number_of_maximizers=0.0
            strict_τ_record = zero(1)
            for i in 1:m
                SITE_PATTERN_DATA=rand(Multinomial(500,model))
                estimates=fourLeafMLE(SITE_PATTERN_DATA)
                number_of_maximizers = number_of_maximizers + length(estimates)
                congruent_topologies = estimates[1][3] # Note we only look at the first MLE - there might be more than one!

                # this next block works only if m=1. If a unique binary topoogy is returned, it records the topology in strict_τ_records.
                if length(estimates) != 1
                    push!(strict_τ_records,"multiple MLEs") # more than one MLE is returned
                else
                    if length(congruent_topologies)==1 # the MLE is unique and displays a quartet split 12|34, 13|24, or 14|23.
                        if (1 in congruent_topologies) push!(strict_τ_records,"unique τ1") end
                        if (2 in congruent_topologies) push!(strict_τ_records,"unique τ2") end
                        if (3 in congruent_topologies) push!(strict_τ_records,"unique τ3") end
                    else
                        push!(strict_τ_records,"unique bdry case") # the unique MLE is a boundary case
                    end
                end

                if 1 in congruent_topologies # test if the MLE is compatible with the true topology τ=1. 
                     τ1_score= τ1_score+1.0
                end
                if 2 in congruent_topologies # test if the MLE is compatible with the LBA  topology τ=2. 
                     τ2_score= τ2_score+1.0
                end
                if 3 in congruent_topologies # test if the MLE is compatible with the LBA  topology τ=3. 
                    τ3_score= τ3_score+1.0
                end
            end
            push!(τ1_scores,τ1_score)
            push!(τ2_scores,τ2_score)
            push!(τ3_scores,τ3_score)
            push!(avg_number_of_maximizers,number_of_maximizers/m)
        end
    end
end # 583 seconds to compute 3610 MLE estimates = about 6 estimates per second per core. Multithreading
    # with 12 cores, we get about 54 per second. (This was based on a test of 44404 trees that took 8182 seconds.)

length(x_values)
# convert branch lengths to evolutionary distances and normalize z

x_dist=-(1/2)*log.(x_values)
y_dist=-(1/2)*log.(y_values)
τ3_scores = τ3_scores / maximum(τ3_scores)
τ2_scores = τ2_scores / maximum(τ2_scores)
τ1_scores = τ1_scores / maximum(τ1_scores)

# # save data for later
# using FileIO, JLD2
# FileIO.save("x_dist.jld2","x_dist",x_dist)
# FileIO.save("y_dist.jld2","y_dist",y_dist)
# FileIO.save("τ3_scores.jld2","τ3_scores",τ3_scores)
# FileIO.save("τ2_scores.jld2","τ2_scores",τ2_scores)
# FileIO.save("τ1_scores.jld2","τ1_scores",τ1_scores)
# FileIO.save("avg_number_of_maximizers.jld2","avg_number_of_maximizers",avg_number_of_maximizers)

# x_dist = FileIO.load("x_dist.jld2","x_dist")
# y_dist = FileIO.load("y_dist.jld2","y_dist")
# τ3_scores = FileIO.load("τ3_scores.jld2","τ3_scores")
# τ2_scores = FileIO.load("τ2_scores.jld2","τ2_scores")
# τ1_scores = FileIO.load("τ1_scores.jld2","τ1_scores")
# avg_number_of_maximizers = FileIO.save("avg_number_of_maximizers.jld2","avg_number_of_maximizers")


plot()
p1 = scatter(x_dist, y_dist, zcolor=τ1_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y",
             title="Estimated topology concordance with 12|34", colorbar=true,
             colorbar_title="Proportion having topology 12|34")


p2 = scatter(x_dist, y_dist, zcolor=τ2_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 13|24",
             colorbar=true)

p3 = scatter(x_dist, y_dist, zcolor=τ3_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 14|23",
             colorbar=true)

p4 = scatter(x_dist, y_dist, zcolor=avg_number_of_maximizers, c = cgrad([:white, :black]), legend=false,
            markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
            ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Avg number of global maxima", colorbar=true)



# Initialize the plot
p5 = plot(legend = :outertopright)

# Define your custom color gradient
color_mapping = Dict("multiple MLEs" => :white, "unique bdry case" => :black, "unique τ1" => :green, "unique τ2" => :red, "unique τ3" => :blue)

# Add each category to the plot separately
for category in ["unique τ1", "unique τ2", "unique τ3", "multiple MLEs", "unique bdry case"]
    mask = strict_τ_records .== category
    scatter!(p5, x_dist[mask], y_dist[mask],
             c = color_mapping[category],
             label = category,
             markershape=:square, markersize=1.2, markerstrokewidth=0)
end

# Set the remaining plot attributes
plot!(p5, aspect_ratio=:equal, xlims=(0, 1.501), ylims=(0, 1.501),
      xlabel="parameter x", ylabel="parameter y", title="Avg number of global maxima")












##
### Plot of how frequently boundary cases are inferred
sequence_length = 200
replicates_per_parameter_regime = 10
@time begin
    x=[]
    y=[]
    z=Float64[]
    d = 0.05*(1:9) # desired distance intervals. Note that taking length zero distances results in an error which I don't understand.
    for a in exp.(-2*d)
        for b in exp.(-2*d)
            push!(x,a)
            push!(y,b)
            model=computeProbabilityVector([b,a,b,a,a],1)
            model=round.(model, digits=16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            print("a=",a,"; b=", b, "\n",model,"\n")
            color_val=0
            for i in 1:replicates_per_parameter_regime
                SITE_PATTERN_DATA=rand(Multinomial(sequence_length,model))
                if fourLeafMLE(SITE_PATTERN_DATA)[1][2] != "R1"
                    color_val=color_val+1.0
                end
            end
            push!(z,color_val)
        end
    end
end


# convert branch lengths to evolutionary distances and normalize z
x_dist=-(1/2)*log.(x)
y_dist=-(1/2)*log.(y)
z = z / maximum(z)

# plot
color_gradient = cgrad([:red, :blue])
plot()
scatter(x_dist, y_dist, zcolor=z, c = color_gradient, legend=false, markershape=:square, markersize=10, markerstrokewidth=.2)
xlabel!("a")
ylabel!("b")
title!("Boundary case frequency plot (evolutionary dist.)")
savefig("boundary-case-frequency-",sequence_length,"bp-",replicates_per_parameter_regime," replicates")




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

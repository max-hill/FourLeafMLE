# import Pkg; Pkg.add("Distributions")
using Distributions,HomotopyContinuation,LinearAlgebra
using Random # for the function randperm
# using ColorSchemes # for colorscheme
# using ColorTypes
using Plots
include("constants.jl")
include("general-auxilliary-functions.jl")
include("aux-functions-for-homotopy-step.jl")
include("aux-functions-for-R3-through-R10.jl")

###__________________________________________________
##
## Generate random data
##___________________________________________________

# SITE_PATTERN_DATA = rand(Multinomial(N,8))
random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,1)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
#SITE_PATTERN_DATA = [125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 125.0, 126.0]

###__________________________________________________
##
## Compute list of MLEs
##___________________________________________________
SITE_PATTERN_DATA=test_model
function fourLeafMLE(SITE_PATTERN_DATA)
    tmp = [compute_R1_MLE(SITE_PATTERN_DATA); compute_R2_MLE(SITE_PATTERN_DATA);
           compute_R3_MLE(SITE_PATTERN_DATA); compute_R4_MLE(SITE_PATTERN_DATA);
           compute_R5_MLE(SITE_PATTERN_DATA); compute_R6_MLE(SITE_PATTERN_DATA);
           compute_R7_MLE(SITE_PATTERN_DATA); compute_R8_MLE(SITE_PATTERN_DATA);
           compute_R9_MLE(SITE_PATTERN_DATA); compute_R10_MLE(SITE_PATTERN_DATA)]
    # keep only the MLEs with the biggest logL:
    maximum_logL = maximum(x[1] for x in tmp)    # shouldn't I be rounding the logL to like 10 decimal places first?
    tmp = [x for x in tmp if x[1] == maximum_logL] 
    # remove any MLEs with at least one Hadamard parameter indistinguishable from zero:
    filter!(x -> all(i -> i >= 100000*eps(), x[2]), tmp)
    return tmp
end


SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
@time fourLeafMLE(SITE_PATTERN_DATA)

###__________________________________________________
##
## True model tests
##___________________________________________________

# FIXME: the following does not return correct result
fourLeafMLE([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])


# Test whether the edge parameters can be correctly estimated from the true model when the true model is a quartet tree with random topology and uniformly-chosen random edge parameters

function test_random_binary_quartet_model()
    random_θ=rand(5)
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    length(output) != 1 && error("MLE is not unique",output) # check if there is a unique MLE
    estimated_θ, estimated_quartet_topology = output[1][2], output[1][3]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_θ,"\n",output) # check if the estimate is consistent with the true topology
    euclidean_distance(random_θ,estimated_θ) > 1000000*eps() && error("Incorrect edge parameter estimate:\n", random_θ,"\n", estimated_θ,"\n", euclidean_distance(random_θ,estimated_θ))
    return true
end
y = [test_random_binary_quartet_model() for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test

# Test whether the estimator returns a unique MLE whose estimated topology is consistent with the true tree, when the model has one ore more infinitely-long branches (and the model is used as the data)
function test_random_infinite_branch_model(number_of_zero_parameters)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be zero
    for i in randperm(5)[1:number_of_zero_parameters]
        random_θ[i] = 0
    end
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return output, random_θ
    end
    estimated_θ, estimated_quartet_topology = output[1][2], output[1][3]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_θ,"\n",output) # check if the estimate is consistent with the true topology
    return true
end
# when one parameter is zero
y = [test_random_infinite_branch_model(1) for i in 1:100];[x for x in y if x != 1] # if this list is empty, it passed the test

# when two parameters are zero
y = [test_random_infinite_branch_model(2) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test

# when three parameters are zero
y = [test_random_infinite_branch_model(3) for i in 1:10]; [x for x in y if x != 1] # if this list is empty, it passed the test




# Test whether the estimator returns a unique MLE whose estimated topology is consistent with the true tree, when the model has one ore more length zero branches (and the model is used as the data)
function test_random_zero_branch_model(number_of_one_parameters)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be zero
    for i in randperm(5)[1:number_of_one_parameters]
        random_θ[i] = 1
    end
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return output, random_θ
    end
    estimated_θ, estimated_quartet_topology = output[1][2], output[1][3]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_θ,"\n",output) # check if the estimate is consistent with the true topology
    return true
end
y = [test_random_zero_branch_model(1) for i in 1:1000]; [x for x in y if x != 1] # if this list is empty, it passed the test
y = [test_random_zero_branch_model(2) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test



# problem case -- possible error with how R3 is computed?
θ=[1.0, 0.07760820166220417, 0.9450935684087735, 1.0, 0.5329313072157882]
test_model=computeProbabilityVector(θ,3)
fourLeafMLE(test_model)

logL = test_model'log.(computeProbabilityVector(θ,3))
computeProbabilityVector([0.5036699508533216, 0.9450935684087723, 1, 0.5329313072157887, 0.07760820166220414],1)


function test_random_boundary_model(k)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be zero
    for i in randperm(5)[1:number_of_parameters_on_boundary]
        random_θ[i] = rand(0:1)
    end
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return output, random_θ
    end
    estimated_θ, estimated_quartet_topology = output[1][2], output[1][3]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_quartet_topology, random_θ,"\n",output) # check if the estimate is consistent with the true topology
    return true
end
y = [test_random_zero_branch_model(2) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test



###__________________________________________________
##
## Felsenstein zone plot
##___________________________________________________
#=
The coloring here isn't good, because there may be more than one MLE type which achieves maximal likelihood.
=#

##
### Plot with hadamard axes
x=[]
y=[]
z=Float64[]

for a in [0.05i for i in 0:19]
    print(" $a")
    for b in [0.05i for i in 0:19]
        push!(x,a)
        push!(y,b)
        model=computeProbabilityVector([b,a,b,a,a],1)
        color_val=0
        for i in 1:10
            SITE_PATTERN_DATA=rand(Multinomial(500,model))
            if 1 in fourLeafMLE(SITE_PATTERN_DATA)[1][3]
                color_val=color_val+1.0
            end
        end
        push!(z,color_val)
    end
end
color_gradient = cgrad([:red, :blue])

z = z / maximum(z) # normalize z

plot()
scatter(x, y, zcolor=z, c = color_gradient, legend=false, markershape=:square, markersize=10, markerstrokewidth=0)
xlabel!("a")
ylabel!("b")
title!("Felsenstein zone plot")
# compare with plot where only binary trees are recorded?

# Include some tests with data from Chor 2003


##
### Plot with evenly-space distance axes
@time begin
    x=[]
    y=[]
    z=Float64[]
    d = 0.05*(1:19) # desired distance intervals. Note that taking length zero distances results in an error which I don't understand.
    for a in exp.(-2*d)
        for b in exp.(-2*d)
            push!(x,a)
            push!(y,b)
            model=computeProbabilityVector([b,a,b,a,a],1)
            model=round.(model, digits = 16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            print("a=",a,"; b=", b, "\n",model,"\n")
            color_val=0
            for i in 1:10
                SITE_PATTERN_DATA=rand(Multinomial(500,model))
                if 1 in fourLeafMLE(SITE_PATTERN_DATA)[1][3]
                    color_val=color_val+1.0
                end
            end
            push!(z,color_val)
        end
    end
end # 583 seconds to compute 3610 MLE estimates = 6.2 estimates per second

# convert branch lengths to evolutionary distances and normalize z
x_dist=-(1/2)*log.(x)
y_dist=-(1/2)*log.(y)
z = z / maximum(z)

# plot
color_gradient = cgrad([:red, :blue])
plot()
scatter(x_dist, y_dist, zcolor=z, c = color_gradient, legend=false, markershape=:square, markersize=10, markerstrokewidth=0)
xlabel!("a")
ylabel!("b")
title!("Felsenstein zone plot (evolutionary distance)")


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
                if fourLeafMLE(SITE_PATTERN_DATA)[1][4] != "R1"
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

###_____________________________________________________________________________
##
## test_functions.jl -- functions used for testing by runtests.jl
##______________________________________________________________________________

"""
   test_random_binary_quartet_model(branch_length_distribution)

# Description

Test whether the edge parameters can be correctly estimated from the true model
when the true model is a quartet tree with random topology and uniformly-chosen
random edge parameters. Success occurs when (1) the estimate is unique, (2) its
topology is compatible with that of the true topology model parameter, (3) a
boundary model is not inferred, and (4) the Euclidean distance between the
estimated edge parameters and true edge parameters is smaller than 1e-10. This
test is probabilistic and occasionally fails due to the limits of approximate
arithmetic on very pathological trees (i.e. trees very close to the boundary of
the model.


# Arguments

- `branch_length_distribution` : the probability distribution from which 
                                 branch lengths are drawn. Defaults to 
                                 uniform.

# Output

If all 4 conditions are met, the test returns 'nothing', otherwise it returns a
vector containing details about what went wrong.

"""
function test_random_binary_quartet_model(branch_length_distribution=:uniform)
    if branch_length_distribution == :uniform
        random_θ=rand(5)
    else
        random_θ=rand(branch_length_distribution,5)
    end
    random_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_topology)
    output=fourLeafMLE(test_model)
    # check if there is a unique MLE:
    if length(output) ≠ 1
        return ["MLE is not unique",random_θ,random_topology,output]
    end
    estimated_quartet_topology, estimated_θ = output[1][3], output[1][4]
    # check if the estimate is consistent with the true topology:
    if !(random_topology in estimated_quartet_topology)
        return ["Incorrect topology estimate",random_θ,random_topology,output]
    end
    tol=1e-10
    # check if a boundary model is inferred
    if length(random_θ) ≠ length(estimated_θ)
        return ["Boundary model inferred. Inferrd topology is compatible with"*
            "true topology",random_θ,random_topology,output]

        if euclidean_distance(random_θ,estimated_θ) > tol
            return ["Inaccurate edge parameter estimate", random_θ,
                    estimated_θ,euclidean_distance(random_θ,estimated_θ)]
        end
        return nothing
    end
end

        


"""
   test_random_infinite_branch_model(N_branches)

# Description

Test whether the estimator returns a unique MLE whose estimated topology is
compatible with the true tree, when the model one or more more infinitely-long
branches (and the model is used as the data)

# Arguments

- `N_branches` : Number of branches of infinite length (the branches are
                 randomly-chosen)

# Output

Returns `nothing` if the test is successful. Otherwise it returns a vector
containing details about what went wrong.

"""
function test_random_infinite_branch_model(N_branches)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be 0
    for i in randperm(5)[1:N_branches] # Random
        random_θ[i] = 0
    end
    random_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return ["MLE not unique", random_θ, random_topology, output]
    end
    estimated_topology, estimated_θ = output[1][3], output[1][4]
    # Check if the estimate is compatible with the true topology
    if (random_topology in estimated_topology)
        return nothing
    else
        return ["Incorrect topology estimate", random_θ, random_topology,
                output]
    end
end


"""
   test_random_zero_branch_models(branch_length_distribution)

# Description

Generates a random tree tree model (with random topology, and random branch
lengths, and with 1-2 randomly-chosen branches having lengths zero
evolutionary distance). Using the model as data, it then tests whether the
estimator returns a unique MLE whose estimated topology is consistent with
the true tree.

# Arguments (optional)

- `branch_length_distribution` : the probability distribution from which 
                                 branch lengths are drawn. Defaults to 
                                 uniform.

# Output

If the MLE is not unique, or if the MLE incorrectly estimated the topology,
a vector is returned which contains details about what happened. Otherwise,
the function returns nothing.

"""
function test_random_zero_branch_model(branch_length_distribution=:uniform)
    # decide how many branches should have length zero:
    number_of_zero_length_branches = rand([1,2])
    # initialize variables outside scope of while loop:
    random_quartet_topology, random_θ, test_model = [], [], []
    # generate a generic model with random branch lengths and 1 or 2
    # randomly-chosen length-zero branches:
    while true
        # generate a random model, using either uniformly-chosen branch
        # lengths or the user-supplied distribution:
        if branch_length_distribution == :uniform
            random_θ=rand(5)
        else
            random_θ=rand(branch_length_distribution,5)
        end
        for i in randperm(5)[1:number_of_zero_length_branches]
            random_θ[i] = 1
        end
        random_quartet_topology=rand([1,2,3])
        test_model=computeProbabilityVector(random_θ,random_quartet_topology)
        # check if the resulting model is generic, and break the loop if it is:
        if test_data_genericity(test_model) break end
    end
    # compute the global MLE for the model:
    output=fourLeafMLE(test_model)
    # Case 1. True model is a star tree
    # if the true model doesn't have a binary tree topology (because the
    # internal branch has length 0), then just check if the MLE exists and
    # is unique:
    if (random_θ[5] == 1)
        if (length(output) == 1)
            return nothing
        else
            # Bad case:
            return [random_quartet_topology, random_θ, output,
                    "Star tree topology. Unique MLE not found."]
        end
    end
    # Case 2. True model is not a star tree
    # check if the estimate(s) have topologies congruent with the true topology:
    if all([(random_quartet_topology in x[3]) for x in output])
        # Return 'nothing' if a unique, topologically correct MLE has been found:
        if length(output) == 1
            return
            nothing
        end
        # If the MLE is not unique, explain what happened:
        if length(output) ≠ 1
            return [random_quartet_topology, random_θ, output,
                    "MLE not unique, but all topology estimates are correct "*
                        "(i.e., compatible with true topology."]
        end
    else
        # Bad case:
        return [random_quartet_topology, random_θ, output,
                "MLE may or may not be unique, and at least one estimate has "*
                    "the wrong topology."]
    end
end




#= Unfinished test:

function test_random_boundary_model(k)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be zero or one (randomly chosen)
    for i in randperm(5)[1:k]
        random_θ[i] = rand(0:1)
    end
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return output, random_θ
    end
    estimated_quartet_topology, estimated_θ = output[1][3], output[1][4]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_quartet_topology, random_θ,"\n",output) # check if the estimate is consistent with the true topology
    return true
end

test_result = [test_random_boundary_model(2) for i in 1:10]
@test all(test_result)

=#

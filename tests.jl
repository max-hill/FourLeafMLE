
###_____________________________________________________________________________
##
## True model tests
##______________________________________________________________________________

# The following line should return an error:
fourLeafMLE([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

# Test whether the edge parameters can be correctly estimated from the true model when the true model is
# a quartet tree with random topology and uniformly-chosen random edge parameters
function test_random_binary_quartet_model()
    random_θ=rand(5)
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    # check if there is a unique MLE:
    if length(output) ≠ 1
        error("MLE is not unique",output)
    end
    estimated_quartet_topology, estimated_θ = output[1][3], output[1][4]
    # check if the estimate is consistent with the true topology:
    if !(random_quartet_topology in estimated_quartet_topology)
        error("Incorrect topology estimate\n", random_θ,"\n",output) 
    end
    tol=1e-10
    if euclidean_distance(random_θ,estimated_θ) > tol
        error("Incorrect edge parameter estimate:\n", random_θ,"\n", estimated_θ,"\n", euclidean_distance(random_θ,estimated_θ))
    end
    return true
end
y = [test_random_binary_quartet_model() for i in 1:1000]; [x for x in y if x != 1] # if this list is empty, it passed the test

# Test whether the estimator returns a unique MLE whose estimated topology is consistent with the true
# tree, when the model has one ore more infinitely-long branches (and the model is used as the data)
function test_random_infinite_branch_model(number_of_infinitely_long_branches)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be 0
    for i in randperm(5)[1:number_of_infinitely_long_branches]
        random_θ[i] = 0
    end
    random_quartet_topology=rand([1,2,3])
    test_model=computeProbabilityVector(random_θ,random_quartet_topology)
    output=fourLeafMLE(test_model)
    if length(output) != 1
        return output, random_θ
    end
    estimated_quartet_topology, estimated_θ = output[1][3], output[1][4]
    (random_quartet_topology in estimated_quartet_topology) || error("Incorrect topology estimate\n", random_θ,"\n",output) # check if the estimate is consistent with the true topology
    return true
end

# when one parameter is zero
y = [test_random_infinite_branch_model(1) for i in 1:100];[x for x in y if x != 1] # if this list is empty, it passed the test

# when two parameters are zero
y = [test_random_infinite_branch_model(2) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test

# when three parameters are zero
y = [test_random_infinite_branch_model(3) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test

"""
   test_random_zero_branch_models()

# Description
Generates a random tree tree model (with random topology, and random branch lengths, and with 1-2
randomly-chosen branches having lengths zero evolutionary distance). Using the model as data, it then
tests whether the estimator returns a unique MLE whose estimated topology is consistent with the true
tree.

# Output
If the MLE is not unique, or if the MLE incorrectly estimated the topology, a vector is returned
which contains details about what happened. Otherwise, the function returns nothing.

"""
function test_random_zero_branch_model(branch_length_distribution = :uniform)
    # decide how many branches should have length zero:
    number_of_zero_length_branches = rand([1,2])
    # initialize variables outside scope of while loop:
    random_quartet_topology, random_θ, test_model = [], [], []
    # generate a generic model with random branch lengths and 1 or 2 randomly-chosen length-zero branches:
    while true
        # generate a random model, using either uniformly-chosen branch lengths or the user-supplied
        # distribution:
        if branch_length_distribution == :uniform
            random_θ=rand(5)
        else
            random_θ=rand(branch_length_distribution,5)
        end
        for i in randperm(5)[1:number_of_zero_length_branches] random_θ[i] = 1 end
        random_quartet_topology=rand([1,2,3])
        test_model=computeProbabilityVector(random_θ,random_quartet_topology)
        # check if the resulting model is generic, and break the loop if it is:
        if test_data_genericity(test_model) break end
    end
    # compute the global MLE for the model:
    output=fourLeafMLE(test_model)
    # Case 1. True model is a star tree
    # if the true model doesn't have a binary tree topology (because the internal branch has length 0),
    # then just check if the MLE exists and is unique:
    if (random_θ[5] == 1)
        if (length(output) == 1)
            return nothing
        else
            # Bad case:
            return [random_quartet_topology, random_θ, output, "Star tree topology. Unique MLE not found."]
        end
    end
    # Case 2. True model is not a star tree
    # check if the estimate(s) have topologies congruent with the true topology:
    if all([(random_quartet_topology in x[3]) for x in output])
        # Return 'nothing' if a unique, topologically correct MLE has been found:
        if length(output) == 1 return nothing end
        # If the MLE is not unique, explain what happened:
        if length(output) ≠ 1
            return [random_quartet_topology, random_θ, output, "MLE not unique, but all estimates have topology consistent with the true topology."]
        end
    else
        # Bad case:
        return [random_quartet_topology, random_θ, output, "MLE may or may not be unique, and at least one of the estimates has the wrong topology."]
    end
end
    

# Run the test with uniformly-chosen branch lengths
test_output=[x for x in [test_random_zero_branch_model() for i in 1:100] if x ≠ nothing] 
test_result = (test_output == Union{}[]) # test_result == true, then the test passed. 

# Next, run the test with branch lengths from a distribution which is skewed towards the center (to make
# the estimation problem harder)
alpha, beta = 0.6, 0.6; skewed_distribution = Beta(alpha,beta)
histogram(rand(skewed_distribution, 1000), label="α=$alpha, β=$beta, (Skewed towards center)", nbins = 30, legend=:top,
          xlabel="Value", ylabel="Frequency", title="Edge parameter distribution (skewed)")
@time test_output=[x for x in [test_random_zero_branch_model(skewed_distribution) for i in 1:1000] if x ≠ nothing]
# Commentary: Failures here appear to occur only with very extreme choice of branch lengths due to
# rounding errors, or where two or more models are so similar that tiebreaking is not possible due to the
# limitation imposed by approximate calculations. For this skewed beta distribution, the error rate
# appears to be on the order of 1 in 1000, and is affected by adjusting tolerances. Error rate increases
# if we decrease parameters α and β, since this leads to many edge parameters which are extremely close
# to zero or to one.

# y=test_output[1]
# y[1]
# y[2]
# y[3]
# y[4]
# bad_model=computeProbabilityVector(y[2],y[1])
# listMaxima(bad_model)
# fourLeafMLE(bad_model)


function test_random_boundary_model(k)
    random_θ=rand(5)
    # Next, set some randomly-chosen Hadamard parameters to be zero or one (randomly chosen)
    for i in randperm(5)[1:number_of_parameters_on_boundary]
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
y = [test_random_zero_branch_model(2) for i in 1:100]; [x for x in y if x != 1] # if this list is empty, it passed the test

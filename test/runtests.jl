###_____________________________________________________________________________
##
## Tests
##______________________________________________________________________________

using Test, Distributions, Random
using FourLeafMLE
include("test_functions.jl")

@testset "Is non-generic data handled correctly?" begin
    # The intended behavior is for non-generic data to result in a DomainError:
    @test_throws DomainError fourLeafMLE([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
end


@testset "Is estimation accurate for random trees in the interior?" begin
    # The easy test
    test_result1 = [test_random_binary_quartet_model(:uniform) for i in 1:1000]

    # The harder test by randomly choosing branch lengths closer to the
    # boudnaryd. (See commentary below).
    alpha, beta = 0.6, 0.6
    skewed_distribution = Beta(alpha,beta)
    test_result2 = [test_random_binary_quartet_model(skewed_distribution) for i in 1:1000]

    # Criterion for passing the tests: >99% success rate
    @test all(==(nothing),test_result1)
    @test sum([x==nothing for x in test_result2])/length(test_result2) > .99

    #= COMMENTARY: 
    In the second test above, the beta distribution was chosen so as to sample
    trees near the boundary, which we expect to be harder to estimate than trees
    whose branch lengths lie well within the interior of the model. To see a
    picture of this distribution, run the following:

    histogram(rand(skewed_distribution, 1000), nbins = 30, legend=:top,
              label="α=$alpha, β=$beta, (Skewed towards center)", 
              xlabel="Value", ylabel="Frequency",
              title="Edge parameter distribution (skewed)")

    Due to its random nature, estimation in the above cases may occasionally
    fail--and more frequently in the harder test--which accounts for the use of
    the 99% cutoff for "pass". Based on my observations, failures here seem to
    occur only with very extreme choice of branch lengths due to rounding
    errors, or where two or more models are so similar that tiebreaking is not
    possible due to the limitation imposed by approximate calculations.

    If failure occurs, one can inspect the bad case with the following code:

    test_failures=[x for x in test_result if x ≠ nothing]
    y=test_failures
    bad_model=computeProbabilityVector(y[2],y[1])
    listMaxima(bad_model)
    fourLeafMLE(bad_model)
    =#
end


@testset "Is estimation accurate for random trees on the boundary?" begin
    # When j branches are infinitely long, and the other 5-j branches have
    # uniformly-chosen branch lengths.
    test_result=Any[0,0,0]
    for j in 1:3
        test_result[j] = [test_random_infinite_branch_model(j) for i in 1:100]
        @test all(==(nothing),test_result[j])
    end
    
    # Uniformly-chosen branch lengths with 1-2 length zero branches
    test_result4 = [test_random_zero_branch_model() for i in 1:100] 
    @test all(==(nothing),test_result4)

    # A more stringent test: run the test on random trees near the boundary (by
    # drawing branch lengths fom a beta distribution which is skewed toward the
    # center).
    alpha, beta = 0.6, 0.6; skewed_distribution = Beta(alpha,beta)
    test_result5 = [test_random_zero_branch_model(skewed_distribution) for i in 1:1000]
    @test sum([x==nothing for x in test_result5])/length(test_result5) > .99
    # For this skewed beta distribution, the error rate appears to be on the
    # order of 1 in 1000 tests, and is affected by adjusting tolerances. Error
    # rate increases if we decrease parameters α and β, since this leads to many
    # edge parameters which are extremely close to zero or to one.   
end

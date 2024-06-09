###_____________________________________________________________________________
##
## experimental-heuristic-method.jl
##______________________________________________________________________________

#= COMMENTARY: This document contains the beginning of code to test whether empirical estimates of
certian semi-algebraic constraints of the model can predict the qualitative behavior of maximum
likelihood estimation, or at least provide a faster heuristic-based method than the current approach
taken in this package (which is to solve the MLE problem on all submodels).

Preliminary tests indicate that the use of the semi-algebraic constraints is method is shockingly
accurate at predicting the MLE, at least for trees with uniformly-distributed edge parameters.

The code and tests in this document are currently under development.

=#


M⁺(i,j,s) = [σ[i]==σ[j] for σ∈SP_const]'s
M⁻(i,j,s) = [σ[i]!=σ[j] for σ∈SP_const]'s

function compute_R1_and_R2_MLE__hopefully(SITE_PATTERN_DATA)
    R1_output=[]
    tol=10e-10
    
    B_12 = (M⁺(1,2,SITE_PATTERN_DATA) - M⁻(1,2,SITE_PATTERN_DATA))
    B_13 = (M⁺(1,3,SITE_PATTERN_DATA) - M⁻(1,3,SITE_PATTERN_DATA))
    B_14 = (M⁺(1,4,SITE_PATTERN_DATA) - M⁻(1,4,SITE_PATTERN_DATA))
    B_23 = (M⁺(2,3,SITE_PATTERN_DATA) - M⁻(2,3,SITE_PATTERN_DATA))
    B_24 = (M⁺(2,4,SITE_PATTERN_DATA) - M⁻(2,4,SITE_PATTERN_DATA))
    B_34 = (M⁺(3,4,SITE_PATTERN_DATA) - M⁻(3,4,SITE_PATTERN_DATA))

    topologies_to_analyze = []
    a1 = B_12*B_34
    a2 = B_13*B_24
    a3 = B_14*B_23
    if a1 ≥ max(a2,a3)
        push!(topologies_to_analyze,1)
    end
    if a2 ≥ max(a1,a3)
        push!(topologies_to_analyze,2)
    end
    if a3 ≥ max(a1,a2)
        push!(topologies_to_analyze,3)
    end
    if a1==a2==a2
        push!(topologies_to_analyze,4)
    end
    
    for τ in topologies_to_analyze
        # First, use homotopy to compute a list of real critical points to the likelihood equations
        cp_list = compute_R1_or_R2_critical_points_with_homotopy(SITE_PATTERN_DATA,τ)
        # Second, remove those critical points which correspond to negative or zero
        # θ values (i.e. critical points not in the interior of the FBM model).
        isValid_for_τ(c)=isValid(c,τ,tol=1.0e-14)
        filter!(i->isValid_for_τ(i),cp_list)
        # valid_cp_list = filter(i->isValid(i),cp_list)
        # Third, if there are no valid critical points, return an exceptional case for topology τ
        if isempty(cp_list)
            empty_result = ((τ in [1,2,3]) ?
                [-Inf, "R1", [], [], "None of the critical points for R1 with topology τ=$τ have only finite branch lengths"] :
                [-Inf, "R2", [], [], "None of the critical points for R2 have only finite branch lengths"])
            push!(R1_output,empty_result)
        end
        # Fourth, save the information of any valid critical points that don't look like they come from a
        # submodel with zero branch lengths:
        for p in cp_list
            logL = SITE_PATTERN_DATA'log.(p)
            θ = compute_hadamard_parameters(p,τ)
            p_result = (if τ in [1,2,3]
                            ((τ ∈ 1:3) && any(≥(1-tol),θ)) ?
                                [-Inf, "R1", [], θ, "Invalid critical point for R1 with τ=$τ due to negative or length-zero branch(s)"] :
                                [logL, "R1", [τ], θ, "θ1, θ2, θ3, θ4, θ5", "binary quartet with topology τ=$τ"]
                        elseif τ == 4
                            (any(≥(1-tol),θ[1:4]) ?
                                [-Inf, "R2", [], θ, "Invalid critical point for R2 due to negative or length-zero branch(s)"] :
                                [logL, "R2", [1,2,3], θ[1:4], "θ₁, θ₂, θ₃, θ₄", [], "star quartet"])
                        end)
            # Save the result for critical point p into R1_output
            push!(R1_output,p_result)
        end
    end
    return R1_output
end




function fourLeafMLE__hopefully(SITE_PATTERN_DATA)
    "Attempt to return the global maxima."
    test_data_genericity(SITE_PATTERN_DATA) || throw(DomainError(SITE_PATTERN_DATA, "Data does not satisfy genericity conditions"))
    tmp = [compute_R1_and_R2_MLE__hopefully(SITE_PATTERN_DATA);
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

random_hadamard_edge_parameters=rand(5)
test_model=computeProbabilityVector(random_hadamard_edge_parameters,2)
N=1000
SITE_PATTERN_DATA = rand(Multinomial(N,test_model))
fourLeafMLE(SITE_PATTERN_DATA) == fourLeafMLE__hopefully(SITE_PATTERN_DATA)

function f(iterates)
    not_equal=0
    results=[]
    for _ in 1:iterates
        random_hadamard_edge_parameters=rand(5)
        model=computeProbabilityVector(random_hadamard_edge_parameters,1)
        N=1000
        SITE_PATTERN_DATA = rand(Multinomial(N,model))
        if fourLeafMLE(SITE_PATTERN_DATA)[1][3] ≠ fourLeafMLE__hopefully(SITE_PATTERN_DATA)[1][3]
            not_equal=not_equal+1
            push!(results, [model,SITE_PATTERN_DATA])
        end
    end
    return results
end
results = f(1000)

model, SITE_PATTERN_DATA = results[7]

# main comparison
fourLeafMLE(SITE_PATTERN_DATA)
fourLeafMLE__hopefully(SITE_PATTERN_DATA)

# additional stuff
listMaxima(SITE_PATTERN_DATA)
compute_hadamard_parameters(model,1)
compute_hadamard_parameters(SITE_PATTERN_DATA/1000,2)

# Check empirical 4 pt condition
B_12,B_13,B_14,B_23,B_24,B_34 = computeB(SITE_PATTERN_DATA)
a1, a2, a3 = B_12*B_34, B_13*B_24, B_14*B_23



fourLeafMLE(model)
fourLeafMLE__hopefully(model)
listMaxima(model)

# Interesting examples

## Star quartet
model = [0.22286206062815297, 0.09753860304499672, 0.17822318651244082, 0.09288433728949778, 0.11116050219014853, 0.09389481087000073, 0.09654453659241308, 0.1068919628723494]

SITE_PATTERN_DATA=[185.0, 112.0, 164.0, 99.0, 119.0, 102.0, 101.0, 118.0]


not_equal=0
τ_fail = 0
τ_fail__hopefully = 0  
iterations=1000
@time for _ in 1:iterations
    random_hadamard_edge_parameters=rand(5)
    model=computeProbabilityVector(random_hadamard_edge_parameters,1)
    N=1000
    SITE_PATTERN_DATA = rand(Multinomial(N,model))
    #print("\n",fourLeafMLE(SITE_PATTERN_DATA) == fourLeafMLE__hopefully(SITE_PATTERN_DATA))
    τ = fourLeafMLE(SITE_PATTERN_DATA)[1][3]
    τh = fourLeafMLE__hopefully(SITE_PATTERN_DATA)[1][3]
    if τ ≠ τh
        not_equal=not_equal+1
        #return [model, SITE_PATTERN_DATA]
    end
    if 1 ∉ τh
        τ_fail__hopefully = τ_fail__hopefully + 1
    end
    if 1 ∉ τ
        τ_fail = τ_fail + 1
    end
end

@time for _ in 1:1000
    random_hadamard_edge_parameters=rand(5)
    model=computeProbabilityVector(random_hadamard_edge_parameters,1)
    N=1000
    SITE_PATTERN_DATA = rand(Multinomial(N,model))
    fourLeafMLE(SITE_PATTERN_DATA)
end

not_equal
τ_fail
τ_fail__hopefully
iterations

SITE_PATTERN_DATA =  [238.0, 65.0, 67.0, 149.0, 166.0, 71.0, 60.0, 184.0]
model = [0.2398126619089718, 0.06852874116615626, 0.07074382791406153, 0.14740855135881822, 0.16730296445488688, 0.06279146652114613, 0.06246438572035237, 0.18094740095560682]*1000
fourLeafMLE(model)
fourLeafMLE(SITE_PATTERN_DATA)
fourLeafMLE__hopefully(SITE_PATTERN_DATA)


listMaxima(SITE_PATTERN_DATA)

function computeB(SITE_PATTERN_DATA)
    N=sum(SITE_PATTERN_DATA)
    M⁺(i,j,s) = [σ[i]==σ[j] for σ∈SP_const]'s
    M⁻(i,j,s) = [σ[i]!=σ[j] for σ∈SP_const]'s
    B_12 = (M⁺(1,2,SITE_PATTERN_DATA) - M⁻(1,2,SITE_PATTERN_DATA))/N
    B_13 = (M⁺(1,3,SITE_PATTERN_DATA) - M⁻(1,3,SITE_PATTERN_DATA))/N
    B_14 = (M⁺(1,4,SITE_PATTERN_DATA) - M⁻(1,4,SITE_PATTERN_DATA))/N
    B_23 = (M⁺(2,3,SITE_PATTERN_DATA) - M⁻(2,3,SITE_PATTERN_DATA))/N
    B_24 = (M⁺(2,4,SITE_PATTERN_DATA) - M⁻(2,4,SITE_PATTERN_DATA))/N
    B_34 = (M⁺(3,4,SITE_PATTERN_DATA) - M⁻(3,4,SITE_PATTERN_DATA))/N
    B = [B_12,B_13,B_14,B_23,B_24,B_34]
    return B
end

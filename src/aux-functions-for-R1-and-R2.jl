__precompile__()

###_____________________________________________________________________________
##
## Compute MLE for R1 and R2 (Binary and star trees)
##______________________________________________________________________________

# COMMENTARY: Since the critical points returned by the homotopy step do not
# necessarily correspond to valid trees (i.e., possessing finite and positive
# branch lengths), it is necessary to throw out the bad critical points. This is
# not so straightforward. The auxilliary functions `compute_x_vector' and
# `isValid' are needed to test whether any branch lengths are infinitely long,
# and this must be done before the branch lengths can be computed directly using
# `compute_hadamard_parameters' (to avoid dividing by zero). The function
# `compute_R1_and_R2_MLE' then filters out any critical points which represent
# submodels with one or more branches with zero length.

"""
   compute_x_vector(p,τ)

# Desccription:
See in-line code commentary.
"""
compute_x_vector(p,τ)=(Mt[τ]*p)[2:8]
# Commentary: the output of this function is a vector x of the form:
# x=[[θ₁*θ₂, θ₃*θ₄, θ₁*θ₃*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₂*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [θ₁*θ₃, θ₂*θ₄, θ₁*θ₂*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₃*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [θ₁*θ₄, θ₃*θ₂, θ₁*θ₂*θ₅, θ₁*θ₃*θ₅, θ₂*θ₅*θ₄, θ₃*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [θ₁*θ₂, θ₃*θ₄, θ₁*θ₃*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₂*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄]][τ]
# Note the case τ=2 and τ=3 have the same vector. This is because we chose M[4]
# to be equal to M[2]. 
#
# The significance of x is as follows. Let y = [1,x]. Then p = (1/8)*M[τ]*y.
# Since M'[τ]*M[τ] = 8*I (by definition), we can recover y (and hence x) with the
# formula y = M'[τ]*p. Then, since x has the form above, we can recover the
# θ-values from x (assuming they are nonzero).
#
# Proof: we'll use the function compute_probability_vector, which explicitly
# implements the formula from Semple and Steel, to show that p = (1/8)*M[τ]*y.
# Run the following code:
# @var θ₁ θ₂ θ₃ θ₄ θ₅ θ
# θ=[θ₁,θ₂,θ₃,θ₄,θ₅]
# Y=[[1, θ₁*θ₂, θ₃*θ₄, θ₁*θ₃*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₂*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [1, θ₁*θ₃, θ₂*θ₄, θ₁*θ₂*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₃*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [1, θ₁*θ₄, θ₃*θ₂, θ₁*θ₂*θ₅, θ₁*θ₃*θ₅, θ₂*θ₅*θ₄, θ₃*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄],
#    [1, θ₁*θ₂, θ₃*θ₄, θ₁*θ₃*θ₅, θ₁*θ₅*θ₄, θ₃*θ₂*θ₅, θ₂*θ₅*θ₄, θ₁*θ₃*θ₂*θ₄]]
# [expand.(compute_probability_vector(θ,i) -(1/8)*M[i]*Y[i]) for i in 1:4]
# These are all zero, so equality holds. This is how we know what x looks
# like, just run the code: 
# [expand.(M'[τ]*compute_probability_vector(θ,τ)) for τ in 1:4]



"""
   isValid(arg1,arg2,optional_arg3)

# Description
Check that the critical point corresponds to nonnegative or positive Hadamard
values. If these conditions are not satisfied, then the model should be thown
out because it's not an FBM critical point (but is rather an IBM critical
point).
"""
function isValid(criticalPoint,τ; tol::Float64=0)
    x = compute_x_vector(criticalPoint,τ)
    return (x[3] > tol && x[4] > tol && x[6] > tol && x[7] > tol
            && x[1] ≥ 0 && x[2] ≥ 0 && x[5] ≥ 0)
end


"""
   compute_hadamard_parameters(p,τ)

# Description
Calculate the Hadamard edge parameters θ such that the CFN process run on a
4-leaf tree with topology τ and edge parameters θ induces the site pattern
distribution given by p.
"""
function compute_hadamard_parameters(p,τ)        
    #x = compute_x_vector(BigFloat.(p),τ)
    x = compute_x_vector(p,τ)
    # We can recover the values θ₁,...,θ₅ from the x-vector using the
    # following calculation:
    if τ in 1:3
        h = sqrt.([x[1]*x[4]/x[6], x[1]*x[6]/x[4], x[2]*x[3]/x[4],
                         x[2]*x[4]/x[3], x[3]*x[6]/x[7]])
        permut  = [1:5, [1,3,2,4,5], [1,3,4,2,5]][τ] 
        return h[permut] # have to permute the elements of h so the values
                         # correspond to the correct θ-values
    elseif τ == 4
        θ₁ = sqrt(x[1]*x[4]/x[6])
        θ₂ = x[1]/θ₁
        θ₃ = x[3]/θ₁
        θ₄ = x[4]/θ₁
        θ₅ = 1
        return [θ₁,θ₂,θ₃,θ₄,θ₅]
    end
end


"""
   compute_R1_and_R2_MLE(SITE_PATTERN_DATA)

Output:
A vector of arrays, providing information about the critical points for classes
R1 and R2. In the case of R1, all three binary quartet topologies are analyzed.
"""

function compute_R1_and_R2_MLE(SITE_PATTERN_DATA)
    R1_output=[]
    tol=10e-10
    for τ in 1:4
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




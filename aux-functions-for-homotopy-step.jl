
"Use HomotopyContinuation to compute the nonnegative real critical points to the
score equations for the CFN model on a quartet tree with topology τ∈{1,2,3}
(corresponding to 13|23, 12|23, 14|23 respectively) given data. The case τ=4
corresponds to the star tree with four leaves."

function compute_FBM_critical_points(SITE_PATTERN_DATA,τ)
    @var p₁ p₂ p₃ p₄ p₅ p₆ p₇ p₈ λ₁ λ₂ λ₃ λ₄ p λ
    @var u[1:8] # data parameters (ie want to solve when u=SITE_PATTERN_DATA)
    p = [p₁,p₂,p₃,p₄,p₅,p₆,p₇,p₈]
    if τ ∈ [1,2,3]
      λ = [λ₁,λ₂,λ₃]
        quadratic_constraints=[
            [p₃*p₅-p₄*p₆+p₂*p₇+p₃*p₇+p₄*p₇+p₅*p₇+p₆*p₇+p₇^2+p₂*p₈+p₇*p₈-p₇,
             p₂*p₅+p₂*p₆+p₃*p₆+p₄*p₆+p₅*p₆+p₆^2-p₄*p₇+p₆*p₇+p₃*p₈+p₆*p₈-p₆],
            [p₃*p₅-p₄*p₆+p₂*p₇+p₃*p₇+p₄*p₇+p₅*p₇+p₆*p₇+p₇^2+p₂*p₈+p₇*p₈-p₇,
             p₂*p₃+p₂*p₄+p₃*p₄+p₄^2+p₄*p₅+p₄*p₆+p₄*p₇-p₆*p₇+p₄*p₈+p₅*p₈-p₄],
            [p₂*p₅+p₂*p₆+p₃*p₆+p₄*p₆+p₅*p₆+p₆^2-p₄*p₇+p₆*p₇+p₃*p₈+p₆*p₈-p₆,
             p₂*p₃+p₂*p₄+p₃*p₄+p₄^2+p₄*p₅+p₄*p₆+p₄*p₇-p₆*p₇+p₄*p₈+p₅*p₈-p₄]][τ]
    elseif τ == 4 # star tree case
        λ = [λ₁,λ₂,λ₃,λ₄]
        quadratic_constraints=[
            -p₁*p₆ + p₁*p₇ + p₂*p₅ - p₂*p₈ + p₃*p₈ - p₄*p₇ - p₅*p₃ + p₆*p₄,
            -p₁*p₄ + p₁*p₇ + p₂*p₃ - p₂*p₈ - p₅*p₃ + p₅*p₈ + p₆*p₄ - p₆*p₇,
            -p₁*p₄ - p₁*p₆ + p₂*p₃ + p₂*p₅ + p₃*p₈ - p₄*p₇ + p₅*p₈ - p₆*p₇]
    end
    C = [p₁+p₂+p₃+p₄+p₅+p₆+p₇+p₈-1; quadratic_constraints]
    Jacobian=differentiate(C,p) # The Jacobian of the constraints
    Λₚ=Diagonal(p)
    scoreEquations = System([u-Λₚ*Jacobian'*λ;C],variables=[p;λ],parameters=u)
    criticalPoints = real_solutions(HomotopyContinuation.solve(
        scoreEquations, solutions(start_solutions[τ]);
        target_parameters = SITE_PATTERN_DATA,
        start_parameters = start_parameters[τ], # parameter homotopy
        threading=true, show_progress = false))
    # if FBM_ONLY_NONSINGULAR == true
    #     # exclude singular critical points
    #     criticalPoints=results(criticalPoints; only_nonsingular=true)
    # end
    # Next, refine the critical points
    F = FixedParameterSystem(scoreEquations,SITE_PATTERN_DATA)
    Fcache = NewtonCache(F)
    refinedSols = refine_critical_points(F,criticalPoints,Fcache)
    # Finally, keep only solutions with nonnegative probabilities.
    filter!(isNonnegative,refinedSols)
    return refinedSols
end


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


"Check that the critical point corresponds to nonnegative or positive Hadamard
values. If these conditions are not satisfied, then the model should be thown
out because it's not an FBM critical point (but is rather an IBM critical
point)."
# function isValid(criticalPoint,τ)
#     x = compute_x_vector(criticalPoint,τ)
#     #machine_zero = eps(Float64)
#     return all(>(FBM_TOLERANCE),x)
# end
function isValid(criticalPoint,τ; tol::Float64=0)
    x = compute_x_vector(criticalPoint,τ)
    return (x[3] > tol && x[4] > tol && x[6] > tol && x[7] > tol
            && x[1] ≥ 0 && x[2] ≥ 0 && x[5] ≥ 0)
end



"Calculate the Hadamard edge parameters θ such that the CFN process run on a
4-leaf tree with topology τ and edge parameters θ induces the site pattern
distribution given by p."
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

"Extract the probability measure from the Newton-refined critical point. Note
this is just the first 8 elements of the critical point (the rest of the entries
are Lagrange multipliers that we don't need anymore). Since the critical points
are known to be real, take only the real part."
get_p(newton_result::NewtonResult)=real(getfield(newton_result,:x))[1:8]

"Refine all critical points by applying Newton's method. Output: a list of
refined probability measures."
function refine_critical_points(
    F::FixedParameterSystem, criticalPoints,
    Fcache::NewtonCache{HomotopyContinuation.MatrixWorkspace{Matrix{ComplexF64}}})
    return [get_p(newton(F,point,atol=1e-128,max_iters=100,Fcache,extended_precision=true))
            for point in criticalPoints]
end


compute_logL_of_model(SITE_PATTERN_DATA,model) = SITE_PATTERN_DATA'log.(model)

"Compute the MLE for SITE_PATTERN_DATA assuming the model is R1: FBM(τ)."
function MLEforFBM(SITE_PATTERN_DATA,τ)
    # First, compute the real critical points to the likelihood equations.
    cp_list = compute_FBM_critical_points(SITE_PATTERN_DATA,τ)
    # Second, remove those critical points which correspond to negative or zero
    # θ values (i.e. critical points not in the interior of the FBM model).
    isValid_for_τ(c)=isValid(c,τ,tol=1.0e-14)
    filter!(i->isValid_for_τ(i),cp_list)
    # valid_cp_list = filter(i->isValid(i),cp_list)
    # Third, if there are no valid critical points, return an exceptional case.
    if isempty(cp_list) return [-9223372036854775808,[],[],[],[], "No critical points for R1/R2 with topology τ=$τ"] end
    # Fourth, identify the critical point p with greatest log-likelihood.
    # (Broadcasts across valid_cp_list, but not over SITE_PATTERN_DATA):
    logL_list = compute_logL_of_model.((SITE_PATTERN_DATA,),cp_list)
    m = argmax(logL_list) # there might be more than one critical point, but
                          # only one is chosen. :( It is possible that this does
                          # not occur for generic data, but I need to check this.
    logL = logL_list[m]
    p = cp_list[m]
    # Fifth, compute the Hadamard parameters which correspond to p.
    θ = compute_hadamard_parameters(p,τ)
    # Finally, we make one final check to make sure the θ-values are in the
    # interval (0,1). If they are greater than 1, then the result is
    # nonsensical. If any θᵢ=0, then the MLE is in a submodel. Since the code
    # for newton's method won't refine more than a certain point (due to
    # hardcoded use of Float64), we sometimes get small errors, so we accept a
    # certain tolerance as to what constitutes greater than 1 or less than 0.
    # This tolerance is defined by the global variable FBM_TOLERANCE.
    # if any(<(0),θ) return [-900000000000000,[θ],[],[],[]] end
    # if FBM_ROUNDING == true
    #     any(<=(FBM_TOLERANCE),θ) && [-9223372036854775807,logL,p,θ,τ]
    # end
    if τ in [1,2,3]
        if (τ ∈ 1:3) && any(≥(1-1000000*eps()),θ)
            return [-911111111111111,[θ],[],[],[]]
        else
            return [logL,θ,[τ],"R1", "binary quartet with topology τ=$τ"]
        end
    elseif τ == 4
        if any(≥(1-1000000*eps()),θ[1:4])
            return [-911111111111111,[θ],[],[],[]]
        else
            return [logL,θ,[1,2,3],"R2", "star quartet"]
        end
    end
end


###__________________________________________________
##
## R1, R2 (Binary trees)
##___________________________________________________
function compute_R1_MLE(SITE_PATTERN_DATA)
    R1_output = []
    for τ in 1:3
        push!(R1_output, MLEforFBM(SITE_PATTERN_DATA,τ))
    end
    return R1_output
end

###__________________________________________________
##
## R2 (star tree)
##___________________________________________________
function compute_R2_MLE(SITE_PATTERN_DATA)
    R2_output= [MLEforFBM(SITE_PATTERN_DATA,4)]
    return R2_output
end

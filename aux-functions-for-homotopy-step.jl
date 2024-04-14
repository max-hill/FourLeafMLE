###_____________________________________________________________________________
##
## Auxilliary functions for the homotopy continuation step
## (used for computing the MLE for classes R1 and R2)
##______________________________________________________________________________

"""
   compute_R1_or_R2_critical_points_with_homotopy(SITE_PATTERN_DATA,τ)

# Arguments
- `SITE_PATTERN_DATA' : length 8 vector of site pattern counts
- `τ' : the tree topology (takes values 1,2,3 or 4). See function description.

# Description:
Use HomotopyContinuation to compute the nonnegative real critical points to the
score equations for the CFN model on a quartet tree with topology τ∈{1,2,3}
(corresponding to 13|23, 12|23, 14|23 respectively) given data. The case τ=4
corresponds to the star tree with four leaves.

# Output:
A list of critical points for the specified data and topology. Each critical
point is a probability measure which may or may not correspond to a tree with
finite and positive branch lengths.

"""
function compute_R1_or_R2_critical_points_with_homotopy(SITE_PATTERN_DATA,τ)
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


"""
   get_p(newton_result)

Extract the probability measure from the Newton-refined critical point. Note
this is just the first 8 elements of the critical point (the rest of the entries
are Lagrange multipliers that we don't need anymore). Since the critical points
are known to be real, take only the real part.
"""
get_p(newton_result::NewtonResult)=real(getfield(newton_result,:x))[1:8]

"Refine all critical points by applying Newton's method. Output: a list of
refined probability measures."
function refine_critical_points(
    F::FixedParameterSystem, criticalPoints,
    Fcache::NewtonCache{HomotopyContinuation.MatrixWorkspace{Matrix{ComplexF64}}})
    return [get_p(newton(F,point,atol=1e-128,max_iters=100,Fcache,extended_precision=true))
            for point in criticalPoints]
end



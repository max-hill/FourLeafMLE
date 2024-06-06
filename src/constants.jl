const SP_const = [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1),
                  (0, 1, 0, 0), (0, 1, 0, 1), (0, 1, 1, 0), (0, 1, 1, 1)]
const all_labels=[[4 3 2 1],[4 3 1 2],[4 2 3 1],[4 2 1 3],[4 1 3 2],[4 1 2 3],
                  [3 4 2 1],[3 4 1 2],[3 2 4 1],[3 2 1 4],[3 1 4 2],[3 1 2 4],
                  [2 4 3 1],[2 4 1 3],[2 3 4 1],[2 3 1 4],[2 1 4 3],[2 1 3 4],
                  [1 4 3 2],[1 4 2 3],[1 3 4 2],[1 3 2 4],[1 2 4 3],[1 2 3 4]]


# Create a list of the elements of the product space {0,1}^n ordered such that
# the elements are increasing in magnitude (when regarded as binary numbers).
S_tmp=Iterators.product(fill([-1,1],4)...)|>collect
const S=map(i->reverse(S_tmp[i]),1:16) # Ordered site pattern list. Global variable.


"Compute the matrix of coefficients used for computing likelihoods for a quartet
with topology τ. (They're Hadamard matrices. See Semple and Steel, chapter 8.6).
"
function compute_coefficient_matrix(τ)
    j₁,j₂,k₁,k₂ = [[1,2,3,4], [1,3,2,4], [1,4,2,3], [1,2,3,4]][τ] # indices depend on τ
    X = Matrix{Int}(undef,8,8)
    for i in 1:8
        row = [1,
               S[i][j₁]*S[i][j₂],
               S[i][k₁]*S[i][k₂],
               S[i][j₁]*S[i][k₁],
               S[i][j₁]*S[i][k₂],
               S[i][j₂]*S[i][k₁],
               S[i][j₂]*S[i][k₂],
               S[i][1]*S[i][2]*S[i][3]*S[i][4]]
        for j in 1:8
            X[i,j] = row[j]
        end
    end
    return X
end

# The entries of M are matrices of coefficients for computing probabilities from
# θ-values. The entries of M' are the corresponding transposes. Since M'[i]=
# 8*inv(M[i]) for i=1,2,3, the transformations in M' are useful for recovering
# Hadamard parameters from a probability measures.
const M = [compute_coefficient_matrix(τ) for τ in 1:4]
const Mt = [compute_coefficient_matrix(τ)' for τ in 1:4]



"Solve a random system. To get a start system for each type."
function compute_all_complex_critical_points(τ)
    (τ ∉ 1:4) && return "bad value of τ"
    SITE_PATTERN_DATA = randn(ComplexF64,8)
    @var  p₁ p₂ p₃ p₄ p₅ p₆ p₇ p₈ λ₁ λ₂ λ₃ λ₄ p λ
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
    scoreEquations = System([SITE_PATTERN_DATA-Λₚ*Jacobian'*λ;C],variables=[p;λ])
    result = solve(scoreEquations, threading=true,
                   show_progress = false)
    # Next check if all the critical points were recovered. If so, save the
    # results and exit the loop.
    (τ ∈ 1:3) && (length(solutions(result)) == 14) && return [SITE_PATTERN_DATA,result]
    (τ == 4) && (length(solutions(result)) == 92) && return [SITE_PATTERN_DATA,result]
end

# For each type of FBM, generate a set of starting parameters and their solution
# for use in parameter homotopy. This may take a few seconds.
tmp_list_of_pairs = [compute_all_complex_critical_points(τ) for τ in 1:4]
const start_parameters = [x[1] for x in tmp_list_of_pairs]
const start_solutions  = [x[2] for x in tmp_list_of_pairs]




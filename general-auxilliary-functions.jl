θ_index_pair(i,j)=return "θ$(min(i,j))$(max(i,j))"

# "Compute the model given a topology τ∈{1,2,3} (= 13|23,12|34,14|23,
# respectively, and vector θ=[θ₁,...,θ₅] of Hadamard branch parameters. Uses
# Formula from Semple and Steel's book on Phylogenetics. See Cor 8.6.6 on p201.
# Try @var θ[1:5]; compute_probability_vector(θ,1)"
function compute_probability_vector(θ,τ)
    j₁,j₂,k₁,k₂ = [[1,2,3,4], [1,3,2,4], [1,4,2,3], [1,2,3,4]][τ] # indices depend on τ
    [(1
      + S[i][j₁]*S[i][j₂]*θ[j₁]*θ[j₂]
      + S[i][k₁]*S[i][k₂]*θ[k₁]*θ[k₂]
      + S[i][j₁]*S[i][k₁]*θ[j₁]*θ[k₁]*θ[5]
      + S[i][j₁]*S[i][k₂]*θ[j₁]*θ[k₂]*θ[5]
      + S[i][j₂]*S[i][k₁]*θ[j₂]*θ[k₁]*θ[5]
      + S[i][j₂]*S[i][k₂]*θ[j₂]*θ[k₂]*θ[5]
      + S[i][1]*S[i][2]*S[i][3]*S[i][4]*θ[1]*θ[2]*θ[3]*θ[4])/8 for i in 1:8]
end

"Auxilliary function. Given a vector θ of Hadamard parameters and a topology
τ∈{1,2,3} (= 13|23,12|34,14|23), output a vector consisting of those products of
the Hadamard parameters (corresponding to certain paths on the quartet with
topology τ) which are used to compute the probability of the site patterns.
Example, run [V(θ,τ) for τ in 1:4]. To see the effect symbolically, run the
command @var θ[1:8]; computeProbabilityVector(θ,1) "
function V(θ,τ)
    j₁,j₂,k₁,k₂ = [[1,2,3,4], [1,3,2,4], [1,4,2,3], [1,2,3,4]][τ] # indices depend on τ
    [1, θ[j₁]*θ[j₂], θ[k₁]*θ[k₂], θ[j₁]*θ[k₁]*θ[5], θ[j₁]*θ[k₂]*θ[5],
     θ[j₂]*θ[k₁]*θ[5], θ[j₂]*θ[k₂]*θ[5], θ[1]*θ[2]*θ[3]*θ[4]]
end

"Compute the model given a topology τ∈{1,2,3} (= 13|23,12|34,14|23),
respectively, and vector θ=[θ₁,...,θ₅] of Hadamard branch parameters. Equivalent
to the (slower) compute_probaiblity_vector in general-functions.jl"
computeProbabilityVector(θ,τ)=M[τ]*V(θ,τ)/8


"Return true iff all elemnts of x are nonnegative."
isNonnegative(x) = all(>=(0),x)

"Compute the Euclidean distance between two vectors x and y."
euclidean_distance(x,y) = sqrt(sum((x-y).^2))


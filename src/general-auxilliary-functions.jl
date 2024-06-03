__precompile__()

θ_index_pair(i,j)=return "θ$(min(i,j))$(max(i,j))"

"""
    compute_probability_vector(θ,τ)

# Arguments

- `θ' : A vector θ = [θ₁,...,θ₅] Hadamard branch parameters.

- `τ' : A number in {1,2,3,4}, corresonding to the topology of the tree. Here, 1,2,3 correspond to
        quartet topologies 13|23, 12|34, and 14|23 respectively, and 4 correspond to the star tree
        topology (a four-leaf tree with internal branch length zero.

# Description

Compute the CFN probability model corresponding to the 4-leaf tree with topology τ and branch lengths θ.
Uses Formula from Semple and Steel's book on Phylogenetics. See Cor 8.6.6 on p201.

# Examples

`compute_probability_vector([.2,.3,.02,.4,.7],1)'

To display the formulas symbolically, run

`@var θ[1:5]; compute_probability_vector(θ,1)'

"""
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

"""
   V(θ,τ)

# Description

Auxilliary function. Given a vector θ of Hadamard parameters and a topology
τ∈{1,2,3} (= 13|23,12|34,14|23), output a vector consisting of those products of
the Hadamard parameters (corresponding to certain paths on the quartet with
topology τ) which are used to compute the probability of the site patterns.


# Example

`[V(θ,τ) for τ in 1:4]'

To see the effect symbolically, run

`@var θ[1:8]; computeProbabilityVector(θ,1)'

"""
function V(θ,τ)
    j₁,j₂,k₁,k₂ = [[1,2,3,4], [1,3,2,4], [1,4,2,3], [1,2,3,4]][τ] # indices depend on τ
    [1, θ[j₁]*θ[j₂], θ[k₁]*θ[k₂], θ[j₁]*θ[k₁]*θ[5], θ[j₁]*θ[k₂]*θ[5],
     θ[j₂]*θ[k₁]*θ[5], θ[j₂]*θ[k₂]*θ[5], θ[1]*θ[2]*θ[3]*θ[4]]
end

"""
   computeProbabilityVector(θ,τ)

# Description

Equivalent to the function `compute_probability_vector', but faster. So see the documentation for that
function.

"""
computeProbabilityVector(θ,τ)=M[τ]*V(θ,τ)/8


"Return true iff all elemnts of x are nonnegative."
isNonnegative(x) = all(>=(0),x)

"Compute the Euclidean distance between two vectors x and y."
euclidean_distance(x,y) = sqrt(sum((x-y).^2))


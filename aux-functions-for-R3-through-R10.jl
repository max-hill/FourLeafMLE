 ###_____________________________________________________________________________
##
## Compute MLE for R3 - R10
##______________________________________________________________________________
#=
COMMENTARY: The three-leaf theorem handles cases R6, R7, R9.

# input: site pattern data of length 8 (for a 4-leaf tree), which satisfies some genericity assumptions. Takes the form SITE_PATTERN_DATA = [xxxx, xxxy, xxyx, xxyy, xyxx, xyxy, xyyx, xyyy]

# desired output: a list of reduced trees which maximize the logL, along with attendant information about those trees. The output is an array with the following organization:

[logL, R-classification, list of compatible binary quartet topologies, branch lengths (excluding those which equal 0 or 1), branch length names, labels, description]


=#

test_D_inclusion(B)=(B[1]>0 && B[2]>0 && B[3]>0 && B[1]*B[2]<B[3] && B[2]*B[3]<B[1] && B[1]*B[3]<B[2])

function compute_B_general_full(SITE_PATTERN_DATA,labels)
    "Self contained function to compute the B=(B_jk,B_jl,B_kl) from the data.
     Input: labels [i j k l] and a frequency vector of the form
            SITE_PATTERN_DATA = [xxxx, xxxy, xxyx, xxyy, xyxx, xyxy, xyyx, xyyy]
     Output: [B_jk, B_jl, B_kl]"
    N=sum(SITE_PATTERN_DATA)
    i,j,k,l = labels
    M⁺(i,j,s) = [σ[i]==σ[j] for σ∈SP_const]'s
    M⁻(i,j,s) = [σ[i]!=σ[j] for σ∈SP_const]'s
    B_ij = (M⁺(i,j,SITE_PATTERN_DATA) - M⁻(i,j,SITE_PATTERN_DATA))/N
    B_ik = (M⁺(i,k,SITE_PATTERN_DATA) - M⁻(i,k,SITE_PATTERN_DATA))/N
    B_il = (M⁺(i,l,SITE_PATTERN_DATA) - M⁻(i,l,SITE_PATTERN_DATA))/N
    B_jk = (M⁺(j,k,SITE_PATTERN_DATA) - M⁻(j,k,SITE_PATTERN_DATA))/N
    B_jl = (M⁺(j,l,SITE_PATTERN_DATA) - M⁻(j,l,SITE_PATTERN_DATA))/N
    B_kl = (M⁺(k,l,SITE_PATTERN_DATA) - M⁻(k,l,SITE_PATTERN_DATA))/N
    return [B_ij,B_ik,B_il,B_jk,B_jl,B_kl]
end


## Compute maximizers in class R3 - DONE
# j is of degree 2. i--j--<{k,l}. k and l form a cherry with an unlabeled parent. edge parameters: h_i, h_k, h_l, h_internal
function compute_R3_MLE(SITE_PATTERN_DATA)
    R3=[[1 2 3 4], [2 1 3 4], [4 3 1 2], [3 4 1 2],
        [1 3 2 4], [2 4 1 3], [3 1 2 4], [4 2 1 3],
        [1 4 2 3], [2 3 1 4], [3 2 1 4], [4 1 2 3]] 
    R3_output=[]
    for labels in R3
        i,j,k,l = labels
        B = compute_B_general_full(SITE_PATTERN_DATA,labels)
        B_ij, B_jk, B_jl, B_kl = B[[1,4,5,6]]
        if B_ij>0 && test_D_inclusion([B_jk,B_jl,B_kl])
            θ=[0,0,0,0,0.0]
            θ[i]=B_ij
            θ[j]=1
            θ[k]=sqrt(B_jk*B_kl/B_jl)
            θ[l]=sqrt(B_jl*B_kl/B_jk)
            θ[5]=sqrt(B_jk*B_jl/B_kl)
            tol=10e-9
            if all(<(1-tol),[θ[i],θ[k],θ[l],θ[5]])
                if issetequal([i,j],[1,2]) || issetequal([k,l],[1,2])
                    τ=1
                elseif issetequal([i,j],[1,3]) || issetequal([k,l],[1,3])
                    τ=2
                elseif issetequal([i,j],[1,4]) || issetequal([k,l],[1,4])
                    τ=3
                end
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,τ))
                maximizer = [logL, "R3", [τ],
                             [B_ij, sqrt(B_jk*B_jl/B_kl), sqrt(B_jk*B_kl/B_jl), sqrt(B_jl*B_kl/B_jk)],
                             "θ$i,θ_internal,θ$k,θ$l",
                             labels, "$i -- $j --< $k, $l"]
                push!(R3_output,maximizer)
            end
        end
    end
    return R3_output
end

## Compute maximizers in class R4 - DONE
# Claw tree with labeled internal node i of degree 3. Leaves are labeled by j, k, and l. There are four elements,
# depending on the label of the internal node.
function compute_R4_MLE(SITE_PATTERN_DATA)
    R4=[[1 2 3 4], [2 1 3 4], [3 1 2 4], [4 1 2 3]]
    R4_output=[]
    for labels in R4
        i,j,k,l = labels
        B_ij, B_ik, B_il = compute_B_general_full(SITE_PATTERN_DATA,labels)[[1,2,3]]
        if B_ij>0 && B_ik>0 && B_il>0 && B_ij<1 && B_ik<1 && B_il<1
            θ=[0,0,0,0,1.0]
            θ[i]=1
            θ[j]=B_ij
            θ[k]=B_ik
            θ[l]=B_il
            if all(<(1-100000*eps()),[θ[j],θ[k],θ[l]])
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,1))
                maximizer = [logL, "R4", [1,2,3],
                             [B_ij, B_ik, B_il],
                             "θ$j,θ$k,θ$l",
                             labels, "$j $k $l tripod with internal node $i"]
                # The maximizer is (θ_j,θ_k,θ_l) = (B_ij,B_ik,B_il), provided that these are all in (0,1)
                push!(R4_output,maximizer)
            end
        end
    end
    return R4_output
end


## Compute maximizers in class R5 - DONE
# [i j k l], picture: i--j--k--l, so there are 3 edge parameters: h_i, h_l, h_{j,k}
function compute_R5_MLE(SITE_PATTERN_DATA)
    R5=[[1 2 3 4], [1 2 4 3], [2 1 3 4], [2 1 4 3],
        [1 3 2 4], [1 3 4 2], [3 1 4 2], [3 1 2 4],
        [1 4 2 3], [1 4 3 2], [4 1 3 2], [4 1 2 3]]
    R5_output=[]
    for labels in R5
        i,j,k,l = labels
        B_ij, B_jk, B_kl = compute_B_general_full(SITE_PATTERN_DATA,labels)[[1,4,6]] 
        if 0<B_ij<1 && 0<B_jk<1 && 0<B_kl<1
            θ=[0,0,0,0,0.0]
            θ[i]=B_ij
            θ[j]=1
            θ[k]=1
            θ[l]=B_kl
            θ[5]=B_jk
            if all(<(1-100000*eps()),[θ[i],θ[l],θ[5]])
                if issetequal([i,j],[1,2])
                    τ=1
                elseif issetequal([i,j],[1,3])
                    τ=2
                elseif issetequal([i,j],[1,4])
                    τ=3
                end
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,τ))
                maximizer = [logL, "R5", [τ],
                             [B_ij, B_jk, B_kl], "θ$i,"*θ_index_pair(j,k)*",θ$l",
                             labels, "$i--$j--$k--$l"]
                # The maximizer is (θ_j,θ_ik,θ_l) = (B_ij,B_ik,B_kl), provided that these are all in (0,1)
                push!(R5_output,maximizer)
            end
        end
    end
    return R5_output
end





## compute maximizers in class R6 - check this
# form [i j k l], 3-leaf tree with leaves j,k,l, and with and i infinitely far
function compute_R6_MLE(SITE_PATTERN_DATA)
    R6=[[1 2 3 4], [2 1 3 4], [3 1 2 4], [4 1 2 3]]
    R6_output=[]
    for labels in R6
        i,j,k,l = labels
        B_jk, B_jl, B_kl = compute_B_general_full(SITE_PATTERN_DATA,labels)[[4,5,6]] # B_jk, B_jl, B_kl
        if test_D_inclusion([B_jk,B_jl,B_kl])
            θ=[0,0,0,0,1.0]        
            θ[j] = sqrt(B_jk*B_jl/B_kl)
            θ[k] = sqrt(B_jk*B_kl/B_jl)
            θ[l] = sqrt(B_jl*B_kl/B_jk)
            if all(<(1-100000*eps()),[θ[j],θ[k],θ[l]])
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,1))
                maximizer=[logL, "R6",  [1,2,3], 
                           [θ[j], θ[k], θ[l]],"θ$j, θ$k, θ$l",
                           labels, "leaf $i infinitely far"]
                push!(R6_output, maximizer)
            end
        end
    end
    return R6_output
end


## compute maximizers in class R7 - DONE
# form [i,j,k,l], i infinitely far (picture: i j--k--l), (different from thesis picture)
function compute_R7_MLE(SITE_PATTERN_DATA)
    R7=[[1 2 3 4], [1 3 2 4], [1 2 4 3], [2 3 4 1], [2 4 3 1], [2 3 1 4],
        [3 2 4 1], [3 4 1 2], [3 1 2 4], [4 3 1 2], [4 1 2 3], [4 2 3 1]] 
    R7_output=[]
    for labels in R7
        i,j,k,l = labels
        B = compute_B_general_full(SITE_PATTERN_DATA,labels)[[4,5,6]]
        B_jk, B_jl, B_kl = B
        if min(B_jk,B_jl)>0
            θ=[0,0,0,0,1.0]        
            θ[j] = 1
            θ[k] = B_jk
            θ[l] = B_jl
            if all(<(1-100000*eps()),[θ[k],θ[l]])
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,1))
                maximizer = [logL, "R7", [1,2,3],
                             [B_jk,B_kl],  θ_index_pair(j,k)*","*θ_index_pair(k,l),
                             labels,  "leaf $i infinitely far and θ$j=1"]
                push!(R7_output,maximizer)
            end
        end
    end
    return R7_output
end


## Compute maximizers in class R8 - DONE
# picture: i--j k--l, with edge parameters: h_ij, h_kl
function compute_R8_MLE(SITE_PATTERN_DATA)
    R8=[[1 2 3 4], [1 3 2 4], [1 4 2 3]] 
    R8_output=[]
    for labels in R8
        i,j,k,l = labels
        B_ij, B_kl = compute_B_general_full(SITE_PATTERN_DATA,labels)[[1,6]]
        if B_ij>0 && B_kl>0 && B_ij<1 && B_kl<1
            θ=[0,0,0,0,0.0]        
            θ[i] = B_ij
            θ[j] = 1
            θ[k] = B_kl
            θ[l] = 1
            if all(<(1-100000*eps()),[θ[i],θ[k]])
                if labels == [1 2 3 4]
                    τ=1
                elseif labels == [1 3 2 4]
                    τ=2
                elseif labels == [1 4 2 3]
                    τ=3
                end
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,τ))
                maximizer = [logL, "R8", [τ],
                             [B_ij, B_kl], θ_index_pair(i,j)*","*θ_index_pair(k,l),
                             labels,  "$i--$j  $k--$l"]
                # The maximizer is (θ_ij,θ_kl) = (B_ij,B_kl), provided that these are in (0,1)
                push!(R8_output,maximizer)
            end
        end
    end
    return R8_output
end



## compute maximizers in class R9 - DONE
# picture: i j k--l, edge parameters: h_kl
function compute_R9_MLE(SITE_PATTERN_DATA)
    R9=[[1 2 3 4], [1 3 2 4], [1 4 2 3], [2 3 1 4], [2 4 1 3], [3 4 1 2]] 
    R9_output=[]
    for labels in R9
        i,j,k,l = labels
        B = compute_B_general_full(SITE_PATTERN_DATA,labels)[[4,5,6]]
        B_jk, B_jl, B_kl = B
        if B_kl>0
            θ=[0,0,0,0,1.0]        
            θ[k] = B_kl
            θ[l] = 1
            if θ[k]<1-100000*eps()
                logL = SITE_PATTERN_DATA'log.(computeProbabilityVector(θ,1))
                maximizer = [logL, "R9", [1,2,3],
                             [B_kl],  θ_index_pair(k,l),
                             labels,  "leaves $i,$j infinitely far"]
                push!(R9_output,maximizer)
            end
        end
    end
    return R9_output
end


# compute maximizers in class R10
function compute_R10_MLE(SITE_PATTERN_DATA)
    logL = SITE_PATTERN_DATA'log.(computeProbabilityVector([0,0,0,0,0],1))
    R10_output=[[logL, "R10", [1,2,3],
                 [], "θ1,θ2,θ3,θ4,θ5", "R10", "all sites independent"]]
    return R10_output
end


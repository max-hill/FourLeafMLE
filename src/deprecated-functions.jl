
###__________________________________________________
##
## Process list of MLEs
##___________________________________________________

# IMPORTANT: need to establish some criterion for determining whether the MLE is achieved by more than onegraph x type.

function round_first_entry!(list_of_arrays)
    for x in list_of_arrays
        x[1] = round(x[1], digits=10)
    end
    return list_of_arrays
end




###__________________________________________________
##
## Deprecated plotting functions
##___________________________________________________

@time begin
    m=1 # number of samples per parameter regime (replicates)
    x_values, y_values, τ1_scores, τ2_scores, τ3_scores, avg_number_of_maximizers = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    #strict_τ_records = Any[]
    sequence_length=1000
    step_size=0.01
    dx_increments, dy_increments = collect(0.01:step_size:1.49), collect(0.01:step_size:1.49) # distance increments
    for x in ProgressBar(exp.(-2*dx_increments))
        for y in exp.(-2*dy_increments)
            push!(x_values,x)
            push!(y_values,y)
            model=computeProbabilityVector([y,x,y,x,x],1)
            model=round.(model, digits = 16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            #print("x=",-log(x)/2,"; y=",-log(y)/2,"\n",model,"\n")
            number_of_maximizers=0.0
            τ_count = zeros(3)
            for i in 1:m
                SITE_PATTERN_DATA=rand(Multinomial(sequence_length,model))
                estimates=fourLeafMLE(SITE_PATTERN_DATA)
                number_of_maximizers = number_of_maximizers + length(estimates)

                # For each topology τ, identify whether there is at least one estimate compatible with τ:
                for τ in 1:3
                    if any([(τ in x[3]) for x in estimates])
                        τ_count[τ] = τ_count[τ] + 1
                    end
                end
        
            end
            push!(τ1_scores,τ_count[1])
            push!(τ2_scores,τ_count[2])
            push!(τ3_scores,τ_count[3])
            push!(avg_number_of_maximizers,number_of_maximizers/m)
        end
    end
end # 583 seconds to compute 3610 MLE estimates = about 6 estimates per second per core. Multithreading
# with 12 cores, we get about 54 per second. (This was based on a test of 44404 trees that took 8182 seconds.)


# convert branch lengths to evolutionary distances and normalize z
x_dist=-(1/2)*log.(x_values)
y_dist=-(1/2)*log.(y_values)

plot()
p1 = scatter(x_values, y_values, zcolor=τ1_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y",
             title="Estimated topology concordance with 12|34", colorbar=true,
             colorbar_title="Proportion having topology 12|34")

p2 = scatter(x_values, y_values, zcolor=τ2_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 13|24",
             colorbar=true)

p3 = scatter(x_values, y_values, zcolor=τ3_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 14|23",
             colorbar=true)



# Initialize the plot
p5 = plot(legend = :outertopright)

# Define your custom color gradient
color_mapping = Dict("multiple MLEs" => :white, "unique bdry case" => :black, "unique τ1" => :green, "unique τ2" => :red, "unique τ3" => :blue)

# Add each category to the plot separately
for category in ["unique τ1", "unique τ2", "unique τ3", "multiple MLEs", "unique bdry case"]
    mask = strict_τ_records .== category
    scatter!(p5, x_dist[mask], y_dist[mask],
             c = color_mapping[category],
             label = category,
             markershape=:square, markersize=1.2, markerstrokewidth=0)
end

# Set the remaining plot attributes
plot!(p5, aspect_ratio=:equal, xlims=(0, 1.501), ylims=(0, 1.501),
      xlabel="parameter x", ylabel="parameter y", title="Avg number of global maxima")





##
### Plot of how frequently boundary cases are inferred
sequence_length = 200
replicates_per_parameter_regime = 1
step_size = 0.02
collect(0.01:step_size:1.49)
@time begin
    x=[]
    y=[]
    z=Float64[]
    d = 0.01*(1:100) # desired distance intervals. 
    for a in exp.(-2*d)
        for b in exp.(-2*d)
            push!(x,a)
            push!(y,b)
            model=computeProbabilityVector([b,a,b,a,a],1)
            model=round.(model, digits=16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            print("a=",a,"; b=", b, "\n",model,"\n")
            color_val=0
            for i in 1:replicates_per_parameter_regime
                SITE_PATTERN_DATA=rand(Multinomial(sequence_length,model))
                if fourLeafMLE(SITE_PATTERN_DATA)[1][2] != "R1"
                    color_val=color_val+1.0
                end
            end
            push!(z,color_val)
        end
    end
end


# convert branch lengths to evolutionary distances and normalize z
x_dist=-(1/2)*log.(x)
y_dist=-(1/2)*log.(y)
z = z / maximum(z)

# plot
color_gradient = cgrad([:red, :blue])
plot()
scatter(x_dist, y_dist, zcolor=z, c = color_gradient, legend=false, markershape=:square, markersize=10, markerstrokewidth=.2)
xlabel!("a")
ylabel!("b")
title!("Boundary case frequency plot (evolutionary dist.)")

# savefig("boundary-case-frequency-",sequence_length,"bp-",replicates_per_parameter_regime," replicates")




#=
The coloring here isn't good, because there may be more than one MLE type which achieves maximal likelihood.
=#

##
### Plot with hadamard axes
x_values, y_values, z_values = Float64[], Float64[], Float64[]
x_increments,y_increments= collect(0:0.05:0.95), collect(0:0.05:0.95)
@time begin
    for x in x_increments
        for y in y_increments
            push!(x_values,x)
            push!(y_values,y)
            print("$x,$y\n")
            model=computeProbabilityVector([y,x,y,x,x],1)
            z=0
            for i in 1:5
                SITE_PATTERN_DATA=rand(Multinomial(500,model))
                if 1 in fourLeafMLE(SITE_PATTERN_DATA)[1][3] # test if the MLE is compatible with the true topology τ=1
                    z=z+1.0
                end
            end
            push!(z_values,z)
        end
    end
end

z_values = z_values / maximum(z_values) # normalize z
plot()
scatter(x_values, y_values, zcolor=z_values, c = cgrad([:white, :black]), legend=false, markershape=:square, markersize=8, markerstrokewidth=0, aspect_ratio=:equal, xlims=(-.05, 1.0), ylims=(-.05, 1.0))
xlabel!("x ")
ylabel!("y")
title!("Felsenstein zone plot (Hadamard)")


p4 = scatter(x_values, y_values, zcolor=avg_number_of_maximizers, c = cgrad([:white, :black]), legend=false,
            markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
            ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Avg number of global maxima", colorbar=true)

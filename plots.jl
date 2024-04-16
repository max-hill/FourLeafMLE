###_____________________________________________________________________________
##
## Felsenstein zone plot (distance version)
##______________________________________________________________________________

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



###_____________________________________________________________________________
##
## Felsenstein zone plot (axes in evolutionary distance)
##______________________________________________________________________________

@time begin
    m=1 # number of samples per parameter regime (replicates)
    x_values, y_values, τ1_scores, τ2_scores, τ3_scores, avg_number_of_maximizers = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    #strict_τ_records = Any[]
    sequence_length=1000
    step_size=0.01
    dx_increments, dy_increments = collect(0.01:step_size:1.49), collect(0.01:step_size:1.49) # distance increments
    for x in exp.(-2*dx_increments)
        for y in exp.(-2*dy_increments)
            push!(x_values,x)
            push!(y_values,y)
            model=computeProbabilityVector([y,x,y,x,x],1)
            model=round.(model, digits = 16) # rounding to eliminate some floating point calculation errors that result in negative probabilities.
            print("x=",-log(x)/2,"; y=",-log(y)/2,"\n",model,"\n")
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
p1 = scatter(x_dist, y_dist, zcolor=τ1_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y",
             title="Estimated topology concordance with 12|34", colorbar=true,
             colorbar_title="Proportion having topology 12|34")


p2 = scatter(x_dist, y_dist, zcolor=τ2_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 13|24",
             colorbar=true)

p3 = scatter(x_dist, y_dist, zcolor=τ3_scores, c = cgrad([:white, :black]), legend=false,
             markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
             ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Estimated topology concordance with 14|23",
             colorbar=true)

p4 = scatter(x_dist, y_dist, zcolor=avg_number_of_maximizers, c = cgrad([:white, :black]), legend=false,
            markershape=:square, markersize=1.2, markerstrokewidth=0, aspect_ratio=:equal, xlims=(0, 1.501),
            ylims=(0, 1.501), xlabel="parameter x", ylabel="parameter y", title="Avg number of global maxima", colorbar=true)



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







###_____________________________________________________________________________
##
## Colored classification plots (with m=1 for each pixel)
##______________________________________________________________________________
import Pkg; Pkg.add("LaTeXStrings")
Pkg.add("ProgressBars")
using Pipe, ProgressBars, LaTeXStrings

"""

   `generate_classification_plot_data(sequence_length, lower_x, lower_y, upper_x, upper_y, number_of_x_values, number_of_y_values; Hadamard_mode=false)'

# Description
For each pair (x,y) with x and y within some user-specified range, do the following:

  (1) generate a single datapoint of length `sequence_length' bp from a quartet tree model with
  topology 12|34 such that the branch lengths of leaves 1 and 3 are y and all other branch lengths
  are x,

  (2) infer the maximum likelihood tree estimate(s) for that datapoint,

  (3) determine which binary quartet topologies are compatible with the resulting maximum likelihood
  estimate (τ=1 is 12|34, τ=2 is 13|23, and τ=3 is 14|23). 

Then save the pair (x,y) together with the inferred compatible topology classification.

# arguments

  - `sequence_length' : length in bp of the data to be generated for each parameter regime (x,y)

  - `lower_x, lower_y, upper_x, upper_y' : these specify the lower and upper bounds for the values
                                           of x and y

  - `number_of_x_values, number_of_y_values' : the number of (evenly space) values to be chosen and
                                               sampled from the interval with specified upper and
                                               lower bounds. This determines the step size and hence
                                               the resolution of the plot; higher values will result
                                               in a plot that is less pixelated but will take longer
                                               to produce.

  - `Hadamard_mode' : The function has two modes. When `Hadamard_mode = false', x and y are
                      interpreted as being branch lengths measured in expected number of mutations
                      per site (i.e., evolutionary distance). When `Hadamard_mode = true', x and y
                      are interpreted as being Hadamard edge parameters, (i.e., which are obtained by
                      applying the function d -> exp(-2d) to the evolutionary distances).

"""
function generate_classification_plot_data(sequence_length, lower_x, lower_y, upper_x, upper_y, number_of_x_values, number_of_y_values; Hadamard_mode=false)
    number_of_x_steps, number_of_y_steps = number_of_x_values-1,number_of_y_values-1
    x_values, y_values, categorical_results = Float64[], Float64[], Any[]
    dx = (upper_x-lower_x)/number_of_x_steps
    dy = (upper_y-lower_y)/number_of_y_steps
    print("\nUsing step sizes: dx = ",dx, ", dy = ",dy,"\n")
    for x in ProgressBar(collect(lower_x:dx:upper_x))
        for y in collect(lower_y:dy:upper_y)
            edge_parameters = (Hadamard_mode ? [y,x,y,x,x] : exp.(-2*[y,x,y,x,x]))
            push!(x_values,x); push!(y_values,y)
            # Compute the estimates. Rounding step is necessary for certain edge cases where
            # approximate arithmetic results probabilities which are negative rather than zero.
            estimates = @pipe (computeProbabilityVector(edge_parameters,1)
                               |> round.(_, digits=16)
                               |> rand(Multinomial(sequence_length,_))
                               |> fourLeafMLE)
            # Identify which topologies the set of MLE points is compatible with
            τ_count = Int[any([(τ in x[3]) for x in estimates]) for τ in 1:3]
            category = ["τ=1", "τ=2", "τ=1,2", "τ=3", "τ=1,3", "τ=2,3", "Any topology"][[1 2 4]*τ_count]
            # Save the result
            append!(categorical_results,category)            
        end
    end
    return [x_values,y_values,categorical_results]
end


"""

   `make_classification_plot(sequence_length, x_values, y_values, categorical_data; marker_size=nothing, Hadamard_mode=false)'

# Description

Makes a plot using data obtained from the function `generate_classification_plot_data'.

# Arguments

  - `sequence_length' : the size of the data, in base pairs, to be generated for each parameter
                       regime (i.e., for each point on the plot)

  - `x_values, y_values, categorical_data' : the output of `generate_classification_plot_data'.

  - `marker_size' : Optional parameter. Determines the size of the squares made by the scatter plot.
                    Depending on the step size (which is automatically deduced from the vectors
                    `x_values' and `y_values'), this may need to be adjusted. Ideally, `marker_size'
                    should be large enough that the squares are visible and are almost or exactly
                    adjacent, but not so large that they overlap. If no value is supplied by the
                    user, a sensible default will be chosen, which is 171 divided by
                    max(number_of_x_values,number_of_y_values).

  - `Hadamard_mode' : Optional parameter. Specifes the plotting mode. If true, then `x_values' and
                      `y_values' are interpreted as Hadamard edge parameters. If false, then they
                      are interpreted to be branch length distances measured in expected number of
                      mutations per site. This affects the axes and labels of the plot. In order to
                      obtain good-looking plots, the value of this variable should be the same as
                      the value used when generating the data (otherwise the boxes won't be evenly
                      spaced).

# Output

A plot like thos shown in the abstract. The plot is saved as an .svg file to
"plots/hadamard-plot.svg".

"""
function make_classification_plot(sequence_length, x_values, y_values, categorical_data; marker_size=nothing, Hadamard_mode=false)
    # Determine plot bounds and step size from the data:
    lower_x, lower_y, upper_x, upper_y = [minimum.([x_values,y_values]); maximum.([x_values,y_values])]
    number_of_y_values, number_of_x_values= [sum([v[i] == v[1] for i in 1:length(v)]) for v in [x_values,y_values]]
    dx, dy = (upper_x-lower_x)/number_of_x_values, (upper_y-lower_y)/number_of_y_values
    # Attempt to pick a sensible value for marker_size if none is supplied:
    if marker_size == nothing
        marker_size = 171/maximum([number_of_x_values, number_of_y_values])
    end
    # Define the color gradient and legend labels:
    color_mapping = Dict("τ=1" => :green, "τ=2" => :red, "τ=3" => :blue, "τ=1,2" => :yellow, "τ=1,3" => :cyan,
                         "τ=2,3" => :magenta, "Any topology" => :white)
    label_mapping = Dict("τ=1" => L"\widehat{\tau}=1", "τ=2" => L"\widehat{\tau}=2", "τ=3" => L"\widehat{\tau}=3",
                         "τ=1,2" => L"\widehat{\tau}=1,2", "τ=1,3" => L"\widehat{\tau}=1,3",
                         "τ=2,3" => L"\widehat{\tau}=2,3", "Any topology" => L"\widehat{\tau} = \textrm{any}\ \textrm{topology}")
    # Begin plotting each color category separately:
    default(fontfamily="Computer Modern")
    p = plot(legend = :outertopright)
    for category in ["τ=1", "τ=2", "τ=3", "Any topology", "τ=1,2", "τ=1,3", "τ=2,3"]
        mask = categorical_results .== category
        scatter!(p, x_values[mask], y_values[mask], c = color_mapping[category], label = label_mapping[category],
                 markershape=:square, markersize=marker_size, markerstrokewidth=0)
    end
    # Set the remaining plot attributes, depending on the type of edge parameters:
    if Hadamard_mode == true
        plot!(p, aspect_ratio=:equal, xlims=(0, 1), ylims=(0, 1),
              xlabel=L"\theta_x \quad \textrm{(%$(number_of_x_values)\ values})",
              ylabel=L"\theta_y \quad \textrm{(%$(number_of_y_values)\ values})")
        savefig("plots/classification-plot-hadamard-k$(sequence_length)-resolution$(number_of_x_values)x$(number_of_y_values)-marker$(marker_size).svg")
    else
        plot!(p, aspect_ratio=:equal, xlims=(lower_x-dx, upper_x+dx), ylims=(lower_y-dy, upper_y+dy),
              xlabel=L"d_x \quad \textrm{(%$(number_of_x_values)\ values})",
              ylabel=L"d_y \quad \textrm{(%$(number_of_y_values)\ values})")
        savefig("plots/classification-plot-distance-k$(sequence_length)-resolution$(number_of_x_values)x$(number_of_y_values)-marker$(marker_size).svg")
    end
end


###_____________________________________________________________________________
##
## Example plots
##______________________________________________________________________________

## Example Hadamard Plot # 1 -- Done
sequence_length=1000
x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, .99, .99, 100, 100, Hadamard_mode=true)
make_classification_plot(sequence_length, x_values, y_values, categorical_results; Hadamard_mode=true) # uses default marker size = 171/100 = 1.71.
make_classification_plot(sequence_length, x_values, y_values, categorical_results; marker_size = 1, Hadamard_mode=true) # observe the effect of decreasing marker_size


## Example Hadamard Plot # 2 -- Done
sequence_length=1000
x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, .99, .99, 200, 200, Hadamard_mode=true)
make_classification_plot(sequence_length, x_values, y_values, categorical_results, Hadamard_mode=true);


## Example Distance Plot # 1 -- Done
sequence_length=1000
x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 100, 100)
make_classification_plot(sequence_length, x_values, y_values, categorical_results)

## Example Distance Plot # 2  -- Done
sequence_length=1000
x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 200, 200)
make_classification_plot(sequence_length, x_values, y_values, categorical_results)


## Example Distance Plot # 3
sequence_length=1000
@time x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 300, 300)
make_classification_plot(sequence_length, x_values, y_values, categorical_results)




###_____________________________________________________________________________
##
## Other ideas (TBD someday)
##______________________________________________________________________________

# Include some tests with data from Chor 2003?

###_____________________________________________________________________________
##
## Felsenstein zone plot (axes in evolutionary distance)
##______________________________________________________________________________
"""
   generate_averaged_classification_plot_data(sequence_length,
                                              m,
                                              lower_x, lower_y, upper_x, upper_y,
                                              number_of_x_values, number_of_y_values; Hadamard_mode=false)


# Description

"""
function generate_averaged_classification_plot_data(sequence_length, m, lower_x, lower_y, upper_x, upper_y, number_of_x_values, number_of_y_values; Hadamard_mode=false)
    # Intialize the vectors which will be used to store the results:
    x_values, y_values, τ1_scores, τ2_scores, τ3_scores, avg_number_of_MLEs = Float64[], Float64[], Any[], Any[], Any[], Float64[]
    number_of_x_steps, number_of_y_steps = number_of_x_values-1,number_of_y_values-1
    # Compute and display step size:
    dx = (upper_x-lower_x)/number_of_x_steps
    dy = (upper_y-lower_y)/number_of_y_steps
    print("\nUsing step sizes: dx = ",dx, ", dy = ",dy,"\n")
    # Main loop (through all values of x and y in the specified range):
    for x in ProgressBar(collect(lower_x:dx:upper_x))
        for y in collect(lower_y:dy:upper_y)
            # Save the x and y coordinates:
            push!(x_values,x); push!(y_values,y)
            # Construct the probability distribution (model) corresponding to pair (x,y):
            model = @pipe ((Hadamard_mode ? [y,x,y,x,x] : exp.(-2*[y,x,y,x,x]))
                           |>computeProbabilityVector(_,1)
                           |>round.(_, digits=16))
            # Intialize variables to save the MLE results for this model:
            topology_counts = zeros(3)
            number_of_MLEs = 0
            for i in 1:m
                # Sample from the model and compute an MLE from the sample:
                estimates = rand(Multinomial(sequence_length,model)) |> fourLeafMLE
                # Identify which topologies the set of ML points is compatible with:
                topology_counts = topology_counts + Int[any([(τ in x[3]) for x in estimates]) for τ in 1:3]
                # Keep a running total of the number of MLEs for samples corresponding to (x,y):
                number_of_MLEs = number_of_MLEs + length(estimates)
            end
            # Save the values (x,y) and the results for the model with those branch lengths:
            push!(τ1_scores,topology_counts[1]/m)
            push!(τ2_scores,topology_counts[2]/m)
            push!(τ3_scores,topology_counts[3]/m)
            push!(avg_number_of_MLEs,number_of_MLEs/m)
        end
    end
    return [x_values,y_values,τ1_scores,τ2_scores,τ3_scores,avg_number_of_MLEs]
end

"""
  make_averaged_classification_plot(sequence_length, 
                                           x_values, y_values, scores;
                                           topology=nothing,
                                           marker_size=nothing,
                                           Hadamard_mode=false)

# Description

Using data from generate_averaged_classification_plot_data, makes a Felsenstein plot with points
representing the propotion of samples whose MLE is concordant with the specified topology.

# Output

Saves a plot in the directory `plots/'
"""
function make_averaged_classification_plot(sequence_length,
                                           m,
                                           x_values, y_values, scores;
                                           topology=nothing,
                                           marker_size=nothing,
                                           Hadamard_mode=false)
    # Determine plot bounds and step size from the data:
    lower_x, lower_y, upper_x, upper_y = [minimum.([x_values,y_values]); maximum.([x_values,y_values])]
    number_of_y_values, number_of_x_values= [sum([v[i] == v[1] for i in 1:length(v)]) for v in [x_values,y_values]]
    dx, dy = (upper_x-lower_x)/number_of_x_values, (upper_y-lower_y)/number_of_y_values
    # Attempt to pick a sensible value for marker_size if none is supplied:
    if marker_size == nothing
        marker_size = 171/maximum([number_of_x_values, number_of_y_values])
    end
    p = scatter(x_values, y_values, zcolor=scores, c = cgrad([:white, :black]), legend=false,
                markershape=:square, markersize=marker_size, markerstrokewidth=0, aspect_ratio=:equal,
                colorbar=true)
    # Identify topology name:
    topology_name=["12|34", "13|24", "14|23"][topology]
    plot!(p, title="MLE concordance with $(topology_name) (m=$(m))",
          colorbar_title="Proportion of estimates consistent with topology $(topology_name)")
    # Begin plotting each color category separately:
    default(fontfamily="Computer Modern")
    # Set the remaining plot attributes, depending on the type of edge parameters:
    if Hadamard_mode == true
        plot!(p, xlims=(0, 1), ylims=(0, 1),
              xlabel=L"\theta_x \quad \textrm{(%$(number_of_x_values)\ values})",
              ylabel=L"\theta_y \quad \textrm{(%$(number_of_y_values)\ values})")
        savefig("plots/averaged-classification-plot-hadamard-τ$(topology)-k$(sequence_length)-resolution$(number_of_x_values)x$(number_of_y_values)-marker$(marker_size).svg")
    else
        plot!(p, xlims=(lower_x-dx, upper_x+dx), ylims=(lower_y-dy, upper_y+dy),
              xlabel=L"d_x \quad \textrm{(%$(number_of_x_values)\ values})",
              ylabel=L"d_y \quad \textrm{(%$(number_of_y_values)\ values})")
        savefig("plots/averaged-classification-plot-distance-τ$(topology)--k$(sequence_length)-resolution$(number_of_x_values)x$(number_of_y_values)-marker$(marker_size).svg")
    end   
end

# # Examples [FIX THESE]
# sequence_length=1000 # in bp
# m=20 # number of replicates per point
# @time x_values, y_values, τ1_scores, τ2_scores, τ3_scores, avg_number_of_MLEs = generate_averaged_classification_plot_data(1000, 0.01, 0.01, 1.5, 1.5, 100, 100, Hadamard_mode=false)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ1_scores; topology=1, Hadamard_mode=false)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ2_scores, topology=2; Hadamard_mode=false)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ3_scores, topology=3; Hadamard_mode=false)

# x_values, y_values, τ1_scores, τ2_scores, τ3_scores, avg_number_of_MLEs = generate_averaged_classification_plot_data(1000, 0.01, 0.01, 0.99, 0.99, 100, 100, Hadamard_mode=true)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ1_scores, topology=1; Hadamard_mode=true)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ2_scores, topology=2; Hadamard_mode=true)
# make_averaged_classification_plot(sequence_length, m, x_values, y_values, τ3_scores, topology=3; Hadamard_mode=true)








###_____________________________________________________________________________
##
## Colored classification plots (with m=1 for each pixel)
##______________________________________________________________________________

"""

   `generate_classification_plot_data(sequence_length,
                                      lower_x, lower_y, upper_x, upper_y,
                                      number_of_x_values, number_of_y_values;
                                      Hadamard_mode=false)'

# Description
For each pair (x,y) with x and y within some user-specified range, do the following:

  (1) generate a single datapoint of length `sequence_length' bp from a quartet tree model with
  topology 12|34 such that the branch lengths of leaves 1 and 3 are y and all other branch lengths
  are x,

  (2) infer the maximum likelihood tree estimate(s) for that datapoint,

  (3) determine which binary quartet topologies are compatible with the resulting maximum likelihood
  estimate (τ=1 is 12|34, τ=2 is 13|23, and τ=3 is 14|23). 

Then save the pair (x,y) together with the inferred compatible topology classification.

# Arguments

  - `sequence_length' : length in bp of the data to be generated for each parameter regime (x,y)

  - `lower_x, lower_y, upper_x, upper_y' : these specify the lower and upper bounds for the values
                                           of x and y

  - `number_of_x_values, number_of_y_values' : the number of (evenly space) values to be chosen and
                                               sampled from the interval with specified upper and
                                               lower bounds. This determines the step size and hence
                                               the resolution of the plot; higher values will result
                                               in a plot that is less pixelated but will take longer
                                               to produce.

  - `Hadamard_mode' : Optional parameter which determines which of 2 modes the function is run in.
                      When `Hadamard_mode = false', x and y are interpreted as being branch lengths
                      measured in expected number of mutations per site (i.e., evolutionary
                      distance). When `Hadamard_mode = true', x and y are interpreted as being
                      Hadamard edge parameters, (i.e., which are obtained by applying the function d
                      -> exp(-2d) to the evolutionary distances).

"""
function generate_classification_plot_data(sequence_length, lower_x, lower_y, upper_x, upper_y,
                                           number_of_x_values, number_of_y_values;
                                           Hadamard_mode=false)
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

   `make_classification_plot(sequence_length,
                             x_values, y_values, categorical_data;
                             marker_size=nothing,
                             Hadamard_mode=false)'

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
function make_classification_plot(sequence_length,
                                  x_values, y_values, categorical_data;
                                  marker_size=nothing, Hadamard_mode=false)
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
    label_mapping = Dict("τ=1" => L"\widehat{\tau}=12|34",
                         "τ=2" => L"\widehat{\tau}=13|24",
                         "τ=3" => L"\widehat{\tau}=14|23",
                         "τ=1,2" => L"\widehat{\tau}=12|34\ \textrm{or}\ 13|24",
                         "τ=1,3" => L"\widehat{\tau}=12|34\ \textrm{or}\ 14|23",
                         "τ=2,3" => L"\widehat{\tau}=13|24\ \textrm{or}\ 14|23",
                         "Any topology" => L"\widehat{\tau} = \textrm{any}\ \textrm{topology}")
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

# ## Example Hadamard Plot # 1 -- Done
# sequence_length=1000
# x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, .99, .99, 100, 100, Hadamard_mode=true)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results; Hadamard_mode=true) # uses default marker size = 171/100 = 1.71.
# make_classification_plot(sequence_length, x_values, y_values, categorical_results; marker_size = 1, Hadamard_mode=true) # observe the effect of decreasing marker_size


# ## Example Hadamard Plot # 2 -- Done
# sequence_length=1000
# x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, .99, .99, 200, 200, Hadamard_mode=true)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results, Hadamard_mode=true);

# ## Example Hadamard Plot # 3 
# sequence_length=1000
# x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, .99, .99, 300, 300, Hadamard_mode=true)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results, Hadamard_mode=true);


# ## Example Distance Plot # 1 -- Done
# sequence_length=1000
# x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 100, 100)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results)

# ## Example Distance Plot # 2  -- Done
# sequence_length=1000
# x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 200, 200)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results)


# ## Example Distance Plot # 3 -- Done
# sequence_length=1000
# @time x_values, y_values, categorical_results = generate_classification_plot_data(sequence_length, .01, .01, 1.5, 1.5, 300, 300)
# make_classification_plot(sequence_length, x_values, y_values, categorical_results)




###_____________________________________________________________________________
##
## Other ideas (TBD someday)
##______________________________________________________________________________

# Include some tests with data from Chor 2003?


###_____________________________________________________________________________
##
## Sample complexity plot 
##______________________________________________________________________________

function generate_sample_complexity_plot_data(lower_f, upper_f, lower_k, upper_k,
                                              number_of_f_values, number_of_k_values,
                                              leaf_length;
                                              Hadamard_mode=false)
    number_of_f_steps, number_of_k_steps = number_of_f_values-1,number_of_k_values-1
    f_values, k_values, categorical_results = Float64[], Float64[], Any[]
    df = (upper_f-lower_f)/number_of_f_steps
    dk = (upper_k-lower_k)/number_of_k_steps
    L=leaf_length
    print("\nUsing step sizes: df = ",df, ", dk = ",dk,"\n")
    for f in ProgressBar(collect(lower_f:df:upper_f))
        for k in collect(lower_k:dk:upper_k)
            edge_parameters = (Hadamard_mode ? [L,L,L,L,f] : exp.(-2*[L,L,L,L,f]))
            sequence_length = Int(round(k, digits=0))
            push!(f_values,f); push!(k_values,sequence_length)
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
    return [f_values,k_values,categorical_results]
end

function make_sample_complexity_plot(f_values, k_values, categorical_data, f_scaling; 
                                  marker_size=nothing, Hadamard_mode=false)
    # Determine plot bounds and step size from the data:
    lower_f, lower_k, upper_f, upper_k = [minimum.([f_values,k_values]); maximum.([f_values,k_values])]
    number_of_k_values, number_of_f_values= [sum([v[i] == v[1] for i in 1:length(v)]) for v in [f_values,k_values]]
    df, dk = (upper_f-lower_f)/number_of_f_values, (upper_k-lower_k)/number_of_k_values
    # Attempt to pick a sensible value for marker_size if none is supplied:
    if marker_size == nothing
        marker_size = 171/maximum([number_of_f_values, number_of_k_values])
    end
    # Define the color gradient and legend labels:
    color_mapping = Dict("τ=1" => :green, "τ=2" => :red, "τ=3" => :blue, "τ=1,2" => :yellow, "τ=1,3" => :cyan,
                         "τ=2,3" => :magenta, "Any topology" => :white)
    label_mapping = Dict("τ=1" => L"\widehat{\tau}=12|34",
                         "τ=2" => L"\widehat{\tau}=13|24",
                         "τ=3" => L"\widehat{\tau}=14|23",
                         "τ=1,2" => L"\widehat{\tau}=12|34\ \textrm{or}\ 13|24",
                         "τ=1,3" => L"\widehat{\tau}=12|34\ \textrm{or}\ 14|23",
                         "τ=2,3" => L"\widehat{\tau}=13|24\ \textrm{or}\ 14|23",
                         "Any topology" => L"\widehat{\tau} = \textrm{any}\ \textrm{topology}")
    # Begin plotting each color category separately:
    default(fontfamily="Computer Modern")
    p = plot(legend = :outertopright)
    for category in ["τ=1", "τ=2", "τ=3", "Any topology", "τ=1,2", "τ=1,3", "τ=2,3"]
        mask = categorical_results .== category
        scatter!(p, f_scaling * f_values[mask], k_values[mask], c = color_mapping[category], label = label_mapping[category],
                 markershape=:square, markersize=marker_size, markerstrokewidth=0)
    end
    # Set the remaining plot attributes, depending on the type of edge parameters:
    if Hadamard_mode == true
        plot!(p, aspect_ratio=:equal, xlims=(0, 1), ylims=(lower_k, upper_k),
              xlabel=L"\textrm{internal\ branch\ length\ }f \quad \textrm{(%$(number_of_f_values)\ values})",
              ylabel=L"\textrm{sequence\ length\ }k \quad \textrm{(%$(number_of_k_values)\ values})")
        savefig("sample-complexity-plot-hadamard.svg")
    else
        plot!(p, aspect_ratio=:equal,
              xlims=((lower_f-df)*f_scaling, upper_f*f_scaling), ylims=(lower_k-dk, upper_k+dk),
              xlabel=L"\textrm{internal\ branch\ length\ }f \quad \textrm{(%$(number_of_f_values)\ values})",
              ylabel=L"\textrm{sequence\ length\ }k \quad \textrm{(%$(number_of_k_values)\ values})")
        # savefig("sample-complexity-plot-distance.svg")
    end
end


# f_values, k_values, categorical_results = generate_sample_complexity_plot_data(.0001, .1, 100, 10000,
#                                                                                50, 50,
#                                                                                .2)
# f_scaling = 100000
# make_sample_complexity_plot(f_values, k_values, categorical_results, f_scaling, marker_size=3)
# f_values

# xticks!(0:.1*f_scaling:10, [string(x/f_scaling) for x in 0:.1*f_scaling:10])
# k_values


# I could run this procedure many times and then blend the colors of the plots

# I need to determine how the picture depends on the leaf length (Also: try choosing random leaf lengths within some biologically plausible region, and blending the results)

# One likelihood correspondend for one topology. A different correspondence for a different topology. Three sets of coordinates. p1's, p2's and data u. When p1's are equal to p2's, what set of data does this occur for?

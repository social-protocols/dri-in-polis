using Plots
using CSV
using DataFrames
using Statistics
using Test

include("prepare-data.jl")
include("dri-calculations.jl")
include("dri-plot.jl")

data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)

function calculate_case_component_scores(case)
	ICs = []
	for stage in 1:2
		data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
		F_df, Pref_df = prepare_data(data_filtered)
		IC = calculate_IC(F_df, Pref_df)
		push!(ICs, IC)
	end

	pre = [p for p in zip(ICs[1][!,:Q], ICs[1][!,:R])]
	post = [p for p in zip(ICs[2][!,:Q], ICs[2][!,:R])]

	cscores = [component_scores(pair[1], pair[2]) for pair in zip(pre, post)]

	h_mean = mean(s[1] for s in cscores)
	v_mean = mean(s[2] for s in cscores)
	h_var = var(s[1] for s in cscores)
	v_var = var(s[2] for s in cscores)

	return (h_mean, v_mean, h_var, v_var)
end

function component_scores_plot(case, case_name)
	# Create a layout with a title plot and two subplots
	title_plot = plot(title = "Component Scores Case $case ($case_name)", grid = false, showaxis = false, bottom_margin = -50Plots.px)

	# Calculate scores
	h_mean, v_mean, h_var, v_var = calculate_case_component_scores(case)
	
	# Get the raw component scores for scatter plot
	ICs = []
	for stage in 1:2
		data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
		F_df, Pref_df = prepare_data(data_filtered)
		IC = calculate_IC(F_df, Pref_df)
		push!(ICs, IC)
	end

	pre = [p for p in zip(ICs[1][!,:Q], ICs[1][!,:R])]
	post = [p for p in zip(ICs[2][!,:Q], ICs[2][!,:R])]
	cscores = [component_scores(pair[1], pair[2]) for pair in zip(pre, post)]

	# Create the scatter plot
	p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)
	
	# Plot the scatter points
	scatter!(p1, 
		[s[1] for s in cscores], 
		[s[2] for s in cscores],
		xlims=(-1,1),
		ylims=(-1,1),
		color=:blue,
		alpha=0.7,
		label="Component Scores"
	)
	
	# Add a diagonal line
	plot!(p1, [-1,1], [-1,1], 
		color=:black, 
		linestyle=:dash, 
		label="Diagonal"
	)
	
	# Add axes at 0
	vline!(p1, [0], color=:black, linestyle=:dash, label="")
	hline!(p1, [0], color=:black, linestyle=:dash, label="")
	
	# Add mean point with error bars
	scatter!(p1, [h_mean], [v_mean], 
		color=:red, 
		markersize=8, 
		label="Mean Score",
		xerr=[sqrt(h_var)],
		yerr=[sqrt(v_var)],
	)
	
	# Add error bars for both dimensions
	# Plots.errorbar!(p1, 
	# 	[h_mean], [v_mean],
	# 	xerr=[sqrt(h_var)],
	# 	yerr=[sqrt(v_var)],
	# 	color=:red,
	# 	label="Standard Deviation"
	# )
	
	# Add labels
	xlabel!(p1, "Horizontal Score")
	ylabel!(p1, "Vertical Score")

	# Combine all plots with the layout
	p = plot(title_plot, p1, layout = @layout([A{0.1h}; B]), size=(1000,500))

	outdir = "local-output/component-scores/"
	mkpath(outdir)

	savefig(p, "$outdir/case-$case-component-scores.png")
end






# # Find first counter-clockwise delta
# function compute_deltas(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64})
# 	# a = (.5, .1)
# 	# b = (-.1, .2)
# 	signs = sign.(b .- a)
# 	# if signs == (1, 1)
# 	# 	x
# 	# if signs == (1, -1)
# 	# 	y
# 	# if signs == (-1, 1)
# 	# 	y
# 	# if signs == (-1, -1)
# 	# 	x
# 	if signs[1] == signs[2]
# 		# horizontal first
# 	      hDistanceDiagonal = abs(a[1] - a[2])

# 	      midpoint = (b[1], a[2]) 
# 	      newHDistanceDiagonal = abs(midpoint[1] - midpoint[2])
# 	      hDelta = newHDistanceDiagonal - hDistanceDiagonal

# 	      # then vertical
# 	      vDistanceDiagonal = abs(midpoint[1] - midpoint[2])
# 	      newVDistanceDiagonal = abs(b[1] - b[2])
# 	      vDelta = newVDistanceDiagonal - vDistanceDiagonal

# 	      return (hDelta, vDelta)
# 	else
# 		# vertical first
# 	      vDistanceDiagonal = abs(a[2] - a[1])

# 	      midpoint = (a[1], b[2])
# 	      newVDistanceDiagonal = abs(midpoint[1] - midpoint[2])
# 	      vDelta = newVDistanceDiagonal - vDistanceDiagonal

# 	      hDistanceDiagonal = abs(midpoint[1] - midpoint[2])
# 	      newHDistanceDiagonal = abs(b[1] - b[2])
# 	      hDelta = newHDistanceDiagonal - hDistanceDiagonal


# 		return (hDelta, vDelta)
# 	end

# end



	# scatter([a,b];xlims=(0,1), ylims=(0,1))
	# # Add a diagonal line
	# plot!([-1,1], [-1,1], 
	# 	color=:black, 
	# 	linestyle=:dash, 
	# 	label="Diagonal"
	# )


function compute_deltas(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64})

	# delta_v_new = abs(b[2] - b[1]) - abs(a[2] - b[1])
	# delta_v_old = abs(b[2] - a[1]) - abs(a[2] - a[1])

	# delta_h_new = abs(b[1] - b[2]) - abs(a[1] - b[2])
	# delta_h_old = abs(b[1] - a[2]) - abs(a[1] - a[2])

	# return (mean([delta_h_old, delta_h_new]), mean([delta_v_old, delta_v_new]))


	# delta_h = abs(b[1] - b[2]) - abs(a[1] - b[2])	

	# delta_v = abs(b[1] - b[2]) - abs(b[1] - a[2])	



	delta_h = abs(b[1] - b[2]) - abs(a[1] - b[2])	
	delta_v = abs(b[1] - b[2]) - abs(b[1] - a[2])


	return (delta_h, delta_v)
end


function component_scores(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64})
	ds = compute_deltas(a, b)
	return ds .* (-1,-1)
end

function test_component_scores(a::Tuple{Float64, Float64}, b::Tuple{Float64, Float64}, expected::Tuple{Float64, Float64})
	scores = component_scores(a, b)

	@show a
	@show b
	@show scores
	@show expected

	@test scores[1] ≈ expected[1]
	@test scores[2] ≈ expected[2]
end




function test() 
	test_component_scores((.9, .1), (.9,.2), (0.0, .1)) # Towards diagonal vertically
	test_component_scores((.9, .2), (.9,.1), (0.0, -.1)) # Away from diagonal vertically
	test_component_scores((.9, .1), (.8,.1), (-.1, 0.0)) # Towards diagonal horizontally
	test_component_scores((.8, .1), (.9,.1), (.1, 0.0)) # Away from diagonal horizontally

	# Top-left Triangle
	test_component_scores((.1, .9), (.1,.8), (0.0, .1)) # Towards diagonal vertically
	test_component_scores((.1, .8), (.1,.9), (0.0, -.1)) # Away from diagonal vertically
	test_component_scores((.1, .9), (.2,.9), (-.1, 0.0)) # Towards diagonal horizontally
	test_component_scores((.2, .9), (.1,.9), (.1, 0.0)) # Away from diagonal horizontally

	# Crossing symmetric
	test_component_scores((1.0, 0.0), (0.0, 1.0), (0.0, 0.0)) # Up and left
	test_component_scores((0.0, 1.0), (1.0, 0.0), (0.0, 0.0)) # Down and right


	# Crossing asymmetric
	test_component_scores((.5, 0.0), (0.0, .2), (-.3, 0.0)) 
		# counter-clockwise 
		# up .2 towards diagonal to (.5,.2)
		#	vertical score is .2
		# left .5 across diagonal to (0, .2)
			# diagonal at 0.2. Distance from diagonal at (.5,.2) was .3. At (0,.2) is .2. Change is -.1
			# horizontal score is -.1

	test_component_scores((0.0, .2), (.5, 0.0), (.3, 0.0)) # Down and right

end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
	# Get unique cases
	cases = unique(data.CaseID)
	
	# Create results DataFrame
	results = DataFrame(
		case = Float64[],
		CaseName = String[],
		h_mean = Float64[],
		v_mean = Float64[],
		h_var = Float64[],
		v_var = Float64[],
		difference = Float64[],  # Add difference column
		dri_delta = Float64[]    # Add DRI delta (post - pre)
	)
	
	# Process each case
	for case in cases
		case_data = filter(r -> r.CaseID == case, data)
		case_name = case_data[1, :Case]
		
		# Calculate component scores
		h_mean, v_mean, h_var, v_var = calculate_case_component_scores(case)
		
		# Calculate difference
		difference = v_mean - h_mean
		
		# Calculate DRI delta
		dri_pre = 0.0
		dri_post = 0.0
		for stage in 1:2
			data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
			F_df, Pref_df = prepare_data(data_filtered)
			IC = calculate_IC(F_df, Pref_df)
			dri = calculate_dri(IC)
			if stage == 1
				dri_pre = dri
			else
				dri_post = dri
			end
		end
		dri_delta = dri_post - dri_pre
		
		# Add to results
		push!(results, (
			case,
			case_name,
			h_mean,
			v_mean,
			h_var,
			v_var,
			difference,
			dri_delta
		))
		
		# Generate plot for each case
		component_scores_plot(case, case_name)
	end
	
	# Save results to CSV
	mkpath("local-output/component-scores")
	CSV.write("local-output/component-scores/component-scores-results.csv", results)
	
	# Create horizontal bar chart comparing all cases
	# Sort results by case number in reverse order
    sort!(results, [:case], by=x -> x, rev=true)
	
	# Find the range of scores
	min_score = minimum([minimum(results.h_mean), minimum(results.v_mean)])
	max_score = maximum([maximum(results.h_mean), maximum(results.v_mean)])
	
	# Add some padding to the range
	score_range = max_score - min_score
	padding = 0.1 * score_range
	x_min = min_score - padding
	x_max = max_score + padding
	
	# Create the plot
	p = plot(
		size=(1000, 400 + 30 * nrow(results)),
		left_margin=60Plots.mm,
		bottom_margin=10Plots.mm,
		top_margin=10Plots.mm,
		right_margin=10Plots.mm,
		legend=:topright,
		title="Component Scores by Case",
		xlims=(x_min, x_max),
		yticks=[]
	)
	
	# Prepare data for side-by-side bars with spacing
	y_positions = 1:1.5:(1.5 * nrow(results))
	bar_width = 0.2
	spacing = bar_width  # Make spacing equal to bar width so bars touch
	
	# Plot each type of bar separately with offset
	bar!(p,
		y_positions .- spacing/2,  # Adjust offset to center the pair of bars
		results.h_mean,
		bar_width=bar_width,
		label="Horizontal Score",
		color=:blue,
		alpha=0.7,
		orientation=:h
	)
	
	bar!(p,
		y_positions .+ spacing/2,  # Adjust offset to center the pair of bars
		results.v_mean,
		bar_width=bar_width,
		label="Vertical Score",
		color=:red,
		alpha=0.7,
		orientation=:h
	)
	
	# Add difference score bar (vertical - horizontal)
	difference_scores = results.v_mean .- results.h_mean
	bar!(p,
		y_positions .+ spacing * 1.5,  # Place to the right of the pair
		difference_scores,
		bar_width=bar_width,
		label="Difference (Vertical - Horizontal)",
		color=:purple,
		alpha=0.7,
		orientation=:h
	)
	
	# Add DRI delta bar
	bar!(p,
		y_positions .+ spacing * 2.5,  # Place to the right of the difference bar
		results.dri_delta,
		bar_width=bar_width,
		label="DRI Delta (Post - Pre)",
		color=:green,
		alpha=0.7,
		orientation=:h
	)
	
	# Add case names as text annotations
	for (i, case_name) in enumerate(results.CaseName)
		annotate!(p, 
			x_min - 0.01 * score_range,
			y_positions[i],
			text(case_name, 8, :right)
		)
	end
	
	# Add a vertical line at x=0
	# vline!(p, [0], color=:black, linestyle=:dash, label="")
	
	# Add labels
	xlabel!(p, "Component Score")
	ylabel!(p, "")
	
	# Save the plot
	savefig(p, "local-output/component-scores/component-scores-comparison.png")
end


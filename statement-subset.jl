using Plots

include("prepare-data.jl")
include("dri-calculations.jl")
include("dri-plot.jl")


function statement_selection_subplot!(p, case, stage; direction=1)
	data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
	# Calculate original DRI
	F_df, Pref_df = prepare_data(data_filtered)
	IC_observed = calculate_IC(F_df, Pref_df)
	observed_dri = calculate_dri(IC_observed)

	vars = var(Matrix(F_df), dims=2)      

	# indices of the half of elements in vars highest variances
	# idx = sortperm(vec(vars))[div(length(vars),5)*4+1:end]
	if direction !== 0
		idx = sortperm(vec(vars) * direction)[1:div(length(vars),4)]
	else
		idx = 1:nrow(F_df)
	end

	IC_observed = calculate_IC(F_df[idx,:], Pref_df)
	observed_dri = calculate_dri(IC_observed)


	title = if stage == 1
		"Pre-Deliberation"
	else
		"Post-Deliberation"
	end

	dri_plot(p,
	    IC_observed[!, :Q], IC_observed[!, :R],
	    title,
	    observed_dri
	)
end

data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)

function statement_selection_sidebyside_plot(case; direction=1)

	title = if direction == 1
		"Least Controversial 25% of Considerations"
	elseif direction == -1
		"Most Controversial 25% of Considerations"
	else
		"All Considerations"
	end

	filename = if direction == 1
		"bottom-quartile"
	elseif direction == -1
		"top-quartile"
	else
		"all-considerations"
	end



    # Create a layout with a title plot and two subplots
    title_plot = plot(title = "DRI Plots Case $case: $title", grid = false, showaxis = false, bottom_margin = -50Plots.px)


    # Create the subplots
    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)
    p2 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

	statement_selection_subplot!(p1[1], case, 1; direction=direction)
	statement_selection_subplot!(p2[1], case, 2; direction=direction)

    # Combine all plots with the layout, using more space for the title section
    p = plot(title_plot, p1, p2, layout = @layout([A{0.1h}; [B C]]), size=(1000,500))


	outdir = "local-output/statement-subset/"
	mkpath(outdir)

	savefig(p, "$outdir/$filename-considerations-case-$case.png")

end

statement_selection_sidebyside_plot(case, direction=-1)
statement_selection_sidebyside_plot(case, direction=0)
statement_selection_sidebyside_plot(case, direction=1)




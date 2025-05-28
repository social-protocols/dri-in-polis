using MultivariateStats
using Statistics
using LinearAlgebra
using StatsBase # For rank correlations
using CSV
using DataFrames
using FactorRotations      
using Plots
using Distances
using GLMNet
using Statistics # For calculating metrics if needed
using Lasso, LinearAlgebra
using Random

include("prepare-data.jl")
include("dri-calculations.jl")
include("dri-plot.jl")
include("dri-calculations.jl")

shuffle_method = :individuals
# dri_variant = :pearson
dri_variant = :standard

function calculate_dri_variant(IC)
    return dri_variant == :standard ? calculate_dri(IC) : dri_variant == :pearson ? calculate_pdri(IC) : error("Unknown DRI variant $dri_variant")
end

# Pearson's DRI is like DRI except that it uses Pearson's correlation instead of scaled distance from the diagonal
function calculate_pdri(IC)
    return cor(IC.Q, IC.R)
end

"""
    create_dri_histograms(case, case_name, DRIs, DRIs_permuted, p_values)

Creates three histograms showing the distribution of DRI values from resampling:
1. Pre-deliberation DRI
2. Post-deliberation DRI
3. Delta DRI (Post - Pre)

Each histogram includes a vertical line showing the actual observed DRI value and the p-value.
"""
function create_dri_histograms(case, case_name, DRIs, DRIs_permuted, p_values)
    # Create main title plot
    title_plot = plot(
        title = "Null Distribution (DRI with Permuted Preferences)\nCase $case ($case_name)",
        grid = false,
        showaxis = false,
        bottom_margin = -50Plots.px
    )

    # Create the subplots
    p1 = plot(layout=(1,1), size=(500,400), dpi=100, margin=5Plots.mm)
    p2 = plot(layout=(1,1), size=(500,400), dpi=100, margin=5Plots.mm)
    p3 = plot(layout=(1,1), size=(500,400), dpi=100, margin=5Plots.mm)

    # Extract DRI values for each stage
    pre_dris = [d[1] for d in DRIs_permuted]
    post_dris = [d[2] for d in DRIs_permuted]
    delta_dris = [d[3] for d in DRIs_permuted]

    # Plot each histogram
    for (i, (p, values, observed, title, p_val)) in enumerate([
        (p1, pre_dris, DRIs[1], "Pre-Deliberation DRI", p_values[1]),
        (p2, post_dris, DRIs[2], "Post-Deliberation DRI", p_values[2]),
        (p3, delta_dris, DRIs[3], "Delta", p_values[3])
    ])
        histogram!(p[1], values,
            alpha=1.0,
            legend=false,
            # title=title,
            # titlefontsize=10,
            margin=5Plots.mm,
            bins=:auto,
            stroke=:none)
        
        vline!(p[1], [observed], color=:red, label="Observed Score")
        mean_random = mean(skipmissing(values))
        vline!(p[1], [mean_random], color=:blue, label="Mean Random Score")
        
        xlabel!(p[1], title, fontsize=10, margin=10Plots.mm)
        ylabel!(p[1], "Count", fontsize=10, margin=10Plots.mm)
        
        # Add Score value, mean, and p-value annotations
        xlims = Plots.xlims(p[1])
        ylims = Plots.ylims(p[1])
        x_ann = xlims[1] + 0.12 * (xlims[2] - xlims[1])
        y_ann = ylims[2] + 0.17 * (ylims[2] - ylims[1])
        annotate!(p[1], x_ann, y_ann, text("Actual = $(round(observed, digits=2))", :red, 9, :left))
        annotate!(p[1], x_ann, y_ann - 0.05 * (ylims[2] - ylims[1]), text("Mean = $(round(mean_random, digits=2))", :blue, 9, :left))
        annotate!(p[1], x_ann, y_ann - 0.1 * (ylims[2] - ylims[1]), text("p = $(round(p_val, digits=3))", :black, 9, :left))
    end

    # Combine all plots
    p = plot(
        title_plot,
        p1, p2, p3,
        layout = @layout([A{0.17h}; [B C D]]),
        size = (1000, 400)
    )

    return p
end

function calculate_ρ(IC) 
    cor(IC.Q,IC.R)
end

function permuted_prefs(data, case, case_name)

    Fs = []
    Prefs = []


    # First get F and Pref matrices for both stages.
    for stage in 1:2
        data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)

        F_df, Pref_df = prepare_data(data_filtered)

        F = Matrix(F_df) .* 1.0
        Pref = Matrix(Pref_df) .* 1.0


        push!(Fs, F)
        push!(Prefs, Pref)
    end

    ### Next do actual DRI calculations. Also save the ICs and a matrix with Spearman's correlations for all pairs.
    ICs = []
    DRIs = []
    ρs = []
    F_cors = []
    Pref_cors = []
    means = []
    vars = []
    for stage in 1:2
        IC = calculate_IC(Fs[stage], Prefs[stage])
        push!(ICs, IC)

        DRI = calculate_dri_variant(IC)
        push!(DRIs, DRI)

        ρ = calculate_ρ(IC)
        push!(ρs, ρ)

        F_cor = pairwise_correlations(Fs[stage]) 
        Pref_cor = pairwise_correlations(Prefs[stage]) 
        push!(F_cors, F_cor)
        push!(Pref_cors, Pref_cor)

        push!(means, (mean(IC.Q), mean(IC.R)))
        push!(vars, (var(IC.Q), var(IC.R)))

    end

    # Add the deltas
    push!(DRIs, DRIs[2] - DRIs[1])
    push!(ρs, ρs[2] - ρs[1])

    nSamples = 10000
    Random.seed!(123)

    DRIs_permuted = []
    ICs_permuted_sample = []
    for i in 1:nSamples
        print("\rCase $case ($case_name): $i/$nSamples")
        DRIs_permuted_sample = []
        ICs_permuted = []
        n_pairs = size(ICs[1])[1]
        n_users = size(Fs[1])[2]

        # shuffled_indices = rand(1:n_users, n_users)
        # shuffled_indices_pairs = rand(1:n_pairs, n_pairs)
        # shuffle
        shuffled_indices = shuffle(1:n_users)
        shuffled_indices_pairs = shuffle(1:n_pairs)


        for stage in 1:2

            F = Fs[stage]

            IC_permuted = copy(ICs[stage])

            if shuffle_method == :individuals
                Pref_cor = Pref_cors[stage]
                IC_permuted.R = [ Pref_cor[shuffled_indices[IC_permuted.P1[k]], shuffled_indices[IC_permuted.P2[k]]] for k in 1:n_pairs]
                # IC_permuted.R = [ Pref_cor[IC_permuted.P1[k], IC_permuted.P2[k]] for k in 1:n_pairs]
            elseif shuffle_method == :pairs
                IC_permuted.R = IC_permuted.R[shuffled_indices_pairs]
            else 
                error("Unknown shuffle method $shuffle_method")
            end

            DRI_permuted = calculate_dri_variant(IC_permuted)

            push!(ICs_permuted, IC_permuted)
            push!(DRIs_permuted_sample, DRI_permuted)
        end
        if i == 1
            ICs_permuted_sample = ICs_permuted
        end

        # Add the delta DRI
        push!(DRIs_permuted_sample, DRIs_permuted_sample[2] - DRIs_permuted_sample[1])
 
        push!(DRIs_permuted, DRIs_permuted_sample)
    end
    println(". Done")

    p_values = []  
    N = length(DRIs_permuted)

    for stage in 1:3
        n_greater = count(DRI_permuted[stage] >= DRIs[stage] for DRI_permuted in DRIs_permuted)
        p_right = n_greater / N
        p_left  = 1.0 - p_right
        p_value = 2 * min(p_left, p_right)
        p_value = min(p_value, 1.0)
        push!(p_values, p_value)
    end

    permuted_means = [mean([DRIs_permuted[j][i] for j in 1:length(DRIs_permuted)]) for i in 1:3]

    # permuted_vars = [
    #     var([DRIs_permuted[j][i] for j in 1:N]) 
    #     for i in 1:3
    # ]

    # Calculate 99% confidence intervals
    confidence_intervals = [
        (quantile([DRIs_permuted[j][i] for j in 1:N], 0.005), quantile([DRIs_permuted[j][i] for j in 1:N], 0.995))
        for i in 1:3
    ]

    p = DRI_Comparison_Plot(case, case_name, ICs, ICs_permuted_sample, "Actual v. Permuted Preferences", "Permuted Preferences (Single Sample)")

    # Create and save histograms
    hist_p = create_dri_histograms(case, case_name, DRIs, DRIs_permuted, p_values)

    # Combine plots side by side
    # layout=(2,1),
    combined_p = plot(p, hist_p, layout = @layout([A{0.70h}; B]), size=(900,1200), margin=0Plots.mm, left_margin=5Plots.mm, right_margin=5Plots.mm, subplot_margins=0Plots.mm)
    
    return combined_p, DRIs, permuted_means, p_values, ρs, confidence_intervals
end

# 2) define Bayesian shrinkage toward the null‐mean (permuted_mean)
# function bayes_shrink(x_actual::Float64, μ0::Float64, σ2::Float64, τ2::Float64)
    # return (τ2 / (τ2 + σ2)) * x_actual + (σ2 / (τ2 + σ2)) * μ0
# end

"""
    create_dri_comparison_chart(results)

Creates a horizontal bar chart comparing DRI values across different cases, showing:
- Pre-deliberation DRI
- Post-deliberation DRI
- Delta DRI
Each with their corresponding permuted means and ρ values.
"""
function create_dri_comparison_chart(results)
    # Find the range of values
    min_val = minimum([minimum(results.DRIPre), minimum(results.MeanPermutedDRIPre),
                      minimum(results.DRIPost), minimum(results.MeanPermutedDRIPost),
                      minimum(results.DRIDelta), minimum(results.MeanPermutedDRIDelta)])
    max_val = maximum([maximum(results.DRIPre), maximum(results.MeanPermutedDRIPre),
                      maximum(results.DRIPost), maximum(results.MeanPermutedDRIPost),
                      maximum(results.DRIDelta), maximum(results.MeanPermutedDRIDelta)])
    
    # Add some padding to the range
    val_range = max_val - min_val
    padding = 0.1 * val_range
    x_min = min_val - padding
    x_max = max_val + padding

    n = nrow(results)
    
    # Create the plot
    p = plot(
        size=(1500, 500 + 30 * (n + 1)),  # Increased height to accommodate titles and means row
        left_margin=70Plots.mm,
        bottom_margin=10Plots.mm,
        top_margin=0Plots.mm,
        right_margin=10Plots.mm,
        xlims=(x_min, x_max),
        yticks=[],
        layout=@layout([A{0.1h}; B{0.05h}; [C D E]])  # Title row, legend row, and three subplots below
    )
    
    # Add main title
    plot!(p[1],
        title="DRI Analysis: Actual vs Permuted Preferences ($dri_variant DRI)",
        grid=false,
        showaxis=false,
        bottom_margin=-140Plots.px
    )

    # Create legend-only plot
    plot!(p[2],
        grid=false,
        showaxis=false,
        legend=:outertop,
        top_margin=-20Plots.px,
        bottom_margin=-40Plots.px,
        # axis=false,
        # framestyle=:none,
        # xaxis=false,
        # yaxis=false,
        # xticks=false,
        # yticks=false,
        # xlims=(0,0),
        # ylims=(0,0)
    )
    # Add dummy series for legend
    bar!(p[2], [], [], label="Actual DRI", color=:blue, alpha=0.7, linewidth=0, bar_width=0.3, xlims=(0,0), ylims=(0,0))
    scatter!(p[2], [], [], label="Permuted DRI (Null Distribution)", marker=:x, color=:black, markerstrokewidth=0.5, markersize=3, xlims=(0,0), ylims=(0,0))

    y_positions = 1:1.5:(1.5 * n)
    bar_width = 0.3
    spacing = 0.4
    
    ### Pre-deliberation
    begin
        plot!(p[3],
            title="Pre-Deliberation",
            legend=false,
            left_margin=0Plots.mm,
            top_margin=0Plots.mm,
        )

        # DRI
        bar!(p[3],
            y_positions,
            results.DRIPre,
            bar_width=bar_width,
            # label="Actual DRI Pre",
            color=:blue,
            alpha=0.7,
            orientation=:h,
            linewidth=0
        )

        # Add DRI value annotations
        for (i, y_pos) in enumerate(y_positions)
            # Calculate x position - always place to the right of the bar
            x_pos = max(0, results.DRIPre[i]) + 0.02 * val_range
            # Check if DRI value is within confidence interval
            p_color = i == n || (results.DRIPre[i] >= results.CIPre[i][1] && results.DRIPre[i] <= results.CIPre[i][2]) ? :black : :darkgreen
            # Only show p-value for case rows (not means row)
            annotation_text = i == n ? "$(round(results.DRIPre[i], digits=2))" : "$(round(results.DRIPre[i], digits=2)) (p=$(round(results.pValuePre[i], digits=3)))"
            annotate!(p[3],
                x_pos,
                y_pos,
                text(annotation_text, 9, p_color, :left)
            )
        end

        # Error Bars
        scatter!(p[3],
            results.MeanPermutedDRIPre[1:n-1], 
            y_positions[1:n-1] .+ spacing,
            marker=:x, color=:black, markerstrokewidth=0.5, markersize=3,
            xerror=[(results.MeanPermutedDRIPre[i] - results.CIPre[i][1], results.CIPre[i][2]  - results.MeanPermutedDRIPre[i]) for i in 1:n-1]
        )

    end
    
    ### Post-deliberation
    begin 
        plot!(p[4],
            title="Post-Deliberation",
            legend=false,
            left_margin=0Plots.mm,
            top_margin=0Plots.mm,
        )

        # DRI
        bar!(p[4],
            y_positions,
            results.DRIPost,
            bar_width=bar_width,
            # label="Actual DRI Post",
            color=:blue,
            alpha=0.7,
            orientation=:h,
            linewidth=0
        )

        # Add DRI value annotations
        for (i, y_pos) in enumerate(y_positions)
            # Calculate x position - always place to the right of the bar
            x_pos = max(0, results.DRIPost[i]) + 0.02 * val_range
            # Check if DRI value is within confidence interval
            p_color = i == n || (results.DRIPost[i] >= results.CIPost[i][1] && results.DRIPost[i] <= results.CIPost[i][2]) ? :black : :darkgreen
            # Only show p-value for case rows (not means row)
            annotation_text = i == n ? "$(round(results.DRIPost[i], digits=2))" : "$(round(results.DRIPost[i], digits=2)) (p=$(round(results.pValuePost[i], digits=3)))"
            annotate!(p[4],
                x_pos,
                y_pos,
                text(annotation_text, 9, p_color, :left)
            )
        end

        # Error Bars
        scatter!(p[4],
            results.MeanPermutedDRIPost[1:n-1], 
            y_positions[1:n-1] .+ spacing,
            marker=:x, color=:black, markerstrokewidth=0.5, markersize=3,
            xerror=[(results.MeanPermutedDRIPost[i] - results.CIPost[i][1], results.CIPost[i][2]  - results.MeanPermutedDRIPost[i]) for i in 1:n-1]
        )

    end 

    ### Delta
    begin
        plot!(p[5],
            title="Delta",
            legend=false,  # Remove legend from Delta plot since it's now in the dedicated legend plot
            left_margin=0Plots.mm,
            top_margin=0Plots.mm,
        )

        # DRI
        bar!(p[5],
            y_positions,
            results.DRIDelta,
            bar_width=bar_width,
            label="Actual DRI",
            color=:blue,
            alpha=0.7,
            orientation=:h,
            linewidth=0
        )

        # Add DRI value annotations
        for (i, y_pos) in enumerate(y_positions)
            # Calculate x position - always place to the right of the bar
            x_pos = max(0, results.DRIDelta[i]) + 0.02 * val_range
            # Check if DRI value is within confidence interval
            p_color = i == n || (results.DRIDelta[i] >= results.CIDelta[i][1] && results.DRIDelta[i] <= results.CIDelta[i][2]) ? :black : :darkgreen
            # Only show p-value for case rows (not means row)
            annotation_text = i == n ? "$(round(results.DRIDelta[i], digits=2))" : "$(round(results.DRIDelta[i], digits=2)) (p=$(round(results.pValueDelta[i], digits=3)))"
            annotate!(p[5],
                x_pos,
                y_pos,
                text(annotation_text, 9, p_color, :left)
            )
        end

        # Error Bars
        scatter!(p[5],
            results.MeanPermutedDRIDelta[1:n-1], 
            y_positions[1:n-1] .+ spacing,
            label="Permuted DRI (Null Distribution)",
            marker=:x, color=:black, markerstrokewidth=0.5, markersize=3,
            xerror=[(results.MeanPermutedDRIDelta[i] - results.CIDelta[i][1], results.CIDelta[i][2]  - results.MeanPermutedDRIDelta[i]) for i in 1:n-1]
        )

    end

    # Add case names and p-values as text annotations (only on the first subplot)
    for (i, case_name) in enumerate(results.CaseName)
        y_pos = y_positions[i]
        if i == n
            # Means row
            annotate!(p[3],
                x_min - 0.01 * val_range,
                y_pos,
                text("Means", 8, :right)
            )
        else
            # Add p-value annotation
            annotate!(p[3],
                x_min - 0.01 * val_range,
                y_pos,
                text("Case $(results.case[i]): $case_name", 8, :right)
            )
        end
    end
    
    # Add vertical lines at x=0 for all subplots
    for i in 3:5
        vline!(p[i], [0], color=:black, label="")
    end
    
    # Axis labels
    # xlabel!(p[2], "DRI Value")
    # xlabel!(p[3], "DRI Value")
    # xlabel!(p[4], "Delta DRI")
    
    return p
end

function permuted_prefs_main()
    # Read input data
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    
    # Get unique cases
    cases = unique(data.CaseID)

    results = DataFrame(
        case = Float64[],
        CaseName = String[],
        ρPre = Float64[],
        ρPost = Float64[],
        ρDelta = Float64[],
        DRIPre = Float64[],
        DRIPost = Float64[],
        DRIDelta = Float64[],
        MeanPermutedDRIPre = Float64[],
        MeanPermutedDRIPost = Float64[],
        MeanPermutedDRIDelta = Float64[],
        pValuePre = Float64[],
        pValuePost = Float64[],
        pValueDelta = Float64[],
        CIPre = Tuple{Float64,Float64}[],
        CIPost = Tuple{Float64,Float64}[],
        CIDelta = Tuple{Float64,Float64}[],
    )

    # Process each case
    for case in cases
        case_data = filter(r -> r.CaseID == case, data)
        case_name = case_data[1, :Case]
            
        print("Case $case ($case_name)")
        (p, DRIs, permuted_means, p_values, ρs, confidence_intervals) = permuted_prefs(case_data, case, case_name)

        # Add to results
        push!(results, (
            case,
            case_name,
            ρs[1],
            ρs[2],
            ρs[3],
            DRIs...,
            permuted_means...,
            p_values...,
            confidence_intervals[1],
            confidence_intervals[2],
            confidence_intervals[3],
        ))

        outdir = "local-output/permuted-preferences/"
        mkpath(outdir)
        savefig(p, "$outdir/case-$case-permuted-preferences-$dri_variant.png")
    end

    # Calculate means and add them as a new row
    mean_row = (
        -1.0,  # Special case number for means
        "Means",
        mean(results.ρPre),
        mean(results.ρPost),
        mean(results.ρDelta),
        mean(results.DRIPre),
        mean(results.DRIPost),
        mean(results.DRIDelta),
        mean(results.MeanPermutedDRIPre),
        mean(results.MeanPermutedDRIPost),
        mean(results.MeanPermutedDRIDelta),
        mean(results.pValuePre),
        mean(results.pValuePost),
        mean(results.pValueDelta),
        (mean([ci[1] for ci in results.CIPre]), mean([ci[2] for ci in results.CIPre])),
        (mean([ci[1] for ci in results.CIPost]), mean([ci[2] for ci in results.CIPost])),
        (mean([ci[1] for ci in results.CIDelta]), mean([ci[2] for ci in results.CIDelta])),
    )
    push!(results, mean_row)

    # Save results to CSV
    CSV.write("local-output/permuted-preferences/permuted-preferences-$dri_variant-results.csv", results)
    
    # Create horizontal bar chart of DRI values
    sort!(results, [:case], by=x -> x, rev=true)
    
    # Create and save the comparison chart
    p = create_dri_comparison_chart(results)
    savefig(p, "local-output/permuted-preferences/dri-comparison-permuted-preferences-$dri_variant.png")
    
    return results
end
permuted_prefs_main()





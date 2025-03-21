using CSV
using DataFrames
using Statistics
using Plots
using GLM

# Import functions from resampling_dri.jl
include("resampling_dri.jl")

"""
    calculate_iai(F_df::DataFrame, Pref_df::DataFrame)

Calculates the Intersubjective Agreement Index (IAI) as the average pairwise Spearman correlation
across all pairs of users, treating considerations and preferences as a single set of responses.
"""
function calculate_iai(F_df::DataFrame, Pref_df::DataFrame)
    # Concatenate F_df and Pref_df vertically
    combined_df = vcat(F_df, Pref_df)
    
    # Calculate pairwise Spearman correlations using reference implementation
    corr_df = pairwise_spearman(combined_df)
    
    # Calculate mean correlation (excluding self-correlations)
    iai = mean(skipmissing(corr_df.Freq))
    
    return iai
end


function analyze_case_stage(data::DataFrame, case_id::String, stage_id::Int)
    # Filter data for specific case and stage
    data_filtered = filter(r -> string(r.CaseID) == case_id && r.StageID == stage_id, data)

    # Calculate DRI
    F_df, Pref_df = prepare_data(data_filtered)
    IC = calculate_IC(F_df, Pref_df)
    dri = calculate_dri(IC)
    
    # Calculate IAI
    iai = calculate_iai(F_df, Pref_df)
    
    return dri, iai
end

function analyze_all_stages(data::DataFrame)
    # Get unique cases and stages
    cases = unique(data.CaseID)
    stages = unique(data.StageID)
    
    # Initialize results DataFrame
    results = DataFrame(
        CaseID = String[],
        StageID = Int[],
        DRI = Float64[],
        IAI = Float64[]
    )
    
    # Analyze each case and stage
    for case in cases
        for stage in stages
            dri, iai = analyze_case_stage(data, string(case), stage)
            push!(results, (string(case), stage, dri, iai))
        end
    end
    
    return results
end

function calculate_changes(results::DataFrame)
    # Sort by CaseID and StageID
    sort!(results, [:CaseID, :StageID])
    
    # Calculate changes
    changes = DataFrame(
        CaseID = String[],
        StageID = Int[],
        DRI_Change = Float64[],
        IAI_Change = Float64[],
        DRI_Change_Pct = Float64[],
        IAI_Change_Pct = Float64[]
    )
    
    # Calculate changes between consecutive stages for each case
    for case in unique(results.CaseID)
        case_data = filter(r -> r.CaseID == case, results)
        for i in 2:nrow(case_data)
            # Skip if starting DRI is negative
            if case_data[i-1, :DRI] < 0
                continue
            end
            
            # Calculate absolute changes
            dri_change = case_data[i, :DRI] - case_data[i-1, :DRI]
            iai_change = case_data[i, :IAI] - case_data[i-1, :IAI]
            
            # Calculate percentage changes
            dri_change_pct = (dri_change / abs(case_data[i-1, :DRI])) * 100
            iai_change_pct = (iai_change / abs(case_data[i-1, :IAI])) * 100
            
            push!(changes, (
                case,
                case_data[i, :StageID],
                dri_change,
                iai_change,
                dri_change_pct,
                iai_change_pct
            ))
        end
    end
    
    return changes
end

function plot_changes(changes::DataFrame, results::DataFrame)
    # Create output directory if it doesn't exist
    output_dir = "Output/iai_analysis"
    mkpath(output_dir)
    
    # Plot 1: Absolute Changes in DRI vs IAI
    p1 = scatter(
        changes.DRI_Change,
        changes.IAI_Change,
        xlabel = "Change in DRI",
        ylabel = "Change in IAI",
        title = "Relationship between Absolute Changes in DRI and IAI",
        legend = false
    )
    
    # Add trend line for absolute changes
    model1 = lm(@formula(IAI_Change ~ DRI_Change), changes)
    x_range1 = range(minimum(changes.DRI_Change), maximum(changes.DRI_Change), length=100)
    y_pred1 = predict(model1, DataFrame(DRI_Change = x_range1))
    plot!(p1, x_range1, y_pred1, color=:red, linewidth=2)
    
    # Add correlation coefficient for absolute changes
    corr1 = cor(changes.DRI_Change, changes.IAI_Change)
    annotate!(p1, 0.02, 0.98, text("r = $(round(corr1, digits=3))", :left, :top))
    
    # Plot 2: Percentage Changes in DRI vs IAI
    p2 = scatter(
        changes.DRI_Change_Pct,
        changes.IAI_Change_Pct,
        xlabel = "Percentage Change in DRI",
        ylabel = "Percentage Change in IAI",
        title = "Relationship between Percentage Changes in DRI and IAI",
        legend = false
    )
    
    # Add trend line for percentage changes
    model2 = lm(@formula(IAI_Change_Pct ~ DRI_Change_Pct), changes)
    x_range2 = range(minimum(changes.DRI_Change_Pct), maximum(changes.DRI_Change_Pct), length=100)
    y_pred2 = predict(model2, DataFrame(DRI_Change_Pct = x_range2))
    plot!(p2, x_range2, y_pred2, color=:red, linewidth=2)
    
    # Add correlation coefficient for percentage changes
    corr2 = cor(changes.DRI_Change_Pct, changes.IAI_Change_Pct)
    annotate!(p2, 0.02, 0.98, text("r = $(round(corr2, digits=3))", :left, :top))
    
    # Plot 3: DRI vs IAI (direct comparison)
    p3 = scatter(
        results.DRI,
        results.IAI,
        xlabel = "DRI",
        ylabel = "IAI",
        title = "Relationship between DRI and IAI",
        legend = false
    )
    
    # Add trend line for direct comparison
    model3 = lm(@formula(IAI ~ DRI), results)
    x_range3 = range(minimum(results.DRI), maximum(results.DRI), length=100)
    y_pred3 = predict(model3, DataFrame(DRI = x_range3))
    plot!(p3, x_range3, y_pred3, color=:red, linewidth=2)
    
    # Add correlation coefficient for direct comparison
    corr3 = cor(results.DRI, results.IAI)
    annotate!(p3, 0.02, 0.98, text("r = $(round(corr3, digits=3))", :left, :top))
    
    # Save all plots
    savefig(p1, "$(output_dir)/dri_iai_changes_abs.png")
    savefig(p2, "$(output_dir)/dri_iai_changes_pct.png")
    savefig(p3, "$(output_dir)/dri_iai_direct.png")
    
    return p1, p2, p3
end

function main()
    # Create output directory if it doesn't exist
    output_dir = "Output/iai_analysis"
    mkpath(output_dir)
    
    # Read data
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    
    # Analyze all stages
    results = analyze_all_stages(data)
    
    # Calculate changes
    changes = calculate_changes(results)
    
    # Create plots
    p1, p2, p3 = plot_changes(changes, results)
    
    # Save results
    CSV.write("$(output_dir)/results.csv", results)
    CSV.write("$(output_dir)/changes.csv", changes)
    
    # Print summary statistics
    println("\nSummary Statistics:")
    println("Number of cases analyzed: ", length(unique(changes.CaseID)))
    println("\nAbsolute Changes:")
    println("Correlation between DRI and IAI changes: ", cor(changes.DRI_Change, changes.IAI_Change))
    println("Mean DRI change: ", mean(changes.DRI_Change))
    println("Mean IAI change: ", mean(changes.IAI_Change))
    println("\nPercentage Changes:")
    println("Correlation between DRI and IAI percentage changes: ", cor(changes.DRI_Change_Pct, changes.IAI_Change_Pct))
    println("Mean DRI percentage change: ", mean(changes.DRI_Change_Pct), "%")
    println("Mean IAI percentage change: ", mean(changes.IAI_Change_Pct), "%")
    println("\nDirect Comparison:")
    println("Correlation between DRI and IAI: ", cor(results.DRI, results.IAI))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end 
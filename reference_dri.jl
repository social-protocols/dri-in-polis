using CSV
using DataFrames
using Statistics
using StatsBase
using HypothesisTests
using Printf
using Dates
using LinearAlgebra
using KernelDensity
using Colors
using Plots
using Random

function dri_calc(df::DataFrame, v1::Symbol, v2::Symbol)
    λ = 1 - sqrt(2)/2
    # Compute the absolute differences divided by sqrt(2) and take the mean (ignoring missing values)
    diffs = abs.(df[!, v1] .- df[!, v2]) ./ sqrt(2)
    mean_diff = mean(skipmissing(diffs))
    # Apply the transformation (scaling the average difference to the -1 to 1 scale)
    dri = 2 * (((1 - mean_diff) - λ) / (1 - λ)) - 1
    return dri
end

function remove_all_missing!(df::DataFrame)
    for name in names(df)
        if all(ismissing, df[!, name])
            select!(df, Not(name))
        end
    end
    return df
end

function pairwise_spearman(df::DataFrame)
    cols = names(df)
    result = DataFrame(Var1 = String[], Var2 = String[], Freq = Union{Missing, Float64}[])
    for i in 1:length(cols)
        for j in 1:(i-1)
            x = df[!, cols[i]]
            y = df[!, cols[j]]
            # Create a vector of indices where neither x nor y is missing.
            valid_inds = [k for k in eachindex(x) if !ismissing(x[k]) && !ismissing(y[k])]
            if length(valid_inds) <= 1
                corr_val = missing
            else
                x_valid = [x[k] for k in valid_inds]
                y_valid = [y[k] for k in valid_inds]
                # Compute the average (tied) ranks for each valid vector.
                rx = float.(tiedrank(x_valid))
                ry = float.(tiedrank(y_valid))
                if length(rx) == 1 || length(ry) == 1
                    corr_val = missing
                else
                    corr_val = cor(rx, ry)
                end
            end
            push!(result, (cols[i], cols[j], corr_val))
        end
    end
    return result
end

function set_names!(df::DataFrame, newnames::Vector{String})
    rename!(df, Dict(old => new for (old, new) in zip(names(df), newnames)))
end


# --- Create output folder if it does not exist ---

outdir = "Output/reference"
mkpath(outdir)


"""
    dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64)

Creates a single DRI plot in the given subplot with:
  - Scatter points with jitter
  - Reference lines (diagonal, horizontal at 0, vertical at 0)
  - DRI annotation
  - Custom axis labels and title
"""
function dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64)
    # Add jitter (similar to geom_jitter with width/height = 0.02)
    jitter_width = 0.02
    x_jitter = xdata .+ jitter_width .* randn(length(xdata))
    y_jitter = ydata .+ jitter_width .* randn(length(ydata))

    # Create the scatter plot
    scatter!(p, x_jitter, y_jitter,
        markersize = 2,
        markercolor = RGBA(0.1, 0, 1, 0.6),
        label = "",
        xlims = (-1.1, 1.1), ylims = (-1.1, 1.1),
        legend = false,
        title = title_str,
        titlefontsize = 10,
        grid = false)
    
    # Add reference lines
    plot!(p, [-1.1, 1.1], [-1.1, 1.1], color = :black, lw = 1, label = "")
    hline!(p, [0], color = :black, lw = 1, label = "")
    vline!(p, [0], color = :black, lw = 1, label = "")
    
    # Add DRI annotation (moved further right)
    x_ann = -1.1 + 0.18 * (1.1 - (-1.1))
    y_ann = -1.1 + 0.9 * (1.1 - (-1.1))
    annotate!(p, x_ann, y_ann, text("DRI = $(round(dri, digits=2))", :red, 12))
    
    # Add axis labels
    xlabel!(p, "Intersubjective Agreement - Considerations")
    ylabel!(p, "Intersubjective Agreement - Preferences")
end

"""
    dri_plot_side_by_side(xdata_pre::Vector, ydata_pre::Vector, xdata_post::Vector, ydata_post::Vector, 
                         main_title::String, dri_pre::Float64, dri_post::Float64)

Creates a side-by-side plot with pre and post deliberation data, including:
  - Main title for the overall figure
  - Subplot titles for pre and post deliberation
  - Same plotting features as dri_plot for each subplot
"""
function dri_plot_side_by_side(xdata_pre::Vector, ydata_pre::Vector, xdata_post::Vector, ydata_post::Vector, 
                             main_title::String, dri_pre::Float64, dri_post::Float64)
    # Create a layout with a title plot and two subplots
    title_plot = plot(title = main_title, grid = false, showaxis = false, bottom_margin = -50Plots.px)
    
    # Create the subplots
    p1 = plot(layout=(1,1), size=(400,400), dpi=100, margin=5Plots.mm)
    p2 = plot(layout=(1,1), size=(400,400), dpi=100, margin=5Plots.mm)
    
    # Create pre-deliberation subplot
    dri_plot(p1[1], xdata_pre, ydata_pre, "Pre-Deliberation", dri_pre)
    
    # Create post-deliberation subplot
    dri_plot(p2[1], xdata_post, ydata_post, "Post-Deliberation", dri_post)
    
    # Combine all plots with the layout, using more space for the title section
    p = plot(title_plot, p1, p2, layout = @layout([A{0.1h}; [B C]]), size=(800,400))
    
    return p
end

function main() 

    # Ensure reproducibility for jittering
    Random.seed!(123)

    # Create output directories if they don't exist
    mkpath("$outdir/Figures")

    # --- Read input data ---
    data_orig = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)

    StudyU = unique(data_orig[:, [:StudyID, :CaseID]])
    sort!(StudyU, :StudyID)

    # Initialize global dataframes to accumulate results.
    IC_Global = DataFrame()
    DRIInd_Global = DataFrame()
    DRI_Global = DataFrame()

    # Loop over each study-case
    for row in eachrow(StudyU)

        study = row.StudyID
        case = row.CaseID
        data_study_case = filter(r -> r.StudyID == study && r.CaseID == case, data_orig)
        
        # Get the maximum stage from the whole dataset
        max_stage = maximum(data_orig.StageID)
        IC = nothing  # will initialize on stage 1
        
        # Loop through each analysis stage
        for stage in 1:max_stage
            # stage = 1
            data_analyse = filter(r -> r.StageID == stage, data_study_case)
            if nrow(data_analyse) > 0
                # Get participant numbers from the analysis data
                PNums = data_analyse.PNum

                # Select columns for Q: F1 to F50.
                # (Assumes column names are like "F1", "F2", … "F50")
                Q_cols = filter(c -> occursin(r"^F\d+$", c) && parse(Int, replace(c, "F" => "")) ≤ 50, names(data_analyse))
                Q = data_analyse[:, Q_cols]
                remove_all_missing!(Q)

                # Select columns for R: Pref1 to Pref10.
                R_cols = filter(c -> occursin(r"^Pref\d+$", c) && parse(Int, replace(c, "Pref" => "")) ≤ 10, names(data_analyse))
                R = data_analyse[:, R_cols]
                remove_all_missing!(R)
                
                # Transpose Q and R.
                # Convert to a matrix, transpose, then convert back to DataFrame.
                Q_mat = Matrix(Q)
                R_mat = Matrix(R)
                Q_trans = permutedims(Q_mat)
                R_trans = permutedims(R_mat)
                Q_df = DataFrame(Q_trans, :auto)
                R_df = DataFrame(R_trans, :auto)
                # Set the column names to the participant numbers (converted to strings)
                set_names!(Q_df, string.(PNums))
                set_names!(R_df, string.(PNums))
                
                # Compute pairwise Spearman correlations for Q and R.
                QWrite = pairwise_spearman(Q_df)
                RWrite = pairwise_spearman(R_df)
                
                # For stage 1, initialize IC with a column of paired participant identifiers.
                if stage == 1
                    P_P = [row.Var1 * "-" * row.Var2 for row in eachrow(QWrite)]
                    # Convert participant identifiers to numbers (if possible); otherwise keep as NaN.
                    to_number(s::String) = try parse(Float64, s) catch; NaN; end
                    P1 = [to_number(row.Var1) for row in eachrow(QWrite)]
                    P2 = [to_number(row.Var2) for row in eachrow(QWrite)]
                    IC = DataFrame(P_P = P_P, P1 = P1, P2 = P2)
                end
                
                # Prepare the current stage's Q and R correlations as new columns.
                Q_colname = Symbol("Q" * string(stage))
                R_colname = Symbol("R" * string(stage))
                Q_stage = DataFrame(Q_colname => QWrite.Freq)
                R_stage = DataFrame(R_colname => RWrite.Freq)
                
                # Append these columns to the existing IC dataframe.
                IC = hcat(IC, Q_stage, R_stage)
            else
                # In case no data is present for this stage, you might choose to skip or add default zeros.
                # (The original R code creates a matrix of zeros for robustness.)
                continue
            end
        end
        
        # Add study-case information to IC.
        nrows_IC = nrow(IC)
        IC[!, :CaseId] = fill(case, nrows_IC)
        IC[!, :StudyID] = fill(study, nrows_IC)
        # Assume the study name and case name are stored in columns "Study" and "Case"
        study_val = data_study_case[1, :Study]
        case_val = data_study_case[1, :Case]
        IC[!, :Study] = fill(study_val, nrows_IC)
        IC[!, :Case] = fill(case_val, nrows_IC)
        
        # Add data status.
        # Initially set all to 1. Then, for any pair where either participant appears in a row with Datacheck == 3,
        # set status to 3.
        IC[!, :Data_Status] = ones(Int, nrows_IC)
        # data_for_status = data_study_case.PNum[data_study_case.Datacheck .== 3]
        # for i in 1:nrows_IC
        #     if (IC[i, :P1] in data_for_status) || (IC[i, :P2] in data_for_status)
        #         IC[i, :Data_Status] = 3
        #     end
        # end

        # --- Analysis of IC points ---
        IC[!, :IC_PRE] = 1 .- abs.((IC[!, :R1] .- IC[!, :Q1]) ./ sqrt(2))
        IC[!, :IC_POST] = 1 .- abs.((IC[!, :R2] .- IC[!, :Q2]) ./ sqrt(2))
        
        # Compute group DRI levels.
        DRI_PRE = dri_calc(IC, :R1, :Q1)
        DRI_POST = dri_calc(IC, :R2, :Q2)
        
        # Create the case-level DRI dataframe.
        DRI_Case = DataFrame(StudyID = study, Study = study_val, CaseID = case, Case = case_val,
                             DRI_PRE = DRI_PRE, DRI_POST = DRI_POST)
        
        # --- Wilcoxon signed-rank tests ---
        # Compare the paired IC_PRE and IC_POST values.
        # wt_twoside = ApproximateWilcoxonSignedRankTest(IC[!, :IC_POST], IC[!, :IC_PRE]; tail = :two_sided)
        # wt_greater = ApproximateWilcoxonSignedRankTest(IC[!, :IC_POST], IC[!, :IC_PRE]; tail = :right)
        # DRI_Case[!, :DRI_one_tailed_p] = pvalue(wt_greater)
        # DRI_Case[!, :DRI_twoside_p] = pvalue(wt_twoside)
        
        # --- Compute individual DRI values ---
        # Get the unique list of participant identifiers from columns P1 and P2.
        Plist = unique(vcat(IC[!, :P1], IC[!, :P2]))
        sort!(Plist)
        DRIInd = DataFrame(Participant = Plist)
        DRIInd[!, :StudyID] = fill(study, length(Plist))
        DRIInd[!, :Study] = fill(study_val, length(Plist))
        DRIInd[!, :CaseID] = fill(case, length(Plist))
        DRIInd[!, :Case] = fill(case_val, length(Plist))
        DRIInd = select(DRIInd, [:StudyID, :Study, :CaseID, :Case, :Participant])
        
        # For each participant, filter the pairs in IC where they appear and compute pre and post DRI.
        DRIPre_vals = Float64[]
        DRIPost_vals = Float64[]
        for p in Plist
            subset_IC = filter(row -> row.P1 == p || row.P2 == p, IC)
            push!(DRIPre_vals, dri_calc(subset_IC, :R1, :Q1))
            push!(DRIPost_vals, dri_calc(subset_IC, :R2, :Q2))
        end
        DRIInd[!, :DRIPre] = DRIPre_vals
        DRIInd[!, :DRIPost] = DRIPost_vals
        
        # --- Append current case results to global dataframes ---
        if nrow(IC_Global) == 0
            IC_Global = IC
            DRIInd_Global = DRIInd
            DRI_Global = DRI_Case
        else
            append!(IC_Global, IC)
            append!(DRIInd_Global, DRIInd)
            append!(DRI_Global, DRI_Case)
        end

        # --- Generate and save side-by-side DRI plot for this case ---
        case_plot = dri_plot_side_by_side(
            IC[!, :Q1], IC[!, :R1],
            IC[!, :Q2], IC[!, :R2],
            "DRI Plots: $(study_val) - $(case_val) (Case $(case))",
            DRI_PRE, DRI_POST
        )
        savefig(case_plot, "$outdir/Figures/Case_$(case).png")
    end

    # --- Produce additional output for non-control cases ---
    # (Here we filter for StageID==1 and CaseID > 0.9; adjust as needed if CaseID is numeric.)
    DRI0 = filter(r -> r.StageID == 1 && r.CaseID > 0.9, data_orig)
    sort!(DRI0, :StudyID)
    DRIInd_Global_1 = filter(r -> r.CaseID > 0.9, DRIInd_Global)
    # Assuming that DRI0 has a column "PNum" and that the rows align as in the original R code:
    DRIInd_Global_2 = hcat(DRIInd_Global_1, DataFrame(PNum = DRI0.PNum))

    # --- Save output CSV files ---
    CSV.write("$outdir/DRIGlobal.csv", DRI_Global)
    CSV.write("$outdir/DRI_Individual_Global_NO_Control.csv", DRIInd_Global_2)
    CSV.write("$outdir/DRI_Individual_Global.csv", DRIInd_Global)

    # =============================================================================
    # Generate and save Figure 2 plots
    # =============================================================================

    # Create side-by-side plot for Figure 2
    fig2_data = filter(row -> row[:StudyID] == 2, IC_Global)
    fig2_dri = filter(row -> row[:StudyID] == 2, DRI_Global)[1,:]
    
    fig2 = dri_plot_side_by_side(
        fig2_data[!, :Q1], fig2_data[!, :R1],
        fig2_data[!, :Q2], fig2_data[!, :R2],
        "Figure 2. DRI Plots: FNQCJ Case (Case $(fig2_dri.CaseID))",
        fig2_dri.DRI_PRE, fig2_dri.DRI_POST
    )
    savefig(fig2, "$outdir/Figures/Fig2.png")

    # =============================================================================
    # Generate and save Figure 3 plots
    # =============================================================================

    # Figure 3 – Part 1: Control Group
    fig3_control_data = filter(row -> row[:StudyID] == 1 && row[:CaseId] == 0.1, IC_Global)
    fig3_control_dri = filter(row -> row[:StudyID] == 1 && row[:CaseID] == 0.1, DRI_Global)[1,:]
    
    fig3_control = dri_plot_side_by_side(
        fig3_control_data[!, :Q1], fig3_control_data[!, :R1],
        fig3_control_data[!, :Q2], fig3_control_data[!, :R2],
        "Figure 3. DRI Plots: Uppsala Speaks (Control) (Case $(fig3_control_dri.CaseID))",
        fig3_control_dri.DRI_PRE, fig3_control_dri.DRI_POST
    )
    savefig(fig3_control, "$outdir/Figures/Fig3_1Control.png")

    # Figure 3 – Part 2: Group Briefing
    fig3_brief_data = filter(row -> row[:StudyID] == 1 && row[:CaseId] == 1, IC_Global)
    fig3_brief_dri = filter(row -> row[:StudyID] == 1 && row[:CaseID] == 1, DRI_Global)[1,:]
    
    fig3_brief = dri_plot_side_by_side(
        fig3_brief_data[!, :Q1], fig3_brief_data[!, :R1],
        fig3_brief_data[!, :Q2], fig3_brief_data[!, :R2],
        "Figure 3. DRI Plots: Uppsala Speaks Study (Group Briefing) (Case $(fig3_brief_dri.CaseID))",
        fig3_brief_dri.DRI_PRE, fig3_brief_dri.DRI_POST
    )
    savefig(fig3_brief, "$outdir/Figures/Fig3_2Brief.png")

    # Figure 3 – Part 3: Group Building Plus
    fig3_building_data = filter(row -> row[:StudyID] == 1 && row[:CaseId] == 2, IC_Global)
    fig3_building_dri = filter(row -> row[:StudyID] == 1 && row[:CaseID] == 2, DRI_Global)[1,:]
    
    fig3_building = dri_plot_side_by_side(
        fig3_building_data[!, :Q1], fig3_building_data[!, :R1],
        fig3_building_data[!, :Q2], fig3_building_data[!, :R2],
        "Figure 3. DRI Plots: Uppsala Speaks Study (Group Building Plus) (Case $(fig3_building_dri.CaseID))",
        fig3_building_dri.DRI_PRE, fig3_building_dri.DRI_POST
    )
    savefig(fig3_building, "$outdir/Figures/Fig3_3Building.png")

end

main()
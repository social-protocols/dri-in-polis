using CSV
using DataFrames
using Statistics
using StatsBase
using Random
using Plots
using Distributions

"""
    prepare_data(data::DataFrame)

Prepares the input data by:
1. Selecting and processing consideration (F) and preference (Pref) columns
2. Transposing the data to get participant-wise matrices
Returns: Tuple of (F_df, Pref_df) DataFrames
"""
function prepare_data(data::DataFrame)
    # Get participant numbers
    PNums = data.PNum

    # Select columns for considerations (F) and preferences (Pref)
    F_cols = filter(c -> occursin(r"^F\d+$", c) && parse(Int, replace(c, "F" => "")) ≤ 50, names(data))
    Pref_cols = filter(c -> occursin(r"^Pref\d+$", c) && parse(Int, replace(c, "Pref" => "")) ≤ 10, names(data))
    
    F = data[:, F_cols]
    Pref = data[:, Pref_cols]
    
    # Remove columns with all missing values
    remove_all_missing!(F)
    remove_all_missing!(Pref)
    
    # Transpose matrices
    F_mat = Matrix(F)
    Pref_mat = Matrix(Pref)
    F_trans = permutedims(F_mat)
    Pref_trans = permutedims(Pref_mat)
    
    # Convert to DataFrames with participant numbers as column names
    F_df = DataFrame(F_trans, :auto)
    Pref_df = DataFrame(Pref_trans, :auto)
    set_names!(F_df, string.(PNums))
    set_names!(Pref_df, string.(PNums))
    
    return F_df, Pref_df
end

"""
    randomize_tags(F_df::DataFrame, Pref_df::DataFrame)

Randomly reassigns rows between F_df and Pref_df while maintaining the original number of rows in each.
Returns: Tuple of (F_df, Pref_df) DataFrames with randomly assigned rows
"""
function randomize_tags(F_df::DataFrame, Pref_df::DataFrame)
    # Combine all rows
    all_rows = vcat(F_df, Pref_df)
    
    # Randomly assign rows back to F and Pref
    n_rows = nrow(all_rows)
    n_F_rows = nrow(F_df)
    random_indices = randperm(n_rows)
    F_indices = random_indices[1:n_F_rows]
    Pref_indices = random_indices[n_F_rows+1:end]
    
    # Create new DataFrames with randomly assigned rows
    F_random = all_rows[F_indices, :]
    Pref_random = all_rows[Pref_indices, :]
    
    return F_random, Pref_random
end

"""
    calculate_IC(F_df::DataFrame, Pref_df::DataFrame)

Calculates the IC DataFrame from F_df and Pref_df by computing pairwise correlations.
Returns: IC DataFrame with Q and R correlations
"""
function calculate_IC(F_df::DataFrame, Pref_df::DataFrame)
    # Compute pairwise correlations
    F_corr = pairwise_spearman(F_df)
    Pref_corr = pairwise_spearman(Pref_df)
    
    # Create IC DataFrame
    P_P = [row.Var1 * "-" * row.Var2 for row in eachrow(F_corr)]
    to_number(s::String) = try parse(Float64, s) catch; NaN; end
    P1 = [to_number(row.Var1) for row in eachrow(F_corr)]
    P2 = [to_number(row.Var2) for row in eachrow(F_corr)]
    
    IC = DataFrame(P_P = P_P, P1 = P1, P2 = P2)
    
    # Add correlation columns directly without stage-specific naming
    IC[!, :Q] = F_corr.Freq
    IC[!, :R] = Pref_corr.Freq
    
    return IC
end

"""
    calculate_dri(IC::DataFrame)

Calculates the DRI using the IC DataFrame.
Returns: DRI value
"""
function calculate_dri(IC::DataFrame)
    return dri_calc(IC, :R, :Q)
end

"""
    calculate_individual_dri(IC::DataFrame)

Calculates individual DRI values for each participant.
Returns: DataFrame with participant IDs and their DRI values
"""
function calculate_individual_dri(IC::DataFrame)
    # Get unique participants
    Plist = unique(vcat(IC[!, :P1], IC[!, :P2]))
    sort!(Plist)
    
    # Create DataFrame for results
    DRIInd = DataFrame(Participant = Plist)
    
    # Calculate DRI for each participant
    dri_vals = Float64[]
    for p in Plist
        subset_IC = filter(row -> row.P1 == p || row.P2 == p, IC)
        push!(dri_vals, calculate_dri(subset_IC))
    end
    DRIInd[!, :DRI] = dri_vals
    
    return DRIInd
end


# Helper functions from original code
function dri_calc(df::DataFrame, v1::Symbol, v2::Symbol)
    λ = 1 - sqrt(2)/2
    diffs = abs.(df[!, v1] .- df[!, v2]) ./ sqrt(2)
    mean_diff = mean(skipmissing(diffs))
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
            valid_inds = [k for k in eachindex(x) if !ismissing(x[k]) && !ismissing(y[k])]
            if length(valid_inds) <= 1
                corr_val = missing
            else
                x_valid = [x[k] for k in valid_inds]
                y_valid = [y[k] for k in valid_inds]
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

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: julia resampling_dri.jl [MODE] [ITERATIONS]")
        println("  MODE: 'original' or 'delta' (default='original')")
        println("  ITERATIONS: number of random iterations (default=100)")
        exit(1)
    end

    # Parse command line arguments
    mode = length(ARGS) > 0 ? ARGS[1] : "original"
    n_iterations = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 100

    if !(mode in ["original", "delta"])
        println("Error: MODE must be either 'original' or 'delta'")
        exit(1)
    end

    # Read input data
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    
    # Get unique cases
    cases = unique(data.CaseID)
    
    # Store results
    if mode == "original"
        results = DataFrame(
            Case = String[],
            CaseName = String[],
            Stage = Int[],
            StageName = String[],
            DRI = Float64[],
            PValue = Float64[]
        )
    else
        results = DataFrame(
            Case = String[],
            CaseName = String[],
            DeltaDRI = Float64[],
            PValue = Float64[]
        )
    end

    # Loop over all cases
    for case in cases
        if mode == "original"
            # Original analysis for each case/stage combination
            for stage in [1, 2]
                data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
                if nrow(data_filtered) == 0
                    println("Warning: No data found for case $case and stage $stage")
                    continue
                end
                
                case_name = data_filtered[1, :Case]
                stage_name = stage == 1 ? "PRE" : "POST"
                
                # Calculate original DRI
                F_df, Pref_df = prepare_data(data_filtered)
                IC_original = calculate_IC(F_df, Pref_df)
                observed_dri = calculate_dri(IC_original)
                
                # Run random DRI calculations
                dri_values = Float64[]
                for i in 1:n_iterations
                    Random.seed!(i)
                    F_random, Pref_random = randomize_tags(F_df, Pref_df)
                    IC_random = calculate_IC(F_random, Pref_random)
                    if IC_random !== nothing
                        dri = calculate_dri(IC_random)
                        push!(dri_values, dri)
                    end
                end
                
                p_value = sum(dri_values .>= observed_dri) / length(dri_values)
                push!(results, (string(case), case_name, Int64(stage), stage_name, observed_dri, p_value))
                
                # Create and save histogram
                p = histogram(dri_values, 
                            title="Distribution of Random DRI Values\nCase $case ($case_name), Stage $stage ($stage_name)",
                            xlabel="DRI",
                            ylabel="Count",
                            legend=false)
                vline!([observed_dri], color=:red, label="Observed DRI")
                
                xlims = Plots.xlims(p)
                ylims = Plots.ylims(p)
                x_ann = xlims[1] + 0.1 * (xlims[2] - xlims[1])
                y_ann = ylims[1] + 0.9 * (ylims[2] - ylims[1])
                annotate!(p, x_ann, y_ann, text("DRI = $(round(observed_dri, digits=2))\np = $(round(p_value, digits=3))", :red, 13))
                
                mkpath("Output/resampling")
                savefig(p, "Output/resampling/dri_distribution_case$(case)_stage$(stage).png")
            end
        else
            # Delta analysis for each case
            data_pre = filter(r -> r.CaseID == case && r.StageID == 1, data)
            data_post = filter(r -> r.CaseID == case && r.StageID == 2, data)
            
            if nrow(data_pre) == 0 || nrow(data_post) == 0
                println("Warning: Missing data for case $case")
                continue
            end
            
            case_name = data_pre[1, :Case]
            
            # Calculate original DRIs
            F_df_pre, Pref_df_pre = prepare_data(data_pre)
            F_df_post, Pref_df_post = prepare_data(data_post)
            IC_pre = calculate_IC(F_df_pre, Pref_df_pre)
            IC_post = calculate_IC(F_df_post, Pref_df_post)
            observed_dri_pre = calculate_dri(IC_pre)
            observed_dri_post = calculate_dri(IC_post)
            observed_delta = observed_dri_post - observed_dri_pre
            
            # Run random delta calculations
            delta_values = Float64[]
            for i in 1:n_iterations
                Random.seed!(i)
                
                # Randomize both stages
                F_random_pre, Pref_random_pre = randomize_tags(F_df_pre, Pref_df_pre)
                F_random_post, Pref_random_post = randomize_tags(F_df_post, Pref_df_post)
                
                IC_random_pre = calculate_IC(F_random_pre, Pref_random_pre)
                IC_random_post = calculate_IC(F_random_post, Pref_random_post)
                
                if IC_random_pre !== nothing && IC_random_post !== nothing
                    dri_pre = calculate_dri(IC_random_pre)
                    dri_post = calculate_dri(IC_random_post)
                    push!(delta_values, dri_post - dri_pre)
                end
            end
            
            # Calculate two-sided p-value (proportion of absolute random values >= absolute observed value)
            p_value = sum(abs.(delta_values) .>= abs(observed_delta)) / length(delta_values)
            
            push!(results, (string(case), case_name, observed_delta, p_value))
            
            # Create and save histogram
            p = histogram(delta_values, 
                        title="Distribution of Random DRI Deltas\nCase $case ($case_name)",
                        xlabel="Delta DRI",
                        ylabel="Count",
                        legend=false)
            vline!([observed_delta], color=:red, label="Observed Delta")
            
            xlims = Plots.xlims(p)
            ylims = Plots.ylims(p)
            x_ann = xlims[1] + 0.1 * (xlims[2] - xlims[1])
            y_ann = ylims[1] + 0.9 * (ylims[2] - ylims[1])
            annotate!(p, x_ann, y_ann, text("Delta = $(round(observed_delta, digits=2))\np = $(round(p_value, digits=3))", :red, 13))
            
            mkpath("Output/resampling")
            savefig(p, "Output/resampling/dri_delta_distribution_case$(case).png")
        end
    end
    
    # Calculate aggregate p-value using Fisher's method
    p_values = results.PValue
    k = length(p_values)
    
    # Handle zero p-values by replacing them with a small number
    p_values = max.(p_values, 1e-10)
    
    # Compute the Fisher's statistic
    T = -2 * sum(log.(p_values))
    
    # Define the chi-squared distribution with 2k degrees of freedom
    chi2_dist = Chisq(2 * k)
    
    # Calculate the combined p-value
    combined_p = 1 - cdf(chi2_dist, T)
    
    println("\nAggregate Results:")
    println("Number of datasets: ", k)
    println("Fisher's statistic: ", T)
    println("Aggregate p-value using Fisher's method: ", combined_p)
    
    # Save results to CSV
    output_file = mode == "original" ? "Output/resampling/dri_results.csv" : "Output/resampling/dri_delta_results.csv"
    CSV.write(output_file, results)
end 
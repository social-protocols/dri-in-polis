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
    randomize_tags(F_df::DataFrame, Pref_df::DataFrame, method::String="method1")

Randomly reassigns rows between F_df and Pref_df while maintaining the original number of rows in each.
Returns: Tuple of (F_df, Pref_df) DataFrames with randomly assigned rows
"""
function randomize_tags(F_df::DataFrame, Pref_df::DataFrame, method::String="method1")
    if method == "method1"
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
    else  # method2
        n_rows = nrow(F_df)
        n_F_rows = n_rows - nrow(Pref_df)
        random_indices = randperm(n_rows)
        F_indices = random_indices[1:n_F_rows]
        Pref_indices = random_indices[n_F_rows+1:end]
        # Create new DataFrames with randomly assigned rows
        F_random = F_df[F_indices, :]
        Pref_random = F_df[Pref_indices, :]
    end
    
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
    global rejected_pairs
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
                    if isnan(corr_val)
                        corr_val = missing
                        rejected_pairs += 1
                    end
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
    if length(ARGS) < 1
        println("Usage: julia resampling_dri.jl [ITERATIONS] [MODE] ")
        println("  ITERATIONS: number of random iterations (default=1000)")
        println("  MODE: 'method1' or 'method2' (default='method1')")
        exit(1)
    end

    # Parse command line arguments
    n_iterations = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1000
    mode = length(ARGS) > 1 ? ARGS[2] : "method1"

    if !(mode in ["method1", "method2"])
        println("Error: MODE must be either 'method1' or 'method2'")
        exit(1)
    end

    # Read input data
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    
    # Get unique cases
    cases = unique(data.CaseID)
    
    # Store results
    results = DataFrame(
        Case = String[],
        CaseName = String[],
        Stage = Int[],
        StageName = String[],
        DRI = Float64[],
        PValue = Float64[]
    )

    # Global counter for rejected pairs
    rejected_pairs = 0

    # Loop over all cases
    for case in cases
        # Create side-by-side plots
        p = plot(layout=(1,2), size=(1200,500), 
                plot_title_fontsize=12,
                margin=20Plots.mm,
                subplot_titles=["Stage 1 (PRE)", "Stage 2 (POST)"])

            
        # Process and plot each stage
        for stage in [1, 2]
            data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
            if nrow(data_filtered) == 0
                println("Warning: No data found for case $case and stage $stage")
                continue
            end

            
            case_name = data_filtered[1, :Case]
            stage_name = stage == 1 ? "PRE" : "POST"
            
            println("\nProcessing Case $case ($case_name) - Stage $stage ($stage_name)")
            
            # Calculate original DRI
            F_df, Pref_df = prepare_data(data_filtered)
            IC_observed = calculate_IC(F_df, Pref_df)
            observed_dri = calculate_dri(IC_observed)

            # Run random DRI calculations
            random_dri_values = Float64[]
            if mode == "method1"
                for i in 1:n_iterations
                    Random.seed!(i)
                    F_random, Pref_random = randomize_tags(F_df, Pref_df, "method1")
                    IC_random = calculate_IC(F_random, Pref_random)
                    random_dri = calculate_dri(IC_random)
                    push!(random_dri_values, random_dri)
                    print("\rIteration $i/$n_iterations")
                end
                println()  # New line after progress
                p_value = sum(random_dri_values .>= observed_dri) / length(random_dri_values)
            else  # method2
                dri_tuples = Tuple{Float64,Float64}[]
                for i in 1:n_iterations
                    Random.seed!(i)
                    F_random, Pref_random = randomize_tags(F_df, Pref_df, "method2")
                    IC_control = calculate_IC(F_random, Pref_df)
                    resampled_dri = calculate_dri(IC_control)
                    IC_random = calculate_IC(F_random, Pref_random)
                    random_dri = calculate_dri(IC_random)
                    push!(dri_tuples, (resampled_dri, random_dri))
                    print("\rIteration $i/$n_iterations")
                end
                println()  # New line after progress
                p_value = count(d -> abs(d[2]) >= abs(d[1]), dri_tuples) / length(dri_tuples)
                mean_resampled_dri = mean([d[1] for d in dri_tuples])

                random_dri_values = [d[2] for d in dri_tuples]
            end

            push!(results, (string(case), case_name, Int64(stage), stage_name, observed_dri, p_value))
            
            # Plot histogram for this stage
            histogram!(p[stage], random_dri_values, 
                alpha=1.0,
                legend=false,
                title="Distribution of Random DRI Values\nCase $case ($case_name) $stage_name",
                titlefontsize=10,
                margin=10Plots.mm,
                bins=:auto,
                stroke=:none)
            
            vline!(p[stage], [observed_dri], color=:red, label="Observed DRI")
            mean_random = mean(skipmissing(random_dri_values))
            vline!(p[stage], [mean_random], color=:blue, label="Mean Random DRI")
            
            xlabel!(p[stage], "DRI", fontsize=10, margin=10Plots.mm)
            ylabel!(p[stage], "Count", fontsize=10, margin=10Plots.mm)
            
            # Add DRI value and p-value annotation
            xlims = Plots.xlims(p[stage])
            ylims = Plots.ylims(p[stage])
            x_ann = xlims[1] + 0.1 * (xlims[2] - xlims[1])
            y_ann = ylims[1] + 0.9 * (ylims[2] - ylims[1])
            annotate!(p[stage], x_ann, y_ann, text("DRI = $(round(observed_dri, digits=2))", :red, 13, :left))
            annotate!(p[stage], x_ann, y_ann - 0.05 * (ylims[2] - ylims[1]), text("Mean = $(round(mean_random, digits=2))", :blue, 13, :left))
            annotate!(p[stage], x_ann, y_ann - 0.1 * (ylims[2] - ylims[1]), text("p = $(round(p_value, digits=3))", :black, 13, :left))
        end
        
        # Save plot
        mkpath("Output/resampling")
        savefig(p, "Output/resampling/dri_distribution_case$(case)_$(mode).png")
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
    println("Aggregate p-value ($mode) using Fisher's method: ", combined_p)
    println("Rejected pairs (zero variance): ", rejected_pairs)
    
    # Save results to CSV
    output_file = "Output/resampling/dri_results_$(mode).csv"
    CSV.write(output_file, results)
end 
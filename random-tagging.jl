using CSV
using DataFrames
using Random
using Plots
using Distributions

include("prepare-data.jl")
include("dri-calculations.jl")

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

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 2
        println("Usage: julia random-tagging.jl [ITERATIONS] [MODE] ")
        println("  ITERATIONS: number of random iterations (default=10000)")
        println("  MODE: 'method1' or 'method2' (default='method1')")
        exit(1)
    end

    # Parse command line arguments
    n_iterations = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 10000
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
    # rejected_pairs = 0

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

                # Calculate two-sided p-value
                N = length(random_dri_values)
                n_greater = count(d -> d >= observed_dri, random_dri_values)
                p_right = n_greater / N
                p_left  = 1.0 - p_right    # fraction that are below observed_dri
                p_value = 2 * min(p_left, p_right)
                p_value = min(p_value, 1.0)

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

                # Calculate two-sided p-value
                random_dri_values = [ d[2] for d in dri_tuples ]
                N = length(random_dri_values)
                n_greater = count(d -> d >= observed_dri, random_dri_values)
                p_right = n_greater / N
                p_left  = 1.0 - p_right
                p_value = 2 * min(p_left, p_right)
                p_value = min(p_value, 1.0)
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
        mkpath("local-output/random-tagging")
        savefig(p, "local-output/random-tagging/dri_distribution_case$(case)_$(mode).png")
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
    # println("Rejected pairs (zero variance): ", rejected_pairs)
    
    # Save results to CSV
    output_file = "local-output/random-tagging/random-tagging_results_$(mode).csv"
    CSV.write(output_file, results)
end 
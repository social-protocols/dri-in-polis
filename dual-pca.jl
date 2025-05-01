using MultivariateStats
using Statistics
using LinearAlgebra
using StatsBase # For rank correlations
using CSV
using DataFrames
using FactorRotations
using Plots
using Distances

const varimax_rotation = true  # Set to true to enable varimax rotation

include("prepare-data.jl")
include("dri-calculations.jl")
include("dri-plot.jl")

"""
    PCAUsers(X, min_variance_explained)

Fit a PCA on the standardized matrix `X` and automatically select
the smallest number of principal components that explain at least
`min_variance_explained` fraction of the variance. Optionally applies
Varimax rotation if `varimax_rotation` is true.

Returns `(X_means, X_stds, W, actual_dimensions, variance_explained)`:
- `X_means`: column means of the original matrix
- `X_stds`: column standard deviations of the original matrix
- `W`: the matrix that projects new data into the reduced dimensions
- `actual_dimensions`: the number of components selected
- `variance_explained`: the fraction of variance explained by the selected components
"""
function PCAUsers(X::Matrix, min_variance_explained::Float64)
    num_items, num_users = size(X)

    # Standardize each user's ranks
    X_means = mean(X, dims=1)
    X_stds = std(X, dims=1, corrected=false)  # population std
    X_standardized_users = (X .- X_means) ./ (X_stds .+ 1e-8)

    # Fit PCA with the maximum possible dimension = min(size(X_standardized_users))
    # We'll then select the number of components that meets the min variance explained
    M_corr_full = fit(PCA, X_standardized_users; maxoutdim=minimum(size(X_standardized_users)), mean=0)

    # Determine how many components are required to reach min_variance_explained
    explained_variances = M_corr_full.prinvars ./ sum(M_corr_full.prinvars)
    cumulative_variances = cumsum(explained_variances)
    # If no dimension reaches min_variance_explained, use all components
    idx = findfirst(≥(min_variance_explained), cumulative_variances)
    actual_dimensions = isnothing(idx) ? length(explained_variances) : idx
    # Ensure at least 2 components are selected
    actual_dimensions = max(2, actual_dimensions)

    # Calculate variance explained by selected components
    variance_explained = sum(explained_variances[1:actual_dimensions])

    # Extract the selected components
    original_loadings_full = projection(M_corr_full)  # ranks x max_possible_dim
    selected_loadings = original_loadings_full[:, 1:actual_dimensions]

    # Compute the PCA scores for the selected components
    Y_corr_full = MultivariateStats.transform(M_corr_full, X_standardized_users)  # users x max_possible_dim
    Y_corr_selected = Y_corr_full[:, 1:actual_dimensions]'

    if varimax_rotation
        # Apply Varimax to the selected loadings
        varimax_result = FactorRotations.rotate(selected_loadings, Varimax())
        rotated_loadings = FactorRotations.loadings(varimax_result)

        # (Optional) print info about top loading items, if desired
        show_top_items = false
        if show_top_items
            for j in 1:actual_dimensions
                component_loadings = rotated_loadings[:, j]
                abs_loadings = abs.(component_loadings)
                sorted_indices = sortperm(abs_loadings, rev=true)
                # Do something with sorted_indices, e.g. print top items
            end
        end

        # Rotate the selected scores as well: Rotated_Scores = T' * Original_Scores
        T = rotation(varimax_result)       # (actual_dimensions x actual_dimensions)
        Y_corr_rotated = T' * Y_corr_selected

        # W is the matrix that, when multiplied by standardized data, gives the reduced representation
        # Original: W_selected = selected_loadings' (size = actual_dimensions x ranks)
        # After rotation: W = T' * selected_loadings'
        W_selected = selected_loadings'
        W = T' * W_selected

        return X_means, X_stds, W, actual_dimensions, variance_explained
    else
        # If no rotation, just return the selected projection matrix
        W = selected_loadings'
        return X_means, X_stds, W, actual_dimensions, variance_explained
    end
end

"""
    reduceDimensions(X, X_means, X_stds, W)

Use the projection matrix `W` (from PCAUsers) to project `X` into
the reduced-dimension space. 
"""
function reduceDimensions(X, X_means, X_stds, W)
    X_standardized_users = (X .- X_means) ./ (X_stds .+ 1e-8)
    return W * X_standardized_users
end

"""
    dualPCA_plot(case, case_name, min_variance_explained, correlation_method)

Performs dual PCA for both F and Pref data, ensuring that at least
`min_variance_explained` fraction of variance is retained.
Calculates DRI before and after deliberation (standard vs. dual PCA),
creates and saves comparison plots, and returns DRI values.
"""
function dualPCA_plot(case, case_name, min_variance_explained::Float64, correlation_method::Symbol)
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    data_filtered = filter(r -> r.CaseID == case && r.StageID == 1, data)
    F_df, Pref_df = prepare_data(data_filtered)

    original_F_dims = size(F_df, 1)
    original_Pref_dims = size(Pref_df, 1)

    F = Matrix(F_df) .* 1.0
    Pref = Matrix(Pref_df) .* 1.0

    # PCA for Stage 1
    F_means, F_stds, WF, F_dimensions, F_variance_explained = PCAUsers(F, min_variance_explained)
    Pref_means, Pref_stds, WPref, Pref_dimensions, Pref_variance_explained = PCAUsers(Pref, 1.0)

    ICs = []
    DRIs = []
    ICs_standard = []
    DRIs_standard = []

    # Loop over 2 stages
    for stage in 1:2
        data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
        F_df, Pref_df = prepare_data(data_filtered)

        # Calculate standard DRI
        IC_standard = calculate_IC(F_df, Pref_df)
        DRI_standard = calculate_dri(IC_standard)
        push!(ICs_standard, IC_standard)
        push!(DRIs_standard, DRI_standard)

        # Calculate dual PCA DRI
        F = Matrix(F_df) .* 1.0
        Pref = Matrix(Pref_df) .* 1.0

        Pref_reduced = reduceDimensions(Pref, Pref_means, Pref_stds, WPref)
        F_reduced = reduceDimensions(F, F_means, F_stds, WF)

        # Optionally center if you want cosine = Pearson
        center = false
        if center
            Pref_reduced = Pref_reduced .- mean(Pref_reduced, dims=1)
            F_reduced = F_reduced .- mean(F_reduced, dims=1)
        end

        IC = calculate_IC(F_reduced, Pref_reduced, correlation_method)
        DRI = calculate_dri(IC)

        push!(ICs, IC)
        push!(DRIs, DRI)
    end

    # Create standard DRI plot
    p_standard = dri_plot_pre_post(
        ICs_standard,
        "Standard DRI\n$original_F_dims/$original_Pref_dims dimensions"
    )

    # Create dual PCA DRI plot
    p_dual = dri_plot_pre_post(
        ICs,
        "Reduced-Dimension DRI (Dual PCA, $correlation_method)\n$F_dimensions/$Pref_dimensions dimensions, $(round(100 * F_variance_explained, digits=0))/$(round(100 * Pref_variance_explained, digits=0))% σ² explained"
    )

    # Create main title plot
    title_plot = plot(
        title = "DRI Plots Case $case ($case_name): Standard vs Dual PCA\n(min $min_variance_explained variance explained, $correlation_method)",
        grid = false,
        showaxis = false,
        bottom_margin = -100Plots.px
    )

    # Combine all plots
    p = plot(
        title_plot,
        p_standard,
        p_dual,
        layout = @layout([A{0.1h}; B; C]),
        size = (1000, 1000)
    )

    # Save the plot
    outdir = "local-output/dual-pca"
    mkpath(outdir)
    savefig(p, "$outdir/case-$case-dual-pca-$F_dimensions-$Pref_dimensions-dims-$correlation_method.png")

    return DRIs[1], DRIs[2], DRIs_standard[1], DRIs_standard[2], F_dimensions, Pref_dimensions, correlation_method
end

# Main execution
function dual_pca_main()
    # Parse command line arguments
    min_variance_explained = 0.99  # Default value
    correlation_method = :pearson  # Default correlation method
    case_number = nothing          # Default to process all cases
    
    if length(ARGS) > 0
        try
            # Parse min_variance_explained
            min_variance_explained = parse(Float64, ARGS[1])
            if min_variance_explained <= 0.0 || min_variance_explained > 1.0
                error("min_variance_explained must be > 0 and ≤ 1")
            end
            
            # Check for correlation method
            if length(ARGS) > 1
                correlation_method = Symbol(ARGS[2])
                if !(correlation_method in [:cosine_similarity, :pearson, :spearman])
                    error("Invalid correlation method. Must be one of: cosine_similarity, pearson, spearman")
                end
            end
            
            # Check for optional case number
            if length(ARGS) > 2
                case_number = parse(Float64, ARGS[3])
            end
        catch e
            println("Error parsing arguments: ", e)
            println("Usage: julia dual-pca.jl [min_variance_explained] [correlation_method] [case_number]")
            println("  min_variance_explained: Fraction of variance to explain (0 < value ≤ 1). Default: 0.99")
            println("  correlation_method: Correlation method to use (default: pearson)")
            println("    Options: cosine_similarity, pearson, spearman")
            println("  case_number: Optional case number to process (default: all cases)")
            exit(1)
        end
    end

    # Read input data
    data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)
    
    # Get unique cases
    cases = unique(data.CaseID)
    
    # Filter cases if a specific case number is provided
    if case_number !== nothing
        if !(case_number in cases)
            println("Error: Case number $case_number not found in the data")
            exit(1)
        end
        cases = [case_number]
    end
    
    # Create results DataFrame
    results = DataFrame(
        case = Float64[],
        CaseName = String[],
        pre_dri = Float64[],
        post_dri = Float64[],
        delta_dri = Float64[],
        standard_pre_dri = Float64[],
        standard_post_dri = Float64[],
        standard_delta_dri = Float64[],
        min_variance_explained_val = Float64[],
        F_dimensions = Int[],
        Pref_dimensions = Int[],
        correlation_method = Symbol[]
    )
    
    # Process each case
    for case in cases
        case_data = filter(r -> r.CaseID == case, data)
        case_name = case_data[1, :Case]
        
        # Calculate DRI values using both methods
        pre_dri, post_dri, standard_pre_dri, standard_post_dri, F_dimensions, Pref_dimensions, corr_method = dualPCA_plot(case, case_name, min_variance_explained, correlation_method)
        delta_dri = post_dri - pre_dri
        standard_delta_dri = standard_post_dri - standard_pre_dri
        
        # Add to results
        push!(results, (
            case,
            case_name,
            pre_dri,
            post_dri,
            delta_dri,
            standard_pre_dri,
            standard_post_dri,
            standard_delta_dri,
            min_variance_explained,
            F_dimensions,
            Pref_dimensions,
            corr_method
        ))
    end
    
    # Calculate means and add them as a new row
    mean_row = (
        -1.0,  # Special case number for means
        "Means",
        mean(results.pre_dri),
        mean(results.post_dri),
        mean(results.delta_dri),
        mean(results.standard_pre_dri),
        mean(results.standard_post_dri),
        mean(results.standard_delta_dri),
        min_variance_explained,
        round(Int, mean(results.F_dimensions)),  # Round to nearest integer
        round(Int, mean(results.Pref_dimensions)),
        correlation_method
    )
    push!(results, mean_row)
    
    # Save results to CSV
    mkpath("local-output/dual-pca")
    CSV.write("local-output/dual-pca/dual-pca-results-$(min_variance_explained)-minvar-$(correlation_method).csv", results)
    
    # Create horizontal bar chart of delta DRI values
    sort!(results, [:case], by=x -> x, rev=true)
    
    # Find the range of delta values
    min_delta = minimum([minimum(results.delta_dri), minimum(results.standard_delta_dri)])
    max_delta = maximum([maximum(results.delta_dri), maximum(results.standard_delta_dri)])
    
    # Add some padding to the range
    delta_range = max_delta - min_delta
    padding = 0.1 * delta_range
    x_min = min_delta - padding
    x_max = max_delta + padding
    
    # Create the plot
    p = plot(
        size=(1500, 500 + 30 * (nrow(results) + 1)),  # Increased height to accommodate titles and means row
        left_margin=70Plots.mm,
        bottom_margin=10Plots.mm,
        top_margin=10Plots.mm,
        right_margin=10Plots.mm,
        legend=:topright,
        xlims=(x_min, x_max),
        yticks=[],
        layout=@layout([A{0.1h}; [B C D]])  # Title row and three subplots below
    )
    
    # Add main title
    plot!(p[1],
        title="DRI Analysis: Standard vs Dual PCA (min $min_variance_explained variance explained)",
        grid=false,
        showaxis=false,
        bottom_margin=-140Plots.px
    )
    
    # Plot bars in each subplot
    y_positions = 1:1.5:(1.5 * nrow(results))
    bar_width = 0.2
    spacing = 0.3
    
    # Pre-deliberation
    plot!(p[2],
        title="Pre-Deliberation",
        legend=false,
        left_margin=0Plots.mm,
    )
    bar!(p[2],
        y_positions .- spacing,
        results.standard_pre_dri,
        bar_width=bar_width,
        label="Standard Pre DRI",
        color=:blue,
        alpha=0.7,
        orientation=:h
    )
    bar!(p[2],
        y_positions,
        results.pre_dri,
        bar_width=bar_width,
        label="Dual PCA Pre DRI",
        color=:red,
        alpha=0.7,
        orientation=:h
    )
    
    # Post-deliberation
    plot!(p[3],
        title="Post-Deliberation",
        legend=false,
        left_margin=0Plots.mm,
    )
    bar!(p[3],
        y_positions .- spacing,
        results.standard_post_dri,
        bar_width=bar_width,
        label="Standard Post DRI",
        color=:blue,
        alpha=0.7,
        orientation=:h
    )
    bar!(p[3],
        y_positions,
        results.post_dri,
        bar_width=bar_width,
        label="Dual PCA Post DRI",
        color=:red,
        alpha=0.7,
        orientation=:h
    )
    
    # Delta
    plot!(p[4],
        title="Delta",
        legend=:topright,
        left_margin=0Plots.mm,
    )
    bar!(p[4],
        y_positions .- spacing,
        results.standard_delta_dri,
        bar_width=bar_width,
        label="Standard",
        color=:blue,
        alpha=0.7,
        orientation=:h
    )
    bar!(p[4],
        y_positions,
        results.delta_dri,
        bar_width=bar_width,
        label="Dual-PCA",
        color=:red,
        alpha=0.7,
        orientation=:h
    )
    
    # Add case names as text annotations (only on the first subplot)
    for (i, case_name) in enumerate(results.CaseName)
        y_pos = y_positions[i]
        if i == nrow(results)
            # Means row
            annotate!(p[2],
                x_min - 0.01 * delta_range,
                y_pos,
                text("Means", 8, :right)
            )
        else
            annotate!(p[2],
                x_min - 0.01 * delta_range,
                y_pos,
                text("Case $(results.case[i]): $case_name", 8, :right)
            )
        end
    end
    
    # Add vertical lines at x=0 for all subplots
    for i in 2:4
        vline!(p[i], [0], color=:black, linestyle=:dash, label="")
    end
    
    # Axis labels
    xlabel!(p[2], "DRI Value")
    xlabel!(p[3], "DRI Value")
    xlabel!(p[4], "Delta DRI")
    
    # Save the plot
    savefig(p, "local-output/dual-pca/dri-comparison-$(min_variance_explained)-minvar-$(correlation_method).png")
end

# Run main if file is invoked directly
if abspath(PROGRAM_FILE) == @__FILE__
    dual_pca_main()
end

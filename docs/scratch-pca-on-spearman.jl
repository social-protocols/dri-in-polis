using MultivariateStats
using Statistics
using LinearAlgebra
using StatsBase # For rank correlations
using FactorRotations # For Varimax

# --- Previous Code (Setup and PCA on Spearman Correlation) ---

# Sample Data (Users as columns, Ranks as rows)
X = [1 2 1 4;
     3 1 2 3;
     2 4 3 1;
     4 3 4 2]
X = Float64.(X)
num_ranks, num_users = size(X)
num_components = 2 # Or choose based on variance explained

# Standardize each user's ranks (column-wise)
X_means = mean(X, dims=1)
X_stds = std(X, dims=1, corrected=false) # Use population std dev
# Add a small epsilon to avoid division by zero if a user gives identical ranks (though unlikely for ranks)
X_standardized_users = (X .- X_means) ./ (X_stds .+ 1e-8)

println("Standardized User Ranks (X_standardized_users):")
display(X_standardized_users)

# Fit PCA on standardized data (equivalent to PCA on user-user Spearman correlation)
# Use mean=0 because data is already centered column-wise
M_corr = fit(PCA, X_standardized_users; maxoutdim=num_components, mean=0)

println("\nPCA Model (M_corr) based on standardized user ranks:")
display(M_corr)

# Principal Component Scores for each user
Y_corr = transform(M_corr, X_standardized_users)
println("\nTransformed User Scores (Y_corr - components x users):")
display(Y_corr)

# --- New Code (Varimax Rotation and Interpretation) ---

# 1. Extract the original (unrotated) loading matrix
# These loadings show how the rank positions (rows) contribute to the components
original_loadings = projection(M_corr) # Dimensions: num_ranks x num_components
println("\nOriginal PCA Loadings (Rank positions x Components):")
display(original_loadings)

# 2. Apply Varimax Rotation
# Ensure FactorRotations is available (using FactorRotations)
try
    varimax_result = rotate(original_loadings, Varimax())
    rotated_loadings = varimax_result.loadings # Dimensions: num_ranks x num_components

    println("\nRotated Loadings (Varimax):")
    display(rotated_loadings)

    # 3. Identify highest loading items (rank positions) for each rotated component
    println("\nHighest Absolute Loading Item (Rank Position) per Rotated Component:")
    for j in 1:num_components
        component_loadings = rotated_loadings[:, j]
        abs_loadings = abs.(component_loadings)
        max_abs_loading, max_idx = findmax(abs_loadings)

        println("  Rotated Component $j:")
        println("    - Highest loading item (rank position): $max_idx")
        println("    - Loading value (absolute): $max_abs_loading")
        println("    - Loading value (original sign): $(component_loadings[max_idx])")

        # Optionally show top N items
        n_top = 3 # Show top 3 items
        sorted_indices = sortperm(abs_loadings, rev=true)
        println("    - Top $n_top items (rank positions) by absolute loading:")
        for i in 1:min(n_top, num_ranks)
            idx = sorted_indices[i]
            loading_val = component_loadings[idx]
            println("      - Rank Position: $idx, Loading: $loading_val")
        end
    end

catch e
    if isa(e, UndefVarError) && e.var == :FactorRotations
        println("\nError: FactorRotations package not loaded.")
        println("Please install (] add FactorRotations) and ensure 'using FactorRotations' is executed.")
    else
        println("\nAn error occurred during Varimax rotation:")
        showerror(stdout, e)
        println()
    end
end
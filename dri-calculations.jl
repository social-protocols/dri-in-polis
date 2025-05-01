using StatsBase
using Statistics
using DataFrames


function calculate_IC(F_df::DataFrame, Pref_df::DataFrame)
    return calculate_IC(Matrix(F_df), Matrix(Pref_df))
end

function calculate_IC(F::Matrix, Pref::Matrix)
    return calculate_IC(F, Pref, :spearman)
end

function calculate_IC(F::Matrix, Pref::Matrix, method)
    nusers = size(F)[2]

    F_corr = pairwise_correlations(F, method)
    Pref_corr = pairwise_correlations(Pref, method)

    pairs = vcat([[(i,j) for j in 1:(i-1)] for i in 1:nusers]...)

    P1 = getindex.(pairs, 1)
    P2 = getindex.(pairs, 2)
    P_P = string.(P1) .* "-" .* string.(P2) 
    Q = [F_corr[pair[1], pair[2]] for pair in pairs]
    R = [Pref_corr[pair[1], pair[2]] for pair in pairs]

    IC = DataFrame(P_P = P_P, P1 = P1, P2 = P2, Q = Q, R = R)
    
    return IC
end



"""
    calculate_dri(IC::DataFrame)

Calculates the DRI using the IC DataFrame.
Returns: DRI value
"""
function calculate_dri(df::DataFrame; v1::Symbol=:R, v2::Symbol=:Q, center=false)

    x = df[!, v1]
    y = df[!, v2]
    if center
        x = ( x .- mean(x) )
        y = ( y .- mean(y) )
    end

    λ = 1 - sqrt(2)/2  
    diffs = abs.(x .- y) ./ sqrt(2)
    mean_diff = mean(skipmissing(diffs))
    dri = 2 * (((1 - mean_diff) - λ) / (1 - λ)) - 1


    return dri
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


function pairwise_correlations(mat::Matrix, method=:spearman)
    nusers = size(mat)[2]

    return [
        correlation(mat[:, i], mat[:, j], method)
        for i in 1:nusers, j in 1:nusers
    ]

    return result
end

function correlation(x, y, method)

    mask = .!ismissing.(x) .& .!ismissing.(y)

    if count(mask) < 1
        return missing
    end

    # Use the mask to filter both vectors so that they are aligned
    col_i::Vector{Float64} = x[mask]
    col_j::Vector{Float64} = y[mask]

    n_original = length(x)

    # Calculate the distance between two columns
    # This is a placeholder function; replace with actual distance calculation logic
    if method == :cosine_similarity
        return 1 - Distances.cosine_dist(col_i, col_j)
    elseif method == :mwcs
        return magnitiude_weighted_cosine_similarity(col_i, col_j)
    elseif method == :dot_product
        return dot(col_i, col_j)
    elseif method == :euclidean_distance
        return Distances.euclidean(col_i, col_j)
    elseif method == :pearson
        if count(mask) < 2
            return missing
        end
        return cor(col_i, col_j)
    elseif method == :spearman
        if length(x) == 1 || length(y) == 1
            return missing
        end
        return spearman(col_i, col_j)
    elseif method == :subset_rank
        return subset_rank_correlation(n_original, col_i, col_j)
    else 
        error("Unsupported method: $method")
    end

end

function spearman(x, y)
    rx = float.(tiedrank(x))
    ry = float.(tiedrank(y))

    return cor(rx, ry)
end

# A rank correlation metric for the special case where we have a subset of the original ranked items. 
function subset_rank_correlation(n_original, x_ranks::Vector{<:Real}, y_ranks::Vector{<:Real}) 
    if length(x_ranks) == 1
        # Special case: subset has exactly one item
        diff = abs(x_ranks[1] - y_ranks[1])
        return 1.0 - 2.0 * diff / (n_original - 1)
    else
        # General case: Pearson correlation on retained ranks
        return cor(x_ranks, y_ranks)
    end
end

function magnitiude_weighted_cosine_similarity(A,B)

    norm_A = norm(A)
    norm_B = norm(B)

    # --- Handle Edge Cases involving Zero Vectors ---
    # If both norms are zero, they are identical zero vectors.
    if norm_A == 0 && norm_B == 0
        return 1.0
    end
    # If exactly one norm is zero, the magnitude similarity factor min/max would be 0.
    # The overall similarity is thus 0.
    if norm_A == 0 || norm_B == 0
        return 0.0
    end

    # --- Calculate Components for Non-Zero Vectors ---
    # Calculate standard cosine similarity
    # Norms norm_A and norm_B are guaranteed non-zero here.
    dot_product = dot(A, B)
    cos_sim = dot_product / (norm_A * norm_B)

    # Calculate magnitude similarity factor
    # max(norm_A, norm_B) is guaranteed non-zero here.
    mag_sim = min(norm_A, norm_B) / max(norm_A, norm_B)

    # --- Combine ---
    final_similarity = cos_sim * mag_sim

    return final_similarity

end



# ----------------------------------------------------------
# 3. Standardization
# ----------------------------------------------------------
function standardize_rows(mat)
    # mat_std = similar(mat, Float64)
    # for j in 1:size(mat, 2)
    #     col = mat[:, j]
    #     mat_std[:, j] = (col .- mean(col)) ./ (std(col) .+ 1e-8)
    # end
    # return mat_std
    means = mean(mat, dims=2)
    stds = std(mat, dims=2)
    return (mat .- means) ./ (stds .+ 1e-8), means, stds
end

function standardize_columns(mat)
    return (mat .- mean(mat, dims=1)) ./ (std(mat, dims=1) .+ 1e-8)
end




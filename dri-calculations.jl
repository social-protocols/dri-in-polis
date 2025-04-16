using StatsBase
using Statistics


"""
    calculate_IC(F_df::DataFrame, Pref_df::DataFrame)

Calculates the IC DataFrame from F_df and Pref_df by computing pairwise correlations.
Returns: IC DataFrame with Q and R correlations
"""
function calculate_IC(F_df::DataFrame, Pref_df::DataFrame)
    # Compute pairwise correlations
    F_corr, rejects = pairwise_spearman(F_df)
    Pref_corr, rejects = pairwise_spearman(Pref_df)
    
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

function calculate_IC(F::Matrix, Pref::Matrix)
    return calculate_IC(DataFrame(F, :auto), DataFrame(Pref,:auto))
end


"""
    calculate_dri(IC::DataFrame)

Calculates the DRI using the IC DataFrame.
Returns: DRI value
"""
function calculate_dri(df::DataFrame; v1::Symbol=:R, v2::Symbol=:Q)
    λ = 1 - sqrt(2)/2
    diffs = abs.(df[!, v1] .- df[!, v2]) ./ sqrt(2)
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



function pairwise_spearman(df::DataFrame)
    rejected_pairs = 0
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
    return (result, rejected_pairs)
end

function calculate_IC(F::Matrix, Pref::Matrix, method)
    # Compute pairwise correlations
    # F_corr, rejects = pairwise_spearman(F_df)
    # Pref_corr, rejects = pairwise_spearman(Pref_df)

    F_corr = pairwise_correlations(F, method)
    Pref_corr = pairwise_correlations(Pref, method)

    
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

function pairwise_correlations(mat::Matrix, method)
    # cols = names(mat)
    (n,m) = size(mat)
    result = DataFrame(Var1 = String[], Var2 = String[], Freq = Union{Missing,Float64}[])
    for i in 1:m
        for j in 1:(i-1)
            col_i = mat[:, i]
            col_j = mat[:, j]

            cval = correlation(col_i, col_j, method)

            push!(result, (string(i), string(j), cval))
        end
    end
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





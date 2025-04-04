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




using CSV
using DataFrames

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

function remove_all_missing!(df::DataFrame)
    for name in names(df)
        if all(ismissing, df[!, name])
            select!(df, Not(name))
        end
    end
    return df
end

function set_names!(df::DataFrame, newnames::Vector{String})
    rename!(df, Dict(old => new for (old, new) in zip(names(df), newnames)))
end


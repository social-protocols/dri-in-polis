using CSV
using DataFrames
using Statistics
using StatsBase
using LinearAlgebra
using Random
using KernelDensity
using Colors
using Plots

###############################################################################
# 1) DRI Calculation and Plot Functions
###############################################################################

"""
    dri_calc(df::DataFrame, v1::Symbol, v2::Symbol)

Computes the DRI statistic for a given DataFrame `df` that has columns `v1` and `v2`
representing correlations in the "preferences" vs. "considerations" dimension for
all pairs of participants. This is the same formula from your original code.
"""
function dri_calc(df::DataFrame, v1::Symbol, v2::Symbol)
    λ = 1 - (sqrt(2)/2)
    # Compute absolute differences divided by sqrt(2)
    diffs = abs.((df[!, v1] .- df[!, v2]) ./ sqrt(2))
    # Remove missing values
    nonmiss = collect(skipmissing(diffs))
    if isempty(nonmiss)
        return missing  # or some default value if preferred
    else
        m = mean(nonmiss)
        return 2 * (((1 - m) - λ) / (1 - λ)) - 1
    end
end

"""
    correlation(x::Vector{<:Real}, y::Vector{<:Real}, method::Symbol)

Returns a correlation coefficient between x and y for the chosen method:
  :pearson, :spearman, :kendall, or :phi.

- If using :phi, it will treat x and y as binary ±1 data (or 0/1) and compute
  the 2×2 "phi coefficient" (equivalent to Pearson correlation on a 2×2 table).
- For small vectors or missing data, returns missing.
"""
function correlation(x::AbstractVector, y::AbstractVector, method::Symbol)
    mask = .!ismissing.(x) .& .!ismissing.(y)
    if count(mask) <= 1
        return missing
    end
    # Use the mask to filter both vectors so that they are aligned
    xclean = x[mask]
    yclean = y[mask]

    if length(unique(xclean)) == 1 || length(unique(yclean)) == 1
        return missing
    end

    return if method == :pearson
        pearson(xclean, yclean)
    elseif method == :pearson_binary
        pearson_binary(xclean, yclean)
    elseif method == :spearman
        spearman(xclean, yclean)
    elseif method == :phi
        phi(xclean, yclean)
    else
        error("Unknown correlation method: $method")
    end
end

function pearson(xclean, yclean)
    return cor(xclean, yclean)
end

function pearson_binary(xclean, yclean)
    xbin = map(v -> v == -1 ? 0 : v, xclean)
    ybin = map(v -> v == -1 ? 0 : v, yclean)
    return cor(xbin, ybin)
end

function spearman(xclean, yclean)
    rx = float.(tiedrank(xclean))
    ry = float.(tiedrank(yclean))
    return cor(rx, ry)
end

function phi(xclean, yclean)
        # Convert -1 to 0 and keep 1 as is (for binary conversion)
        xbin = map(v -> v == -1 ? 0 : v, xclean)
        ybin = map(v -> v == -1 ? 0 : v, yclean)
        mask2 = [xb in (0, 1) && yb in (0, 1) for (xb, yb) in zip(xbin, ybin)]
        if count(mask2) <= 1
            return missing
        end
        xbin = xbin[mask2]
        ybin = ybin[mask2]
        a = sum(xb == 1 && yb == 1 for (xb, yb) in zip(xbin, ybin))
        b = sum(xb == 1 && yb == 0 for (xb, yb) in zip(xbin, ybin))
        c = sum(xb == 0 && yb == 1 for (xb, yb) in zip(xbin, ybin))
        d = sum(xb == 0 && yb == 0 for (xb, yb) in zip(xbin, ybin))
        denom = sqrt((a+b)*(c+d)*(a+c)*(b+d))
        return denom == 0 ? missing : (a*d - b*c) / denom
end

"""
    pairwise_correlation(df::DataFrame; method=:spearman)

Given a DataFrame `df` whose columns are participants, returns a DataFrame of pairwise
correlations with columns (Var1, Var2, Freq). The correlation is chosen by `method`.
"""
function pairwise_correlation(df::DataFrame; method=:spearman)
    cols = names(df)
    result = DataFrame(Var1 = String[], Var2 = String[], Freq = Union{Missing,Float64}[])
    for i in 1:length(cols)
        for j in 1:(i-1)
            col_i = df[!, cols[i]]
            col_j = df[!, cols[j]]
            cval = correlation(col_i, col_j, method)
            push!(result, (string(cols[i]), string(cols[j]), cval))
        end
    end
    return result
end

"""
    dri_plot(data::DataFrame, x::Symbol, y::Symbol, title_str::String, DRIval::Float64)

Similar to your original ggplot-like scatter:
  - Jittered points in the Q vs. R space
  - Axes from -1.1..1.1
  - Reference lines
  - An annotation with the DRI value
Returns a Plots.jl Plot object.
"""
function dri_plot(data::DataFrame, x::Symbol, y::Symbol, title_str::String, DRIval::Float64)
    xdata = data[!, x]
    ydata = data[!, y]

    # Jitter
    jitter_width = 0.02
    # Random.seed!(123)  # for reproducibility
    x_jitter = xdata .+ jitter_width .* randn(length(xdata))
    y_jitter = ydata .+ jitter_width .* randn(length(ydata))

    p = scatter(
        x_jitter, y_jitter,
        markersize = 2,
        markerstrokewidth=0,
        markercolor = RGBA(0.1, 0, 1, 0.6),
        label = "",
        xlims = (-1.1, 1.1),
        ylims = (-1.1, 1.1),
        legend = false,
        grid = false,
        title = title_str,
        xlabel = "Intersubjective Agreement - Considerations",
        ylabel = "Intersubjective Agreement - Preferences",
        dpi = 300,
    )

    # Add black reference lines
    plot!(p, [-1.1, 1.1], [-1.1, 1.1], color = :black, lw = 1, label = "")
    hline!(p, [0], color = :black, lw = 1, label = "")
    vline!(p, [0], color = :black, lw = 1, label = "")

    # Annotate with DRI
    x_ann = -0.88
    y_ann = 0.88
    annotate!(p, x_ann, y_ann, text("DRI = $(round(DRIval, digits=2))", :red, 13))

    return p
end

###############################################################################
# 2) Main function: single-case pol.is data
###############################################################################

function main_polis_dri(case::String; correlation_method::Symbol)

    votes_file = "Input/polis/$case/votes.csv"

    votes_df = CSV.read(votes_file, DataFrame)

    # Standardize column names
    rename!(votes_df, Dict(Symbol("comment-id") => :comment_id, Symbol("voter-id") => :voter_id, Symbol("vote") => :vote))

    # drop any rows that have missing or zero votes
    df = filter(r -> r.vote in (-1, 1), votes_df)

    if nrow(df) == 0
        throw(ArgumentError("No data for either considerations or preferences. Exiting."))
    end

    # We want each column to be a single participant's votes.
    function pivot_votes_to_matrix(df::DataFrame)
        # We'll do wide pivot with 'comment_id' as row identifier, 'voter_id' as column, 'vote' as value
        unstacked = unstack(df, :comment_id, :voter_id, :vote, combine=last)
        sort!(unstacked, :comment_id)

        # Drop comment ID. Replace missing values with 0.
        select(coalesce.(unstacked, 0), Not(:comment_id))
    end


    df_wide = pivot_votes_to_matrix(df)

    # Currently we don't actually have the 'tag' column. For now, just randomly tag items as considerations or preferences
    # Random.seed!(9362)
    tags = rand(["c", "p"], nrow(df_wide))

    consider_counts = [ sum((df_wide[:,j] .* (tags .== "c")) .!== 0) for j in 1:ncol(df_wide)]
    pref_counts = [ sum((df_wide[:,j] .* (tags .== "p")) .!== 0) for j in 1:ncol(df_wide)]


    # get indices of users that have at least min_answers preferences and min_answers considerations
    consider_indices = findall(consider_counts .>= min_answers)
    pref_indices = findall(pref_counts .>= min_answers)
    usabel_indices = intersect(consider_indices, pref_indices)

    consider_mat = df_wide[findall(tags .== "c"), usabel_indices] 
    pref_mat = df_wide[findall(tags .== "p"), usabel_indices] 
    
    # --- Compute pairwise correlations
    QWrite = pairwise_correlation(consider_mat; method=correlation_method)
    RWrite = pairwise_correlation(pref_mat; method=correlation_method)

    # Both QWrite and RWrite have columns (Var1, Var2, Freq).
    # We want to merge them by (Var1, Var2) so we can store Q in e.g. Q1, R in R1
    # But first, rename the Freq column in each
    rename!(QWrite, :Freq => :Q1)
    rename!(RWrite, :Freq => :R1)

    # We'll do an outerjoin or innerjoin on (Var1, Var2).
    # But remember QWrite only has pairs i>j, and RWrite the same. They should match up.
    # A simple approach is to create a DataFrame with all pairs from QWrite, then leftjoin RWrite.
    IC = leftjoin(QWrite, RWrite, on=[:Var1, :Var2])

    # --- Compute group-level DRI
    group_DRI = dri_calc(IC, :R1, :Q1)
    # println("Group-level DRI = $group_DRI")

    # --- Compute individual-level DRI
    # 1) gather unique participant IDs from all pairs
    partlist = unique(vcat(IC.Var1, IC.Var2))
    # 2) for each participant, filter rows of IC that contain them
    DRIInd = DataFrame(Participant = partlist, DRI = zeros(length(partlist)))
    for (i,p) in enumerate(partlist)
        sub = filter(r -> (r.Var1 == p || r.Var2 == p), IC)
        @show i, p
        if size(sub)[1] > 0
            dri = dri_calc(sub, :R1, :Q1)
            if !ismissing(dri)
                DRIInd[i, :DRI] = dri
            end
        end
    end

    outdir = "Output/polis/$case"

    mkpath("$outdir/Figures")

    # --- Save outputs
    CSV.write("$outdir/IC_polis-$(correlation_method).csv", IC)
    CSV.write("$outdir/DRIInd_polis-$(correlation_method).csv", DRIInd)

    # --- Make one scatter plot for the entire group
    p = dri_plot(IC, :Q1, :R1, "DRI for pol.is data (method=$(correlation_method))", group_DRI)
    savefig(p, "$outdir/Figures/polis_dri_plot-$(correlation_method).png")
    println("Wrote .csv and .png files to $outdir/")
end

###############################################################################
# 3) Run it
###############################################################################

min_answers = 5

if abspath(PROGRAM_FILE) == @__FILE__

    if length(ARGS) < 1
        println("Usage: julia polis_dri.jl CASE [method]")
        println("  CASE: name of the case (e.g. american-assembly.bowling-green)")
        println("  method: correlation method (default=:phi)")
        exit(1)
    end

    # first argument is name of case
    case = ARGS[1]

    correlation_method=:phi
    if length(ARGS) > 1
        correlation_method = Symbol(ARGS[2])
    end

    @show correlation_method


    # Example usage: adapt to your actual filenames
    main_polis_dri(case; correlation_method=correlation_method)
end

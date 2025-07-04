using Random
using CSV
include("dri-plot.jl")
include("prepare-data.jl")

outdir = "docs/validity-test"

function complete_agreement(n, m)
    hcat([1:n for i in 1:m]...)
end

function completely_random(n, m)
     hcat([shuffle(1:n) for i in 1:m]...)
end

function completely_polarized(n, m)
    hcat(
        [1:n for i in 1:m/2]..., 
        [n+1 .- collect(1:n) for i in (m/2+1):m]...
    )
end


function random_diffuse()
    Random.seed!(123)
    nUsers = 100
    nX = 10 
    nY = 10

    # random numbers between

    X = completely_random(nX, nUsers)
    Y = completely_random(nY, nUsers)

    IC = calculate_IC(X, Y)
    DRI = calculate_dri(IC)

    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    dri_plot(p1[1], IC.Q, IC.R, "DRI Plot: Random Rankings\n$nUsers users, $nX considerations, $nY preferences", DRI; show_pearsons=true)
    p1
end


# random but large range
function random_concentrated()
    Random.seed!(123)
    nUsers = 100
    nX = 100 
    nY = 100

    # random numbers between

    X = completely_random(nX, nUsers)
    Y = completely_random(nY, nUsers)


    IC = calculate_IC(X, Y)
    DRI = calculate_dri(IC)

    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    dri_plot(p1[1], IC.Q, IC.R, "DRI Plot: Random Rankings\n$nUsers users, $nX considerations, $nY preferences", DRI; show_pearsons=true)
    p1
end


# random but large range
function random_corresponding_agreement()
    Random.seed!(123)
    nUsers = 100 
    nRandomX = 10
    nRandomY = 10
    nConsensusx = 10
    nConsensusY = 10

    # random numbers between

    # Half random, half complete agreement
    X = [ completely_random(nRandomX, nUsers);
          nRandomX .+ complete_agreement(nConsensusx, nUsers) ]
    Y = [ completely_random(nRandomY, nUsers);
          nRandomY .+ complete_agreement(nConsensusY, nUsers) ]



    IC = calculate_IC(X, Y)
    DRI = calculate_dri(IC)

    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)
    dri_plot(p1[1], IC.Q, IC.R, "DRI Plot: Mix Random / Complete Agreement\n$nUsers users\nConsiderations: $nRandomX random, $nConsensusx complete agreement\nPreferences: $nRandomY random, $nConsensusY complete agreement", DRI; show_pearsons=true)
    p1
end


function concentrated_top()
    Random.seed!(123)
    nUsers = 100
    nX = 100
    nY = 100

    X = completely_random(nX, nUsers)
    Y = complete_agreement(nY, nUsers)

    IC = calculate_IC(X, Y)
    DRI = calculate_dri(IC)

    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    dri_plot(p1[1], IC.Q, IC.R, "DRI Plot: Random Considerations, Complete Agreement on Preferences\n$nUsers users, $nX considerations, $nY preferences", DRI; show_pearsons=true)
    p1
end


function concentrated_top_bottom()
    nUsers = 100
    nX = 100
    nY = 100

    X = hcat([shuffle(1:nX) for i in 1:nUsers]...)
    X = completely_random(nX, nUsers)
    Y = completelyPolarized(nX, nUsers)

    IC = calculate_IC(X, Y)
    DRI = calculate_dri(IC)

    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    dri_plot(p1[1], IC.Q, IC.R, "DRI Plot: Random Considerations, Completely Polarized Preferences\n$nUsers users, $nX considerations, $nY preferences", DRI; show_pearsons=true)
    p1
end


function shuffle_case(data, case)
    Random.seed!(123)

    n_users = size(filter(r -> r.CaseID == case && r.StageID == 1, data))[1]
    shuffle_indices = shuffle(1:n_users)

	ICs = []
    for stage in 1:2
        data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
        F_df, Pref_df = prepare_data(data_filtered)
        IC = calculate_IC(F_df, Pref_df)

        n_pairs = size(IC.R, 1)
        Pref_cor = pairwise_correlations(Matrix(Pref_df)) 
        IC.R = [ Pref_cor[shuffle_indices[IC.P1[i]], shuffle_indices[IC.P2[i]]] for i in 1:n_pairs]

        push!(ICs, IC)
    end
    return ICs
end


function get_ICs(data, case)
    Random.seed!(123)

	ICs = []
    for stage in 1:2
        data_filtered = filter(r -> r.CaseID == case && r.StageID == stage, data)
        F_df, Pref_df = prepare_data(data_filtered)
        IC = calculate_IC(F_df, Pref_df)
        push!(ICs, IC)
    end
    return ICs
end

function frankenstudy(data, case1, case2)
    Random.seed!(123)
    n_users1 = size(filter(r -> r.CaseID == case1 && r.StageID == 1, data))[1]
    n_users2 = size(filter(r -> r.CaseID == case2 && r.StageID == 2, data))[1]
    # If n_users2 < n_users1, we need multiple shuffles to get enough indices
    shuffle_indices = Int[]
    while length(shuffle_indices) < n_users1
        append!(shuffle_indices, shuffle(1:n_users2))
    end
    # Take just the first n_users1 indices
    shuffle_indices = shuffle_indices[1:n_users1]

	ICs = []
    for stage in 1:2
        data_filtered1 = filter(r -> r.CaseID == case1 && r.StageID == stage, data)
        F_df, _ = prepare_data(data_filtered1)

        data_filtered2 = filter(r -> r.CaseID == case2 && r.StageID == stage, data)
        _, Pref_df = prepare_data(data_filtered2)

        Pref_df_sampled = Matrix(Pref_df)[:,shuffle_indices]

	    IC = calculate_IC(Matrix(F_df), Pref_df_sampled)
        push!(ICs, IC)
	end

	return ICs
end


function frankenstudy_method2(data, case1, case2)
    Random.seed!(123)

    ICs = []
    for stage in 1:2
        data_filtered1 = filter(r -> r.CaseID == case1 && r.StageID == stage, data)
        F_df, Pref_df = prepare_data(data_filtered1)
        IC_case1 = calculate_IC(F_df, Pref_df)

        data_filtered2 = filter(r -> r.CaseID == case2 && r.StageID == stage, data)
        F_df, Pref_df = prepare_data(data_filtered2)
        IC_case2 = calculate_IC(F_df, Pref_df)

        n_pairs1 = size(IC_case1)[1]
        n_pairs2 = size(IC_case2)[1]


        shuffle_indices = Int[]
        while length(shuffle_indices) < n_pairs1
            append!(shuffle_indices, shuffle(1:n_pairs2))
        end
        # Take just the first n_pairs1 indices
        shuffle_indices = shuffle_indices[1:n_pairs1]

        IC = copy(IC_case1)
        IC.R = IC_case2.R[shuffle_indices]

        push!(ICs, IC)
    end

    return ICs
end


function shuffled_plot_pre_post(data, case)
    case_name = get_case_name(data, case)
	ICs = shuffle_case(data, case)
	dri_plot_pre_post(ICs, "DRI Plot (Shuffled): $case_name")
end

function regular_plot(data, case, stage)
    case_name = get_case_name(data, case)
    ICs = get_ICs(data, case)
    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    stage_name = stage == 1 ? "Pre-Deliberation" : "Post-Deliberation"
    dri_plot(p1[1], ICs[stage].Q, ICs[stage].R, "DRI Plots: $case_name\n$stage_name", calculate_dri(ICs[stage]))
end


function shuffled_plot(data, case, stage)
    case_name = get_case_name(data, case)
    ICs = shuffle_case(data, case)
    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    stage_name = stage == 1 ? "Pre-Deliberation" : "Post-Deliberation"
    dri_plot(p1[1], ICs[stage].Q, ICs[stage].R, "DRI Plots (Shuffled): $case_name $stage_name", calculate_dri(ICs[stage]))
end

function shuffled_against_standard_plot(data, case, stage)
    case_name = get_case_name(data, case)
    ICs_shuffled = shuffle_case(data, case)
    ICs = get_ICs(data, case)

    stage_name = stage == 1 ? "Pre-Deliberation" : "Post-Deliberation"
    dri_plot_side_by_side([ICs[stage], ICs_shuffled[stage]], "DRI Plot: $case_name $stage_name", ["Actual", "Shuffled"]; show_pearsons=true)
end


function frankenstudy_plot(data, case1, case2)

    case1_name = get_case_name(data, case1)
    case2_name = get_case_name(data, case2)
	ICs = frankenstudy_method2(data, case1, case2)
	dri_plot_pre_post(ICs, "DRI Plots (Cross-Case): $case1_name against $case2_name"; show_pearsons=true)
end

function shuffled_vs_standard_pre_post_plot(data, case)

    case_name = get_case_name(data, case)
    ICs_shuffled = shuffle_case(data, case)
    ICs = get_ICs(data, case)

    DRI_Comparison_Plot(case, case_name, ICs, ICs_shuffled, "Standard v. Shuffled", "Shuffled"; show_pearsons=true)

end

function get_case_name(data, case)
    filter(r -> r.CaseID == case, data)[1, :Case]
end

data = CSV.read("Input/Data1_Raw_Input.csv", DataFrame)


p = random_diffuse()
# save plot
savefig(p, "$outdir/random-diffuse.png")

p = random_concentrated()
# save plot
savefig(p, "$outdir/random-concentrated.png")

p = random_corresponding_agreement()
# save plot
savefig(p, "$outdir/random-corresponding-agreement.png")


p = shuffled_plot(data, 3.0, 2)
# save plot
savefig(p, "$outdir/shuffled-3.0.png")

p = shuffled_plot(data, 18.0, 2)
# save plot
savefig(p, "$outdir/shuffled-18.0.png")


p = regular_plot(data, 3.0, 2)
# save plot
savefig(p, "$outdir/regular-3.0.png")


p = shuffled_against_standard_plot(data, 3.0, 2)
# save plot
savefig(p, "$outdir/shuffled-against-standard-3.0.png")

p = shuffled_against_standard_plot(data, 18.0, 2)
# save plot
savefig(p, "$outdir/shuffled-against-standard-18.0.png")


p = shuffled_against_standard_plot(data, 12.0, 2)
savefig(p, "$outdir/shuffled-against-standard-12.0.png")


p = shuffled_against_standard_plot(data, 13.0, 2)
savefig(p, "$outdir/shuffled-against-standard-13.0.png")


p = shuffled_plot_pre_post(data, 3.0)
# save plot
savefig(p, "$outdir/shuffled-pre-post-3.0.png")

p = shuffled_plot_pre_post(data, 18.0)
# save plot
savefig(p, "$outdir/shuffled-pre-post-18.0.png")

p = frankenstudy_plot(data, 18.0, 3.0)
# save plot
savefig(p, "$outdir/frankenstudy-18.0-3.0.png")

p = shuffled_vs_standard_pre_post_plot(data, 18.0)
savefig(p, "$outdir/shuffled-vs-standard-pre-post-18.0.png")

p = shuffled_vs_standard_pre_post_plot(data, 12.0)
savefig(p, "$outdir/shuffled-vs-standard-pre-post-12.0.png")


p = shuffled_vs_standard_pre_post_plot(data, 3.0)
savefig(p, "$outdir/shuffled-vs-standard-pre-post-3.0.png")




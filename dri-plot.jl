using Plots
using LinearAlgebra

include("dri-calculations.jl")

"""
    dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64; error_bars::Bool=false, show_pearsons::Bool=false)

Creates a single DRI plot in the given subplot with:
  - Scatter points with jitter
  - Reference lines (diagonal, horizontal at 0, vertical at 0)
  - DRI annotation
  - Optional correlation coefficient (ρ) annotation
  - Custom axis labels and title
  - Optional error bars showing standard deviation
"""
function dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64; error_bars::Bool=false, show_pearsons::Bool=false)
    # Add jitter (similar to geom_jitter with width/height = 0.02)
    jitter_width = 0.00
    x_jitter = xdata .+ jitter_width .* randn(length(xdata))
    y_jitter = ydata .+ jitter_width .* randn(length(ydata))

    # Create the scatter plot
    scatter!(p, x_jitter, y_jitter,
        markersize = 1,
        markercolor = RGBA(0.1, 0, 1, 0.6),
        label = "",
        xlims = (-1.1, 1.1), ylims = (-1.1, 1.1),
        legend = false,
        title = title_str,
        titlefontsize = 10,
        grid = false)
    
    # Add reference lines
    plot!(p, [-1.1, 1.1], [-1.1, 1.1], color = :black, lw = 1, label = "")
    hline!(p, [0], color = :black, lw = 1, label = "")
    vline!(p, [0], color = :black, lw = 1, label = "")
    
    # Add DRI annotation (moved further right)
    x_ann = -1.1 + 0.05 * (1.1 - (-1.1))  # Position at 5% from left edge
    y_ann = -1.1 + 0.95 * (1.1 - (-1.1))  # Moved up to 95% from bottom
    annotate!(p, x_ann, y_ann, text("DRI = $(round(dri, digits=2))", :red, 12, :left))
    
    # Add correlation coefficient (ρ) annotation if requested
    if show_pearsons
        rho = cor(xdata, ydata)
        y_ann_rho = y_ann - 0.08 * (1.1 - (-1.1))  # Increased spacing to 8% of plot height
        annotate!(p, x_ann, y_ann_rho, text("ρ = $(round(rho, digits=2))", :red, 12, :left))
    end
    
    # Add error bars if requested
    if error_bars
        x_mean = mean(xdata)
        y_mean = mean(ydata)
        x_std = std(xdata)
        y_std = std(ydata)
        
        # Add horizontal error bar
        plot!(p, [x_mean - x_std, x_mean + x_std], [y_mean, y_mean], 
            color = :red, lw = 1, label = "")
        # Add vertical error bar
        plot!(p, [x_mean, x_mean], [y_mean - y_std, y_mean + y_std], 
            color = :red, lw = 1, label = "")
    end
    
    # Add axis labels
    xlabel!(p, "Intersubjective Agreement - Considerations")
    ylabel!(p, "Intersubjective Agreement - Preferences")
end



"""
Creates a side-by-side plot with pre and post deliberation data, including:
  - Main title for the overall figure
  - Subplot titles for pre and post deliberation
  - Same plotting features as dri_plot for each subplot
  - Optional error bars showing standard deviation
  - Optional correlation coefficient (ρ) annotation
"""
function dri_plot_pre_post(ICs, main_title::String; error_bars::Bool=false, show_pearsons::Bool=false)
    # Call dri_plot_side_by_side with constructed ICs
    dri_plot_side_by_side(ICs, main_title, ["Pre-Deliberation", "Post-Deliberation"]; error_bars=error_bars, show_pearsons=show_pearsons)
end

function dri_plot_side_by_side(ICs, main_title, subtitles; error_bars::Bool=false, show_pearsons::Bool=false)
    # Create a layout with a title plot and two subplots
    title_plot = plot(title = main_title, grid = false, showaxis = false, bottom_margin = -50Plots.px)
    
    # Create the subplots
    subplots = []
    for i in 1:2       
        subplot = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)
        dri_plot(subplot[1], ICs[i].Q, ICs[i].R, subtitles[i], calculate_dri(ICs[i]); error_bars=error_bars, show_pearsons=show_pearsons)
        push!(subplots, subplot)  
    end
    plot(title_plot, subplots[1], subplots[2], layout = @layout([A{0.1h}; [B C]]), size=(800,400))
end



# x_vals = [4,5,6.2]
# y_vals = [10,11,12.2]
# p = plot(title="fo")
# a, b = scatterWithRegression(p, x_vals, y_vals; title="foo")
# p
function scatterWithRegression(p, x_vals, y_vals; title="title")

    m = mean(x_vals, dims=1)
    x_vals_centered = x_vals .- m

    # Add column of ones for intercept
    X = hcat(ones(length(x_vals_centered)), x_vals_centered)
    
    # Calculate coefficients (intercept and slope)
    coeffs = X \ y_vals
    
    # Generate points for regression line across full x range

    minx = min(minimum(x_vals),-1.0)
    miny = min(minimum(y_vals),-1.0)
    maxx = max(maximum(x_vals),1.0)
    maxy = max(maximum(y_vals),1.0)

    x_line = [minx, maxx]
    # y_line = coeffs[1] .+ coeffs[2] .* (x_line)
    y_line = coeffs[1] .+ coeffs[2] .* (x_line .- m)


    # Plot the scatter points
    scatter!(p,
        x_vals,
        y_vals,
        # [s[1] for s in cscores], 
        # [s[2] for s in cscores],
        xlims=(minx,maxx),
        ylims=(miny,maxy),
        title=title
    )

    # Plot regression line
    plot!(p, x_line, y_line,
        color=:black,
        linestyle=:solid,
        linewidth=2,
        label="Regression Line"
    )

    y_hat = coeffs[1] .+ coeffs[2] .* x_vals_centered

    # norm(y_vals .- mean(y_vals))

    return coeffs[1], coeffs[2], norm(y_vals .- y_hat)
end


function DRI_Comparison_Plot(case, case_name, ICs_standard, ICs, title, subtitle; error_bars::Bool=false, show_pearsons::Bool=false)
    # Create standard DRI plot
    p_standard = dri_plot_pre_post(ICs_standard, "Standard DRI"; error_bars=error_bars, show_pearsons=show_pearsons)

    # Create dual PCA DRI plot
    p_dual = dri_plot_pre_post(ICs, subtitle; error_bars=error_bars, show_pearsons=show_pearsons)

    # Create main title plot
    title_plot = plot(
        title = "DRI Plots Case $case ($case_name): $title",
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

    return p
end


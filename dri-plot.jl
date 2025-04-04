using Plots


"""
    dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64)

Creates a single DRI plot in the given subplot with:
  - Scatter points with jitter
  - Reference lines (diagonal, horizontal at 0, vertical at 0)
  - DRI annotation
  - Custom axis labels and title
"""
function dri_plot(p::Plots.Subplot, xdata::Vector, ydata::Vector, title_str::String, dri::Float64)
    # Add jitter (similar to geom_jitter with width/height = 0.02)
    jitter_width = 0.02
    x_jitter = xdata .+ jitter_width .* randn(length(xdata))
    y_jitter = ydata .+ jitter_width .* randn(length(ydata))

    # Create the scatter plot
    scatter!(p, x_jitter, y_jitter,
        markersize = 2,
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
    x_ann = -1.1 + 0.18 * (1.1 - (-1.1))
    y_ann = -1.1 + 0.9 * (1.1 - (-1.1))
    annotate!(p, x_ann, y_ann, text("DRI = $(round(dri, digits=2))", :red, 12))
    
    # Add axis labels
    xlabel!(p, "Intersubjective Agreement - Considerations")
    ylabel!(p, "Intersubjective Agreement - Preferences")
end


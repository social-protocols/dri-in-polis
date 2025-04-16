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


"""
    dri_plot_side_by_side(xdata_pre::Vector, ydata_pre::Vector, xdata_post::Vector, ydata_post::Vector, 
                         main_title::String, dri_pre::Float64, dri_post::Float64)

Creates a side-by-side plot with pre and post deliberation data, including:
  - Main title for the overall figure
  - Subplot titles for pre and post deliberation
  - Same plotting features as dri_plot for each subplot
"""
function dri_plot_side_by_side(xdata_pre::Vector, ydata_pre::Vector, xdata_post::Vector, ydata_post::Vector, 
                             main_title::String, dri_pre::Float64, dri_post::Float64)
    # Create a layout with a title plot and two subplots
    title_plot = plot(title = main_title, grid = false, showaxis = false, bottom_margin = -50Plots.px)
    
    # Create the subplots
    p1 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)
    p2 = plot(layout=(1,1), size=(500,500), dpi=100, margin=5Plots.mm)

    # Create pre-deliberation subplot
    dri_plot(p1[1], xdata_pre, ydata_pre, "Pre-Deliberation", dri_pre)
    
    # Create post-deliberation subplot
    dri_plot(p2[1], xdata_post, ydata_post, "Post-Deliberation", dri_post)
    
    # Combine all plots with the layout, using more space for the title section
    p = plot(title_plot, p1, p2, layout = @layout([A{0.1h}; [B C]]), size=(800,400))
    
    return p
end

#=

    Functions for quick plotting using CairoMakie with pre-defined settings.

=#

module QuickPlots
using CairoMakie
using ..PlotSettings

# plot settings 
CairoMakie.activate!(;px_per_unit=4)
PlotSettings.set_global_themes!()
text_xloc = 0.05
text_yloc = 0.95
xalign = :left
yalign = :top;
fig_width_1, fig_height_1, fig_width_1_2, fig_height_1_2, fig_width_2_2, fig_height_2_2, fontsize, xlabelsize, ylabelsize, xticklabelsize, yticklabelsize, col_gap, xgridvisible, ygridvisible = PlotSettings.load_settings();

default_colors = [:tomato, :royalblue, :aquamarine4, :black, :rebeccapurple, :lightpink2]
default_labels = ["", "", "", "", "", "", ""]
default_linestyles = [:solid, :solid, :solid, :solid, :solid, :solid]
default_linewidths = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
default_alphas = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

alt_linestyles = [:solid, :dash, :dashdot, :dashdotdot, :dot, :dash]
alt_linewidths = [2.0, 1.9, 1.8, 1.7, 1.6, 1.5]
alt_alphas = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75]

function plot11(x, y;
    colors = default_colors,
    labels = default_labels,
    linestyles = default_linestyles,
    linewidths = default_linewidths,
    alphas = default_alphas,
    xlabel = "",
    ylabel = "",
    lim_x_min = nothing,
    lim_x_max = nothing,
    lim_y_min = nothing,
    lim_y_max = nothing,
    xscale = identity,
    yscale = identity,
    text_xloc = nothing, # loacation of text overlayed on plot
    text_yloc = nothing, # loacation of text overlayed on plot
    text = nothing, # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :rt, # position of legend
    labelsize = fontsize, # size of legend and axis labels
    framevisible = true, # whether to draw box around legend
    hlines = nothing, # y-values to draw horizontal lines at in the form [y1, y2, ...]
    vlines = nothing, # x-values to draw vertical lines at in the form [x1, x2, ...]
    save_plot = false, # whether to save the plot to file
    fname = "", # file name to save the plot to
    scatter_lines = false # whether to plot scatter lines instead of normal lines
    )

    with_theme(theme_latexfonts()) do
    f = Figure(size = (fig_width_1, fig_height_1));
    ax = Axis(f[1,1], 
        limits = (lim_x_min, lim_x_max, lim_y_min, lim_y_max),
        xscale = xscale,
        yscale = yscale,
        ylabel=ylabel,
        xlabel=xlabel,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,
        # xticks = xticks,
        # yticks = yticks,
    )


    for i in eachindex(x)
        if scatter_lines
            scatterlines!(ax, x[i], y[i], color=colors[i], label=labels[i])
        else
            lines!(ax, x[i], y[i], color=colors[i], label=labels[i], linestyle=linestyles[i])
        end
    end

    if typeof(text) <: Nothing
        nothing
    else
        text!(ax, text_xloc, text_yloc; text = text, space=:relative, align = (xalign, yalign), fontsize = fontsize)
    end

    if vlines !== nothing
        for vv in vlines
            vlines!(ax, [vv], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if hlines !== nothing
        for hh in hlines
            hlines!(ax, [hh], color = :black, linestyle = :dash, alpha = 0.3)
        end
    end

    if legend
        axislegend(ax, position = position, labelsize = labelsize, framevisible = framevisible)
    end

    resize_to_layout!(f)
    if save_plot
        save(fname, f)
    end
    display(f)
    end
end


end
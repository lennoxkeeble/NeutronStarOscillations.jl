#=

    Module for creating basic animations of time-domain simulation results using CairoMakie.

=#

module Animations
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

@views function compute_y(var, t, time)
    return var[argmin(abs.(t .- time)), :]
end

@views function compute_y_normalized(var, t, time)
    y = var[argmin(abs.(t .- time)), :]
    return y / maximum(abs.(y))
end

@views function compute_y_max(var, t, time)
    return maximum(var[argmin(abs.(t .- time)), :])
end

@views function compute_y_min(var, t, time)
    return minimum(var[argmin(abs.(t .- time)), :])
end

round_up(x::Float64, digits::Int) = ceil(x * 10.0^digits) / 10.0^digits
round_down(x::Float64, digits::Int) = floor(x * 10.0^digits) / 10.0^digits

# label_time(time::Float64)::AbstractString = L"t = %$(round(Int, time))"
label_time(time::Float64; sigdigits::Int64 = 3)::AbstractString = L"t = %$(round(time; sigdigits=sigdigits))\,\mathrm{ms}"

function make_plain_animation(anim_name, timestamps, framerate, fig_width, fig_height, xlabel, ylabel, t, x, y, label1, text_xloc, text_yloc, xalign, yalign, colors, linestyles, alphas, linewidths, labels; labelsize = fontsize, framevisible = true, legend = true, position = :cc, normalize=false, label2 = "", fix_ylims = false, lim_y_min = -100.0, lim_y_max = -200.0)
    with_theme(theme_latexfonts()) do
        f = Figure(size = (fig_width, fig_height))
        ax1 = Axis(f[1, 1], ylabel=ylabel, xlabel = xlabel,
        xlabelsize = ylabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xticksvisible = true, yticksvisible = true,
        xticklabelsvisible = true, yticklabelsvisible = true,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,)
        
        time = Observable(0.0)
        
        lines!(ax1, x, @lift(compute_y(y, t, $time)), color=colors[1], linestyle=linestyles[1], alpha=alphas[1], linewidth=linewidths[1], label=labels[1])
        text!(ax1, text_xloc, text_yloc; text =  @lift(label_time($time; sigdigits = 2)), space=:relative,
            align = (xalign, yalign), fontsize = fontsize)
        isequal(label1, "") ? nothing : text!(ax1, text_xloc, text_yloc - 0.07; text = label1, space=:relative,
            align = (xalign, yalign), fontsize = fontsize)
        if label2 != ""
            text!(ax1, text_xloc, text_yloc - 0.15; text = label2, space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
        end

        if legend
            # axislegend(ax1, position = position, labelsize = labelsize, framevisible = framevisible)
            # Add a centered x-label in a new row below both axes
            Legend(f[0,1], ax1, position = position, labelsize = labelsize, framevisible = framevisible, orientation=:horizontal)
        end
        resize_to_layout!(f)

        record(f, anim_name, timestamps; framerate = framerate) do tt
            time[] = tt

            ymin = compute_y_min(y, t, tt)
            ymax = compute_y_max(y, t, tt)
            ylims = [ymin, ymax]
            y_offset_factor = 0.05
            y_range = ylims[2] - ylims[1]
            if y_range == 0.0
                y_range = 1e-12
            end
            ylims = [ylims[1] - y_offset_factor * y_range, ylims[2] + y_offset_factor * y_range]
            fix_ylims ? ylims!(ax1, lim_y_min, lim_y_max) : ylims!(ax1, ylims[1], ylims[2])
        end
    end

end

function make_plain_animation(anim_name, timestamps, framerate, fig_width, fig_height, xlabel, ylabel, t, x, y1, y2, y3, label1, text_xloc, text_yloc, xalign, yalign, colors, linestyles, alphas, linewidths, labels; labelsize = fontsize, framevisible = true, legend = true, position = :cc, normalize=false, label2 = "", fix_ylims = false, ymin = -100.0, ymax = -200.0)
    with_theme(theme_latexfonts()) do
        f = Figure(size = (fig_width, fig_height))
        ax1 = Axis(f[1, 1], ylabel=ylabel, xlabel = xlabel,
        xlabelsize = ylabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xticksvisible = true, yticksvisible = true,
        xticklabelsvisible = true, yticklabelsvisible = true,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,)
        
        time = Observable(0.0)
        
        lines!(ax1, x, @lift(compute_y(y1, t, $time)), color=colors[1], linestyle=linestyles[1], alpha=alphas[1], linewidth=linewidths[1], label=labels[1])
        lines!(ax1, x, @lift(compute_y(y2, t, $time)), color=colors[2], linestyle=linestyles[2], alpha=alphas[2], linewidth=linewidths[2], label=labels[2])
        lines!(ax1, x, @lift(compute_y(y3, t, $time)), color=colors[3], linestyle=linestyles[3], alpha=alphas[3], linewidth=linewidths[3], label=labels[3])
        text!(ax1, text_xloc, text_yloc; text =  @lift(label_time($time; sigdigits = 2)), space=:relative,
            align = (xalign, yalign), fontsize = fontsize)
        isequal(label1, "") ? nothing : text!(ax1, text_xloc, text_yloc - 0.07; text = label1, space=:relative,
            align = (xalign, yalign), fontsize = fontsize)
        if label2 != ""
            text!(ax1, text_xloc, text_yloc - 0.15; text = label2, space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
        end

        if legend
            # axislegend(ax1, position = position, labelsize = labelsize, framevisible = framevisible)
            # Add a centered x-label in a new row below both axes
            Legend(f[0,1], ax1, position = position, labelsize = labelsize, framevisible = framevisible, orientation=:horizontal)
        end
        resize_to_layout!(f)

        record(f, anim_name, timestamps; framerate = framerate) do tt
            time[] = tt

            ymin = minimum([compute_y_min(y1, t, tt), compute_y_min(y2, t, tt), compute_y_min(y3, t, tt)])
            ymax = maximum([compute_y_max(y1, t, tt), compute_y_max(y2, t, tt), compute_y_max(y3, t, tt)])
            ylims = [ymin, ymax]
            y_offset_factor = 0.05
            y_range = ylims[2] - ylims[1]
            if y_range == 0.0
                y_range = 1e-12
            end
            ylims = [ylims[1] - y_offset_factor * y_range, ylims[2] + y_offset_factor * y_range]
            fix_ylims ? ylims!(ax1, ymin, ymax) : ylims!(ax1, ylims[1], ylims[2])
        end
    end

end



function make_animation(anim_name, timestamps, framerate, fig_width, fig_height, xlabel, ylabel, t, x, y, label1, text_xloc, text_yloc, xalign, yalign; color=:tomato, linestyle=:solid, alpha=0.8, linewidth=1.5, normalize=false, label2 = "")
    with_theme(theme_latexfonts()) do
        f = Figure(size = (fig_width, fig_height))
        ax1 = Axis(f[1, 1], ylabel=ylabel, xlabel = xlabel,
        xlabelsize = ylabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
        xticksvisible = true, yticksvisible = true,
        xticklabelsvisible = true, yticklabelsvisible = true,
        xgridvisible = xgridvisible,
        ygridvisible = ygridvisible,
        ytickformat = values -> ["$(round(value, digits=2))" for value in values],
        xtickformat = values -> ["$(round(Int, value))" for value in values])
        
        time = Observable(0.0)
        y_offset_factor = 0.05
        if normalize
            lines!(ax1, x, @lift(compute_y_normalized(y, t, $time)), color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)
            text!(ax1, text_xloc, text_yloc; text =  @lift(label_time($time; sigdigits = 2)), space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
            text!(ax1, text_xloc, text_yloc - 0.07; text = label1, space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
            if label2 != ""
                text!(ax1, text_xloc, text_yloc - 0.20; text = label2, space=:relative,
                    align = (xalign, yalign), fontsize = fontsize)
            end
            resize_to_layout!(f)

            record(f, anim_name, timestamps; framerate = framerate) do t
                time[] = t
                xlims = (0, x[end])
                ylims = [-1.1, 1.1]

                xlims!(ax1, xlims[1], xlims[2])
                ylims!(ax1, ylims[1], ylims[2])
                
                ax1.xticks = range(xlims[1], xlims[2], length=2)
                ax1.yticks = [-1, -0.5, 0.0, 0.5, 1.0]
            end
        else
            lines!(ax1, x, @lift(compute_y(y, t, $time)), color=color, linestyle=linestyle, alpha=alpha, linewidth=linewidth)
            text!(ax1, text_xloc, text_yloc; text =  @lift(label_time($time; sigdigits = 2)), space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
            text!(ax1, text_xloc, text_yloc - 0.07; text = label1, space=:relative,
                align = (xalign, yalign), fontsize = fontsize)
            if label2 != ""
                text!(ax1, text_xloc, text_yloc - 0.20; text = label2, space=:relative,
                    align = (xalign, yalign), fontsize = fontsize)
            end
            resize_to_layout!(f)

            record(f, anim_name, timestamps; framerate = framerate) do tt
                time[] = tt
                xlims = (0, x[end])
                ylims = [-compute_y_max(y, t, tt), compute_y_max(y, t, tt)]
                y_range = ylims[2] - ylims[1]
                ylims = [ylims[1] - y_offset_factor * y_range, ylims[2] + y_offset_factor * y_range]

                xlims!(ax1, xlims[1], xlims[2])
                ylims!(ax1, ylims[1], ylims[2])
                
                ax1.xticks = range(xlims[1], xlims[2], length=2)
                ax1.yticks = range(ylims[1], ylims[2], length=2)
            end
        end
    end

    end
end
#=

    Global plot settings used for plotting with CairoMakie.
=#

module PlotSettings
using CairoMakie

const fig_width_1 = 500;
const fig_height_1 = 400;
const fig_width_1_2 = fig_width_1 * 2
const fig_height_1_2 = fig_height_1
const fig_width_2_2 = fig_width_1 * 2
const fig_height_2_2 = fig_height_1 * 2

const fontsize = 22
const xlabelsize = 32   # font size for x-axis label
const ylabelsize = 32   # font size for y-axis label
const xticklabelsize = 22  # font size for x-axis tick labels
const yticklabelsize = 22   # font size for y-axis tick labels
const col_gap = 40.0
const xgridvisible = false
const ygridvisible = false

function set_global_themes!()
    fontsize_theme = Theme(fontsize = PlotSettings.fontsize)
    update_theme!(fontsize_theme)
    rasterize_theme = Theme(rasterize = true)
    set_theme!(rasterize_theme)
end

function load_settings()
    set_global_themes!()
    return fig_width_1, fig_height_1, fig_width_1_2, fig_height_1_2, fig_width_2_2, fig_height_2_2, fontsize, xlabelsize, ylabelsize, xticklabelsize, yticklabelsize, col_gap, xgridvisible, ygridvisible
end

end
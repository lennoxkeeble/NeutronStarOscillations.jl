#=

    Module comprising of functions which solve the perfect fluid and BDNK equations of motion governing linear radial perturbations of polytropic neutron stars in the
    time domain in the Cowling approximation. (Note we don't provide separate perfect fluid functions; one takes η = ζ = 0 in the Eckart functions). Both sets of equations are
    integrated in time using an explicit RK4 scheme. Both numerical schemes are second-order accurate since we use
    second-order finite differences for spatial derivatives.

=#

module CowlingTimeDomain
using LaTeXStrings
using CairoMakie
using ..PlotSettings
using ..QuickPlots
using ..Animations
using NeutronStarOscillations
using Printf
using HDF5


one_norm(x::AbstractArray) = sum(abs, x) / length(x)
function td_fname(star::NeutronStarOscillations.Star, h::Float64)::String
    if star.mode == -1
        return star.data_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h);
    else
        return star.data_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_mode_%s_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, star.mode, h);
    end
end

function plot_fname(star::NeutronStarOscillations.Star, h::Float64, var::String, time::Float64)::String
    if star.mode == -1
        return star.fig_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e_%s_time_%.2f.png", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h, var, time);
    else
        return star.fig_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_mode_%s_h_%.1e_%s_time_%.2f.png", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, star.mode, h, var, time);
    end
end

function animation_fname(star::NeutronStarOscillations.Star, h::Float64, var::String)::String
    if star.mode == -1
        return star.fig_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e_%s.mp4", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h, var);
    else
        return star.fig_path * @sprintf("Cowling_time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_mode_%s_h_%.1e_%s.mp4", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, star.mode, h, var);
    end
end

function convergence_plot_fname(star::NeutronStarOscillations.Star, var::String)::String
    if star.mode == -1
        return star.fig_path * @sprintf("Cowling_time_domain_%s_convergence_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s.png", var, star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T);
    else
        return star.fig_path * @sprintf("Cowling_time_domain_%s_convergence_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_mode_%s.png", var, star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, star.mode);
    end
end


function load_solution(star::NeutronStarOscillations.Star, h::Float64)
    filename = td_fname(star, h);
    file = h5open(filename, "r")
    finalizer(file) do f
        close(f)  # Automatically close when the file object is garbage collected
    end
    return file 
end

label_time(time::Float64; sigdigits::Int64 = 3)::AbstractString = L"t = %$(round(time; sigdigits=sigdigits))\,\mathrm{ms}"
get_time_idx(time::Float64, t::Vector{Float64})::Int64 = argmin(@. (abs(t - time)))
const static_vars::Vector{String} = ["cs", "cs_prime", "m", "p", "ε", "ν"];
const eckart_IRs::Vector{String} = ["IR"]; # independent residuals for Eckart stars
const BDNK_IRs::Vector{String} = ["IR1", "IR2"]; # independent residuals for BDNK stars

# plot single solution
function plot_var(star::NeutronStarOscillations.Star, h::Float64, var::String, time::Float64;
    desample_factor::Int64=1,
    annotate_time::Bool=true,
    colors = NeutronStarOscillations.QuickPlots.default_colors,
    labels = NeutronStarOscillations.QuickPlots.default_labels,
    linestyles = NeutronStarOscillations.QuickPlots.default_linestyles,
    linewidths = NeutronStarOscillations.QuickPlots.default_linewidths,
    alphas = NeutronStarOscillations.QuickPlots.default_alphas,
    xlabel = "",
    ylabel = "",
    lim_x_min = nothing,
    lim_x_max = nothing,
    lim_y_min = nothing,
    lim_y_max = nothing,
    yscale = identity,
    text_xloc = nothing, # loacation of text overlayed on plot
    text_yloc = nothing, # loacation of text overlayed on plot
    text = nothing, # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :rt, # position of legend
    labelsize = NeutronStarOscillations.QuickPlots.fontsize, # size of legend and axis labels
    framevisible = true, # whether to draw box around legend
    hlines = nothing, # y-values to draw horizontal lines at in the form [y1, y2, ...]
    vlines = nothing, # x-values to draw vertical lines at in the form [x1, x2, ...]
    save_plot = false, # whether to save the plot to file
    scatter_lines = false # whether to plot scatter lines instead of normal lines
    )
    
    sol = load_solution(star, h);
    if !(var ∈ keys(sol["solution"]))
        throw(ArgumentError("Variable $var not found in solution file! Available variables are: $(keys(sol["solution"]))"))
    end

    t = sol["solution/t"][:];
    x = [sol["solution/r"][1:desample_factor:end]];
    if var ∈ static_vars
        y = [sol["solution"][var][1:desample_factor:end]];
        plot_text = ""
    else
        idx = get_time_idx(time, t);
        y = [sol["solution"][var][idx, 1:desample_factor:end]];
        true_time = t[idx];
        plot_text = label_time(true_time; sigdigits = 3)
    end

    close(sol)

    if yscale == log10
        for i in eachindex(y)
            y[i] = abs.(y[i])
        end
    end


    if annotate_time
        text_xloc = 0.05;
        text_yloc = 0.95;
        text = plot_text;
        xalign = :left;
        yalign = :top;
    end

    fname = plot_fname(star, h, var, true_time);

    NeutronStarOscillations.QuickPlots.plot11(x, y;
    colors = colors,
    labels = labels,
    linestyles = linestyles,
    linewidths = linewidths,
    alphas = alphas,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    yscale = yscale,
    text_xloc = text_xloc,
    text_yloc = text_yloc,
    text = text,
    xalign = xalign,
    yalign = yalign,
    legend = legend,
    position = position,
    labelsize = labelsize,
    framevisible = framevisible,
    hlines = hlines,
    vlines = vlines,
    save_plot = save_plot,
    fname = fname, 
    scatter_lines = scatter_lines
    )
end

# plot one norm of independent residuals
function plot_IR_one_norms(star::NeutronStarOscillations.Star, h::Vector{Float64}, var::String, time::Float64;
    desample_factor::Int64=1,
    colors = NeutronStarOscillations.QuickPlots.default_colors,
    labels = NeutronStarOscillations.QuickPlots.default_labels,
    linestyles = NeutronStarOscillations.QuickPlots.default_linestyles,
    linewidths = NeutronStarOscillations.QuickPlots.default_linewidths,
    alphas = NeutronStarOscillations.QuickPlots.default_alphas,
    xlabel = "",
    ylabel = "",
    lim_x_min = nothing,
    lim_x_max = nothing,
    lim_y_min = nothing,
    lim_y_max = nothing,
    yscale = identity,
    text_xloc = nothing, # loacation of text overlayed on plot
    text_yloc = nothing, # loacation of text overlayed on plot
    text = nothing, # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :rt, # position of legend
    labelsize = NeutronStarOscillations.QuickPlots.fontsize, # size of legend and axis labels
    framevisible = true, # whether to draw box around legend
    hlines = nothing, # y-values to draw horizontal lines at in the form [y1, y2, ...]
    vlines = nothing, # x-values to draw vertical lines at in the form [x1, x2, ...]
    save_plot = false, # whether to save the plot to file
    scatter_lines = false # whether to plot scatter lines instead of normal lines
    )

    if star.τε == star.τP == star.τQ == 0.0
        if !(var ∈ eckart_IRs)
            throw(ArgumentError("Invalid independent residual for Eckart star. Choose from: $(eckart_IRs)"))
        end
    else
        if !(var ∈ BDNK_IRs)
            throw(ArgumentError("Invalid independent residual for BDNK star. Choose from: $(BDNK_IRs)"))
        end
    end
    
    sols = [load_solution(star, hh) for hh in h];


    x = [sol["solution/t"][:] for sol in sols];
    y = Vector{Float64}[];

    for i in eachindex(sols)
        tt = sols[i]["solution/t"][:];
        yy = zero(tt)
        for j in eachindex(tt)
            yy[j] = one_norm(sols[i]["solution"][var][j, :])
        end
        push!(y, yy);
    end

    if desample_factor != 1
        for i in eachindex(x)
            x[i] = x[i][1:desample_factor:end];
            y[i] = y[i][1:desample_factor:end];
        end
    end

    fname = convergence_plot_fname(star, var);

    NeutronStarOscillations.QuickPlots.plot11(x, y;
    colors = colors,
    labels = labels,
    linestyles = linestyles,
    linewidths = linewidths,
    alphas = alphas,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    yscale = yscale,
    text_xloc = text_xloc,
    text_yloc = text_yloc,
    text = text,
    xalign = xalign,
    yalign = yalign,
    legend = legend,
    position = position,
    labelsize = labelsize,
    framevisible = framevisible,
    hlines = hlines,
    vlines = vlines,
    save_plot = save_plot,
    fname = fname, 
    scatter_lines = scatter_lines
    )
end

# animate single solution
function animate_var(star::NeutronStarOscillations.Star, h::Float64, var::String;
    desample_factor::Int64=1,
    stop_time::Float64=0.1, # stop time of animation in ms
    animation_length::Float64=10.0, # length of animation in seconds
    framerate::Int64=10, # frames per second
    colors = NeutronStarOscillations.QuickPlots.default_colors,
    labels = NeutronStarOscillations.QuickPlots.default_labels,
    linestyles = NeutronStarOscillations.QuickPlots.default_linestyles,
    linewidths = NeutronStarOscillations.QuickPlots.default_linewidths,
    alphas = NeutronStarOscillations.QuickPlots.default_alphas,
    xlabel = "",
    ylabel = "",
    lim_y_min = nothing,
    lim_y_max = nothing,
    fix_ylims = false,
    yscale = identity,
    text_xloc = 0.05, # loacation of text overlayed on plot
    text_yloc = 0.95, # loacation of text overlayed on plot
    text = "", # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :cc, # position of legend
    width = 1.3 * NeutronStarOscillations.QuickPlots.fig_width_1,
    height = 1.3 * NeutronStarOscillations.QuickPlots.fig_height_1,
    )
    
    sol = load_solution(star, h);
    if !(var ∈ keys(sol["solution"]))
        throw(ArgumentError("Variable $var not found in solution file! Available variables are: $(keys(sol["solution"]))"))
    end

    anim_fname = animation_fname(star, h, var)

    start_time = 0.0;
    num_plots = Int(animation_length * framerate)

    times = range(start = start_time, stop = stop_time, length = num_plots) |> collect
    timestamps = times

    t = sol["solution/t"][:];
    x = sol["solution/r"][1:desample_factor:end];
    y = sol["solution/"*var][:, 1:desample_factor:end];

    Animations.make_plain_animation(anim_fname, timestamps, framerate, width, height,
        xlabel, ylabel, t, x, y, text, text_xloc, text_yloc, xalign, yalign, colors,
        linestyles, alphas, linewidths, labels; legend = legend, position = position, fix_ylims = fix_ylims, lim_y_min = lim_y_min, lim_y_max = lim_y_max)
end

function plot_characteristic_speeds(star::NeutronStarOscillations.Star, h::Float64)
    if star.τε == star.τP == star.τQ == 0.0
        throw(ArgumentError("Characteristic speeds can only be computed for BDNK stars!"))
    end
    
    sol = load_solution(star, h);
    Λ0 = sol["solution/Λ0"][:];
    Λ1 = sol["solution/Λ1"][:];
    Λ2 = sol["solution/Λ2"][:];
    cp2 = zero(Λ0);
    cm2 = zero(Λ0);
    try
        cp2 .= @. (Λ1 + sqrt(Λ1^2 - Λ0)) / (Λ2);
        cm2 .= @. (Λ1 - sqrt(Λ1^2 - Λ0)) / (Λ2);
    catch e
        if e isa DomainError
            throw(DomainError("Characteristic speeds are complex ——— an acausal frame has been chosen!"))
        else
            rethrow(e)
        end
    end
    
    x = [sol["solution/r"][:], sol["solution/r"][:]];
    y = [cm2, cp2];
    labels = [L"c_{-}^{2}", L"c_{+}^{2}"];

    xlabel = L"r\,[\mathrm{km}]";
    ylabel = L"c^{2}";
    
    close(sol)

    NeutronStarOscillations.QuickPlots.plot11(
        x, y;
        labels = labels,
        xlabel = xlabel,
        ylabel = ylabel,
        legend = true,
        position = :lc,
        framevisible = true,
        )
end

module PerfectFluid
using Printf
using HDF5
using NeutronStarOscillations
using ...HDF5Helper
using ...TOV
using ..CowlingTimeDomain
using ...CowlingTDEqs
using ...FiniteDiffOrder4

const c = 299792458;
const sec_to_km = c * 1e-3;
# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function initialize_solution_file(filename::String, m::AbstractVector{Float64}, p::AbstractVector{Float64}, ε::AbstractVector{Float64}, ν::AbstractVector{Float64}, r::AbstractVector{Float64}, cs::AbstractVector{Float64}, cs_prime::AbstractVector{Float64}, u1::AbstractVector{Float64}, v1::AbstractVector{Float64}, ds::Int64, N::Int, num_time_steps::Int, chunk_size::Int)
    fmode = "w"
    file = h5open(filename, fmode); # filename should include the path

    # create group for solution data
    solution_group_name = "solution"
    HDF5Helper.Cowling.create_file_group!(file, solution_group_name);

    # first store simulation parameters
    # param_group_name = "parameters"
    # HDF5Helper.Cowling.create_file_group!(file, param_group_name);


    # save data fixed at run time
    file["solution"]["r"] = r[1:ds:end];
    file["solution"]["m"] = m[1:ds:end];
    file["solution"]["p"] = p[1:ds:end];
    file["solution"]["ε"] = ε[1:ds:end];
    file["solution"]["ν"] = ν[1:ds:end];
    file["solution"]["cs"] = cs[1:ds:end];
    file["solution"]["cs_prime"] = cs_prime[1:ds:end];

    # create datasets for the solution
    dataspace = ((num_time_steps+1, N), (num_time_steps+1, N))
    chunk = (chunk_size, N)
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "xi", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "xi_dot", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "IR", Float64, dataspace, chunk);


    # manually save initial data
    file[solution_group_name]["xi"][1, :] = u1[1:ds:end];
    file[solution_group_name]["xi_dot"][1, :] = v1[1:ds:end];
    file[solution_group_name]["IR"][1, :] = zero(u1[1:ds:end]) * NaN; # don't compute IR at initial time
    return file
end

@views function update_solution_array!(sol::AbstractArray{Float64}, u1::Vector{Float64}, v1::Vector{Float64}, IR::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = u1[1:ds_fact:nMax];
    sol[time_idx, 2, :] = v1[1:ds_fact:nMax];
    sol[time_idx, 3, :] = IR[1:ds_fact:nMax];
end

@views function update_solution_array!(sol::AbstractArray{Float64}, sol_three_levels::AbstractArray{Float64}, IR::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = sol_three_levels[2, 1, 1:ds_fact:nMax];
    sol[time_idx, 2, :] = sol_three_levels[2, 2, 1:ds_fact:nMax];
    sol[time_idx, 3, :] = IR[1:ds_fact:nMax];
end

function load_solution(filename::String)
    file = h5open(filename, "r")
    finalizer(file) do f
        close(f)  # Automatically close when the file object is garbage collected
    end
    return file 
end

# compute the k1, k2, k3, k4 terms for the RK4 time integration
@views function compute_ks!(
    u1::Vector{Float64}, v1::Vector{Float64}, 
    u1_k1::Vector{Float64}, v1_k1::Vector{Float64},
    u1_k2::Vector{Float64}, v1_k2::Vector{Float64},
    u1_k3::Vector{Float64}, v1_k3::Vector{Float64},
    k1::Matrix{Float64}, k2::Matrix{Float64}, k3::Matrix{Float64}, k4::Matrix{Float64},
    m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64},
    r::Vector{Float64}, h::Float64, k::Float64, nPoints::Int64)
    ### BCs at r=0 for regularity: ξ[t,0] = u1[t, 0] = 0, ∂_{t}ξ[t,0] = v1[t, 0] = 0, so the equations are already solved at r = 0 ###

    ## COMPUTE K1 ##
    @inbounds for i = 2:(nPoints-1)
        k1[1, i] = k * v1[i]
        k1[2, i] = k * CowlingTDEqs.RK4.PerfectFluid.dv1_dt(i, u1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

        u1_k1[i] = u1[i] + 0.5 * k1[1, i]
        v1_k1[i] = v1[i] + 0.5 * k1[2, i]
    end

    i = nPoints
    k1[1, i] = k * v1[i]
    k1[2, i] = k * CowlingTDEqs.RK4Surface.PerfectFluid.dv1_dt(i, u1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

    u1_k1[i] = u1[i] + 0.5 * k1[1, i]
    v1_k1[i] = v1[i] + 0.5 * k1[2, i]

    ## COMPUTE K2 ##
    @inbounds for i = 2:(nPoints-1)
        k2[1, i] = k * v1_k1[i]
        k2[2, i] = k * CowlingTDEqs.RK4.PerfectFluid.dv1_dt(i, u1_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

        u1_k2[i] = u1[i] + 0.5 * k2[1, i]
        v1_k2[i] = v1[i] + 0.5 * k2[2, i]
    end

    i = nPoints
    k2[1, i] = k * v1_k1[i]
    k2[2, i] = k * CowlingTDEqs.RK4Surface.PerfectFluid.dv1_dt(i, u1_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

    u1_k2[i] = u1[i] + 0.5 * k2[1, i]
    v1_k2[i] = v1[i] + 0.5 * k2[2, i]


    ## COMPUTE K3 ##
    @inbounds for i = 2:(nPoints-1)
        k3[1, i] = k * v1_k2[i]
        k3[2, i] = k * CowlingTDEqs.RK4.PerfectFluid.dv1_dt(i, u1_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

        u1_k3[i] = u1[i] + k3[1, i]
        v1_k3[i] = v1[i] + k3[2, i]
    end

    i = nPoints
    k3[1, i] = k * v1_k2[i]
    k3[2, i] = k * CowlingTDEqs.RK4Surface.PerfectFluid.dv1_dt(i, u1_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)

    u1_k3[i] = u1[i] + k3[1, i]
    v1_k3[i] = v1[i] + k3[2, i]

    ## COMPUTE K4 ##
    @inbounds for i = 2:(nPoints-1)
        k4[1, i] = k * v1_k3[i]
        k4[2, i] = k * CowlingTDEqs.RK4.PerfectFluid.dv1_dt(i, u1_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h)
    end

    i = nPoints
    k4[1, i] = k * v1_k3[i]
    k4[2, i] = k * CowlingTDEqs.RK4Surface.PerfectFluid.dv1_dt(i, u1_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h) 
end

# compute independent residual
@views function compute_IR!(sol_three_levels::AbstractArray, IR::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, v1_tm1::Vector{Float64}, u1_t0::Vector{Float64}, v1_t0::Vector{Float64}, v1_tp1::Vector{Float64}, h::Float64, k::Float64, nPoints::Int64)
    # fill arrays
    v1_tm1 .= sol_three_levels[1, 2, :];

    u1_t0 .= sol_three_levels[2, 1, :];
    v1_t0 .= sol_three_levels[2, 2, :];

    v1_tp1 .= sol_three_levels[3, 2, :];

    # equations solved exactly at center
    IR[1] = 0.0

    for i = 2:nPoints-1
        IR[i] = CowlingTDEqs.Leapfrog.PerfectFluid.Eq1(i, v1_tm1, v1_tp1, u1_t0, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h, k)
    end
    i = nPoints
    IR[i] = CowlingTDEqs.LeapfrogSurface.PerfectFluid.Eq1(i, v1_tm1, v1_tp1, u1_t0, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], h, k)
end

# apply Kreiss-Oliger dissipation
function kreiss_oliger(var::AbstractVector{Float64}, j::Int, coef::Float64, nPoints::Int64)
    if j == 1
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[1]
        u_plus_1 = var[2]
        u_plus_2 = var[3]
    elseif j == 2
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[2]
        u_plus_1 = var[3]
        u_plus_2 = var[4]
    elseif j == nPoints-1
        u_minus_2 = var[nPoints-3]
        u_minus_1 = var[nPoints-2]
        u = var[nPoints-1]
        u_plus_1 = var[nPoints]
        u_plus_2 = 0.0
    elseif j == nPoints
        u_minus_2 = var[nPoints-2]
        u_minus_1 = var[nPoints-1]
        u = var[nPoints]
        u_plus_1 = 0.0
        u_plus_2 = 0.0
    else
        u_minus_2 = var[j-2]
        u_minus_1 = var[j-1]
        u = var[j]
        u_plus_1 = var[j+1]
        u_plus_2 = var[j+2]
    end

    return coef * (u_plus_2 - 4.0 * u_plus_1 + 6.0 * u - 4.0 * u_minus_1 + u_minus_2) / 16.0
end

function solve(star::NeutronStarOscillations.Star, h::Float64; print_progress::Bool=true)
    fname = CowlingTimeDomain.td_fname(star, h)
    if abs(star.ξ_ID(0.0)) > 1e-16
        error("|ξ(r=0)| > 1e-16. ξ must be zero at the center by regularity. ξ(0.0) = $(star.ξ_ID(0.0))")
    end

    # convert time from ms to km
    km_to_ms = 1e3 / (sec_to_km)
    target_dt_save = star.dt_save / km_to_ms
    total_time = star.T / km_to_ms

    k = star.CFL * h; # time step from CFL condition
    total_num_time_steps = Int(ceil(total_time / k));
    time_steps = range(start = k, step = k, stop = total_time) |> collect
    dt_save_idx = argmin(@. abs(time_steps - target_dt_save)); dt_save = time_steps[dt_save_idx]; # actual dt_save that will be used

    save_times = range(start = dt_save, step = dt_save, stop = total_time) |> collect
    save_time_idx = [argmin(@. abs(time_steps - save_times[i])) for i in eachindex(save_times)];

    ################ SOLVE FOR TOV BACKGROUND ################
    h_TOV = 1e-4 < h ? 1e-4 : h;
    r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; TD=true, save_to_file=false);

    # Rs = r[end];
    TOV_length = length(m);
    cs = [sqrt(star.dp_dε(ε[i])) for i in 1:TOV_length];
    cs_prime = [cs_prime_func(m[i], p[i], ε[i], cs[i], r[i], star.d2p_dε2) for i in 1:TOV_length];
    cs_prime[1] = 0.0

    # after running TOV at high resolution, downsampled to desired h
    ds_fact = argmin(@. abs(r - h)) - 1; # downsampling factor for saving data
    r = r[1:ds_fact:end];
    m = m[1:ds_fact:end];
    p = p[1:ds_fact:end];
    ε = ε[1:ds_fact:end];
    ν = ν[1:ds_fact:end];
    cs = cs[1:ds_fact:end];
    cs_prime = cs_prime[1:ds_fact:end];
    TOV_length = length(m);

    ################ SET UP SOLUTION ARRAYS ################
    # 2 perturbation variables and 1 independent residual (IR) we want to compute. To compute the IR, we must retain three time levels of the perturbations. 
    ds_fact = argmin(@. abs(r - star.h_save)) - 1; # downsampling factor for saving data
    nPointsSpace_save = 1 + (TOV_length-1) ÷ ds_fact; # number of points to save
    sol_save = zeros(Float64, star.save_every+1, 4, nPointsSpace_save); # 2 perturbation variables, 1 IR (we do save_every + 1 to ensure that at the final save there is always sufficient space to tack on the solution at the final time step)
    spatial_N = TOV_length
    sol_three_levels = zeros(Float64, 3, 2, spatial_N);  # 2 perturbation variables

    # set up arrays for IR computation
    v1_tm1 = zeros(Float64, spatial_N);
    u1_t0 = zeros(Float64, spatial_N);
    v1_t0 = zeros(Float64, spatial_N);
    v1_tp1 = zeros(Float64, spatial_N);

    ################ SET UP INITIAL DATA ################
    nPointsSpace = TOV_length;
    u1_0 = @. star.ξ_ID(r); # ξ
    v1_0 = @. star.ξ_dt_ID(r); # ∂ξ/∂t

    if isfile(fname)
        print_progress ? println("File $fname already exists and will be overwritten.") : nothing
        rm(fname)
    end

    num_saved_steps = length(save_time_idx);
    file = initialize_solution_file(fname, m, p, ε, ν, r, cs, cs_prime, u1_0, v1_0, ds_fact, nPointsSpace_save, num_saved_steps, star.save_every)

    ################ SET UP RK4 STEP ARRAYS ################
    k1 = zeros(Float64, 2, nPointsSpace);
    k2 = zeros(Float64, 2, nPointsSpace);
    k3 = zeros(Float64, 2, nPointsSpace);
    k4 = zeros(Float64, 2, nPointsSpace);

    u1_n = zeros(nPointsSpace);
    v1_n = zeros(nPointsSpace);

    u1_np1 = zeros(nPointsSpace);
    v1_np1 = zeros(nPointsSpace);

    u1_n[1:nPointsSpace] .= u1_0;
    v1_n[1:nPointsSpace] .= v1_0;

    u1_k1 = zero(u1_n);
    v1_k1 = zero(v1_n);

    u1_k2 = zero(u1_n);
    v1_k2 = zero(v1_n);

    u1_k3 = zero(u1_n);
    v1_k3 = zero(v1_n);

    IR = zero(r);

    # time integration loop
    saved_idx_counter = 0;
    save_idx0 = 2; # first index to save to in solution array (we skip the initial data which is already saved)
    saved_time = Float64[0.0];
    for i = 1:total_num_time_steps
        print_string = "Completion: $(round(100 * i/(total_num_time_steps); digits=5))%\r"  
        print_progress ? print(print_string) : nothing

        # Add KO dissipation. Utilize the advanced time step arrays to add the KO dissipation (this is because we do not want to do something like u1_n[j] += -kreiss_oliger(u1_n, j, coef) since we do not want to edit in place)
        u1_np1 .= u1_n;
        v1_np1 .= v1_n;

        @inbounds for j in 3:nPointsSpace-2 # don't apply at all points
            u1_np1[j] += -kreiss_oliger(u1_n, j, star.KO, nPointsSpace)
            v1_np1[j] += -kreiss_oliger(v1_n, j, star.KO, nPointsSpace)
        end

        u1_n .= u1_np1;
        v1_n .= v1_np1;


        compute_ks!(u1_n, v1_n, u1_k1, v1_k1, u1_k2, v1_k2, u1_k3, v1_k3, k1, k2, k3, k4, m, p, ε, ν, cs, cs_prime, r, h, k, nPointsSpace)
 

        for ii in 1:nPointsSpace
            u1_np1[ii] = u1_n[ii] + (k1[1, ii] + 2k2[1, ii] + 2k3[1, ii] + k4[1, ii]) / 6.0
            v1_np1[ii] = v1_n[ii] + (k1[2, ii] + 2k2[2, ii] + 2k3[2, ii] + k4[2, ii]) / 6.0
        end

        # check BCs are still satisfied
        if abs(u1_np1[1]) < 1e-16 && abs(v1_np1[1]) == 0.0
            nothing
        else
            throw(ArgumentError("BCs at r=0 not satisfied! u1_np1[1] = $(u1_np1[1]), v1_np1[1] = $(v1_np1[1]), k1[1,1] = $(k1[1,1]), k2[1,1] = $(k2[1,1]), k3[1,1] = $(k3[1,1]), k4[1,1] = $(k4[1,1])"))
        end
    
        # update three-level solution array
        sol_three_levels[1, :, :] .= sol_three_levels[2, :, :];
        sol_three_levels[2, :, :] .= sol_three_levels[3, :, :];
        sol_three_levels[3, 1, :] .= u1_np1[:];
        sol_three_levels[3, 2, :] .= v1_np1[:];

        # compute IR and send solution to file if at a save time (note that we use i - 1 instead of i because we can only compute the IR with the advanced time step)
        if i - 1 ∈ save_time_idx
            append!(saved_time, time_steps[i - 1])
            compute_IR!(sol_three_levels, IR, m, p, ε, ν, cs, cs_prime, r, v1_tm1, u1_t0, v1_t0, v1_tp1, h, k, nPointsSpace)
            saved_idx_counter += 1;
            update_solution_array!(sol_save, sol_three_levels, IR, saved_idx_counter, ds_fact, nPointsSpace)
            if saved_idx_counter == star.save_every
                HDF5Helper.Cowling.save_to_file_PF!(file, sol_save, save_idx0, save_idx0 + star.save_every - 1)
                saved_idx_counter = 0;
                save_idx0 += star.save_every;
            end
        end

        # initialization for next time step
        u1_n .= u1_np1;
        v1_n .= v1_np1;

        if i == total_num_time_steps
            # do not compute IR at final time step so save solution here manually
            if i ∈ save_time_idx
                saved_idx_counter += 1;
                append!(saved_time, time_steps[save_time_idx[end]])
                IR .= NaN;
                update_solution_array!(sol_save, u1_np1, v1_np1, IR, saved_idx_counter, ds_fact, nPointsSpace)
            end
            HDF5Helper.Cowling.save_final_chunk_to_file_PF!(file, sol_save, save_idx0, save_idx0 + saved_idx_counter - 1, saved_idx_counter, saved_time .* km_to_ms)
        end
    end
end

end

module BDNK
using Printf
using HDF5
using ...HDF5Helper
using ...TOV
using ..CowlingTimeDomain
using ...CowlingTDEqs
using ...FiniteDiffOrder4
using ...BDNKCharacteristicSpeeds
using NeutronStarOscillations

const c = 299792458;
const sec_to_km = c * 1e-3;

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function initialize_solution_file(filename::String, m::AbstractVector{Float64}, p::AbstractVector{Float64}, ε::AbstractVector{Float64}, ν::AbstractVector{Float64}, r::AbstractVector{Float64}, cs::AbstractVector{Float64}, cs_prime::AbstractVector{Float64}, cs_prime_prime::AbstractVector{Float64}, Λ0::AbstractVector{Float64}, Λ1::AbstractVector{Float64}, Λ2::AbstractVector{Float64}, u1::AbstractVector{Float64}, u2::AbstractVector{Float64}, w1::AbstractVector{Float64}, w2::AbstractVector{Float64}, v1::AbstractVector{Float64}, v2::AbstractVector{Float64}, ds::Int64, N::Int, num_time_steps::Int, chunk_size::Int)
    fmode = "w"
    file = h5open(filename, fmode); # filename should include the path

    # create group for solution data
    solution_group_name = "solution"
    HDF5Helper.Cowling.create_file_group!(file, solution_group_name);

    # first store simulation parameters
    # param_group_name = "parameters"
    # HDF5Helper.Cowling.create_file_group!(file, param_group_name);


    # save data fixed at run time
    file["solution"]["r"] = r[1:ds:end];
    file["solution"]["m"] = m[1:ds:end];
    file["solution"]["p"] = p[1:ds:end];
    file["solution"]["ε"] = ε[1:ds:end];
    file["solution"]["ν"] = ν[1:ds:end];
    file["solution"]["cs"] = cs[1:ds:end];
    file["solution"]["cs_prime"] = cs_prime[1:ds:end];
    file["solution"]["cs_prime_prime"] = cs_prime_prime[1:ds:end];
    file["solution"]["Λ0"] = Λ0[1:ds:end];
    file["solution"]["Λ1"] = Λ1[1:ds:end];
    file["solution"]["Λ2"] = Λ2[1:ds:end];

    # create datasets for the solution
    dataspace = ((num_time_steps+1, N), (num_time_steps+1, N))
    chunk = (chunk_size, N)
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "du", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "deps", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "du_dr", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "deps_dr", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "du_dot", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "deps_dot", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "IR1", Float64, dataspace, chunk);
    HDF5Helper.Cowling.create_dataset!(file, solution_group_name, "IR2", Float64, dataspace, chunk);


    # manually save initial data
    file[solution_group_name]["du"][1, :] = u1[1:ds:end];
    file[solution_group_name]["deps"][1, :] = u2[1:ds:end];
    file[solution_group_name]["du_dr"][1, :] = w1[1:ds:end];
    file[solution_group_name]["deps_dr"][1, :] = w2[1:ds:end];
    file[solution_group_name]["du_dot"][1, :] = v1[1:ds:end];
    file[solution_group_name]["deps_dot"][1, :] = v2[1:ds:end];
    file[solution_group_name]["IR1"][1, :] = zero(u1[1:ds:end]) * NaN; # don't compute IR at initial time
    file[solution_group_name]["IR2"][1, :] = zero(u1[1:ds:end]) * NaN;
    return file
end

@views function update_solution_array!(sol::AbstractArray{Float64}, u1::Vector{Float64}, u2::Vector{Float64}, w1::Vector{Float64}, w2::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, IR1::Vector{Float64}, IR2::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = u1[1:ds_fact:nMax];
    sol[time_idx, 2, :] = u2[1:ds_fact:nMax];
    sol[time_idx, 3, :] = w1[1:ds_fact:nMax];
    sol[time_idx, 4, :] = w2[1:ds_fact:nMax];
    sol[time_idx, 5, :] = v1[1:ds_fact:nMax];
    sol[time_idx, 6, :] = v2[1:ds_fact:nMax];
    sol[time_idx, 7, :] = IR1[1:ds_fact:nMax];
    sol[time_idx, 8, :] = IR2[1:ds_fact:nMax];
end

@views function update_solution_array!(sol::AbstractArray{Float64}, sol_three_levels::AbstractArray{Float64}, IR1::Vector{Float64}, IR2::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = sol_three_levels[2, 1, 1:ds_fact:nMax];
    sol[time_idx, 2, :] = sol_three_levels[2, 2, 1:ds_fact:nMax];
    sol[time_idx, 3, :] = sol_three_levels[2, 3, 1:ds_fact:nMax];
    sol[time_idx, 4, :] = sol_three_levels[2, 4, 1:ds_fact:nMax];
    sol[time_idx, 5, :] = sol_three_levels[2, 5, 1:ds_fact:nMax];
    sol[time_idx, 6, :] = sol_three_levels[2, 6, 1:ds_fact:nMax];
    sol[time_idx, 7, :] = IR1[1:ds_fact:nMax];
    sol[time_idx, 8, :] = IR2[1:ds_fact:nMax];
end

function load_solution(filename::String)
    file = h5open(filename, "r")
    finalizer(file) do f
        close(f)  # Automatically close when the file object is garbage collected
    end
    return file 
end

# compute the k1, k2, k3, k4 terms for the RK4 time integration
@views function compute_ks!(
    u1::Vector{Float64}, u2::Vector{Float64}, w1::Vector{Float64}, w2::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, 
    u1_k1::Vector{Float64}, u2_k1::Vector{Float64}, w1_k1::Vector{Float64}, w2_k1::Vector{Float64}, v1_k1::Vector{Float64}, v2_k1::Vector{Float64},
    u1_k2::Vector{Float64}, u2_k2::Vector{Float64}, w1_k2::Vector{Float64}, w2_k2::Vector{Float64}, v1_k2::Vector{Float64}, v2_k2::Vector{Float64}, 
    u1_k3::Vector{Float64}, u2_k3::Vector{Float64}, w1_k3::Vector{Float64}, w2_k3::Vector{Float64}, v1_k3::Vector{Float64}, v2_k3::Vector{Float64},
    k1::Matrix{Float64}, k2::Matrix{Float64}, k3::Matrix{Float64}, k4::Matrix{Float64},
    λ::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64},
    cs_prime_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64, k::Float64, nPoints::Int64)
    
    ## COMPUTE K1 ##
    # use taylor expansion at center
    CenterMax = 1;
    for i = 1:CenterMax
        du1_dt, du2_dt, dw1_dt, dw2_dt, dv1_dt, dv2_dt = CowlingTDEqs.CenterEqs.BDNK.compute_all(λ, cs, p[1], ε[1], ν[1], u1, u2, w1, w2, v1, v2, η, ζ, τε, τP, τQ, L, h)

        k1[1, i] = k * du1_dt
        k1[2, i] = k * du2_dt
        k1[3, i] = k * dw1_dt
        k1[4, i] = k * dw2_dt
        k1[5, i] = k * dv1_dt
        k1[6, i] = k * dv2_dt

        u1_k1[i] = u1[i] + 0.5 * k1[1, i]
        u2_k1[i] = u2[i] + 0.5 * k1[2, i]
        w1_k1[i] = w1[i] + 0.5 * k1[3, i]
        w2_k1[i] = w2[i] + 0.5 * k1[4, i]
        v1_k1[i] = v1[i] + 0.5 * k1[5, i]
        v2_k1[i] = v2[i] + 0.5 * k1[6, i]
    end

    @inbounds for i = (CenterMax+1):(nPoints-1)
        k1[1, i] = k * v1[i]
        k1[2, i] = k * v2[i]
        k1[3, i] = k * (v1[i+1] - v1[i-1])/(2*h)
        k1[4, i] = k * (v2[i+1] - v2[i-1])/(2*h)
        k1[5, i] = k * CowlingTDEqs.RK4.BDNK.dv1_dt(i, u1, u2, w1, w2, v1, v2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
        k1[6, i] = k * CowlingTDEqs.RK4.BDNK.dv2_dt(i, u1, u2, w1, w2, v1, v2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)

        u1_k1[i] = u1[i] + 0.5 * k1[1, i]
        u2_k1[i] = u2[i] + 0.5 * k1[2, i]
        w1_k1[i] = w1[i] + 0.5 * k1[3, i]
        w2_k1[i] = w2[i] + 0.5 * k1[4, i]
        v1_k1[i] = v1[i] + 0.5 * k1[5, i]
        v2_k1[i] = v2[i] + 0.5 * k1[6, i]
    end

    i = nPoints
    k1[1, i] = k * v1[i]
    k1[2, i] = k * v2[i]
    k1[3, i] = k * (v1[-2+i]-4v1[-1+i]+3v1[i])/(2h)
    k1[4, i] = k * (v2[-2+i]-4v2[-1+i]+3v2[i])/(2h)
    k1[5, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv1_dt(i, u1, u2, w1, w2, v1, v2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    k1[6, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv2_dt(i, u1, u2, w1, w2, v1, v2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    

    u1_k1[i] = u1[i] + 0.5 * k1[1, i]
    u2_k1[i] = u2[i] + 0.5 * k1[2, i]
    w1_k1[i] = w1[i] + 0.5 * k1[3, i]
    w2_k1[i] = w2[i] + 0.5 * k1[4, i]
    v1_k1[i] = v1[i] + 0.5 * k1[5, i]
    v2_k1[i] = v2[i] + 0.5 * k1[6, i]


    ## COMPUTE K2 ##
    # use taylor expansion at center
    for i = 1:CenterMax
        du1_dt, du2_dt, dw1_dt, dw2_dt, dv1_dt, dv2_dt = CowlingTDEqs.CenterEqs.BDNK.compute_all(λ, cs, p[1], ε[1], ν[1], u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, η, ζ, τε, τP, τQ, L, h)

        k2[1, i] = k * du1_dt
        k2[2, i] = k * du2_dt
        k2[3, i] = k * dw1_dt
        k2[4, i] = k * dw2_dt
        k2[5, i] = k * dv1_dt
        k2[6, i] = k * dv2_dt

        u1_k2[i] = u1[i] + 0.5 * k2[1, i]
        u2_k2[i] = u2[i] + 0.5 * k2[2, i]
        w1_k2[i] = w1[i] + 0.5 * k2[3, i]
        w2_k2[i] = w2[i] + 0.5 * k2[4, i]
        v1_k2[i] = v1[i] + 0.5 * k2[5, i]
        v2_k2[i] = v2[i] + 0.5 * k2[6, i]
    end

    @inbounds for i = (CenterMax+1):(nPoints-1)
        k2[1, i] = k * v1_k1[i]
        k2[2, i] = k * v2_k1[i]
        k2[3, i] = k * (v1_k1[i+1] - v1_k1[i-1])/(2*h)
        k2[4, i] = k * (v2_k1[i+1] - v2_k1[i-1])/(2*h)
        k2[5, i] = k * CowlingTDEqs.RK4.BDNK.dv1_dt(i, u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
        k2[6, i] = k * CowlingTDEqs.RK4.BDNK.dv2_dt(i, u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)

        u1_k2[i] = u1[i] + 0.5 * k2[1, i]
        u2_k2[i] = u2[i] + 0.5 * k2[2, i]
        w1_k2[i] = w1[i] + 0.5 * k2[3, i]
        w2_k2[i] = w2[i] + 0.5 * k2[4, i]
        v1_k2[i] = v1[i] + 0.5 * k2[5, i]
        v2_k2[i] = v2[i] + 0.5 * k2[6, i]
    end

    i = nPoints
    k2[1, i] = k * v1_k1[i]
    k2[2, i] = k * v2_k1[i]
    k2[3, i] = k * (v1_k1[-2+i]-4v1_k1[-1+i]+3v1_k1[i])/(2h)
    k2[4, i] = k * (v2_k1[-2+i]-4v2_k1[-1+i]+3v2_k1[i])/(2h)
    k2[5, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv1_dt(i, u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    k2[6, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv2_dt(i, u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)

    u1_k2[i] = u1[i] + 0.5 * k2[1, i]
    u2_k2[i] = u2[i] + 0.5 * k2[2, i]
    w1_k2[i] = w1[i] + 0.5 * k2[3, i]
    w2_k2[i] = w2[i] + 0.5 * k2[4, i]
    v1_k2[i] = v1[i] + 0.5 * k2[5, i]
    v2_k2[i] = v2[i] + 0.5 * k2[6, i]


    ## COMPUTE K3 ##
    # use taylor expansion at center
    for i = 1:CenterMax
        du1_dt, du2_dt, dw1_dt, dw2_dt, dv1_dt, dv2_dt = CowlingTDEqs.CenterEqs.BDNK.compute_all(λ, cs, p[1], ε[1], ν[1], u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, η, ζ, τε, τP, τQ, L, h)

        k3[1, i] = k * du1_dt
        k3[2, i] = k * du2_dt
        k3[3, i] = k * dw1_dt
        k3[4, i] = k * dw2_dt
        k3[5, i] = k * dv1_dt
        k3[6, i] = k * dv2_dt

        u1_k3[i] = u1[i] + k3[1, i]
        u2_k3[i] = u2[i] + k3[2, i]
        w1_k3[i] = w1[i] + k3[3, i]
        w2_k3[i] = w2[i] + k3[4, i]
        v1_k3[i] = v1[i] + k3[5, i]
        v2_k3[i] = v2[i] + k3[6, i]
    end

    @inbounds for i = (CenterMax+1):(nPoints-1)
        k3[1, i] = k * v1_k2[i]
        k3[2, i] = k * v2_k2[i]
        k3[3, i] = k * (v1_k2[i+1] - v1_k2[i-1])/(2*h)
        k3[4, i] = k * (v2_k2[i+1] - v2_k2[i-1])/(2*h)
        k3[5, i] = k * CowlingTDEqs.RK4.BDNK.dv1_dt(i, u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
        k3[6, i] = k * CowlingTDEqs.RK4.BDNK.dv2_dt(i, u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)

        u1_k3[i] = u1[i] + k3[1, i]
        u2_k3[i] = u2[i] + k3[2, i]
        w1_k3[i] = w1[i] + k3[3, i]
        w2_k3[i] = w2[i] + k3[4, i]
        v1_k3[i] = v1[i] + k3[5, i]
        v2_k3[i] = v2[i] + k3[6, i]
    end

    i = nPoints
    k3[1, i] = k * v1_k2[i]
    k3[2, i] = k * v2_k2[i]
    k3[3, i] = k * (v1_k2[-2+i]-4v1_k2[-1+i]+3v1_k2[i])/(2h)
    k3[4, i] = k * (v2_k2[-2+i]-4v2_k2[-1+i]+3v2_k2[i])/(2h)
    k3[5, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv1_dt(i, u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    k3[6, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv2_dt(i, u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)

    u1_k3[i] = u1[i] + k3[1, i]
    u2_k3[i] = u2[i] + k3[2, i]
    w1_k3[i] = w1[i] + k3[3, i]
    w2_k3[i] = w2[i] + k3[4, i]
    v1_k3[i] = v1[i] + k3[5, i]
    v2_k3[i] = v2[i] + k3[6, i]

    ## COMPUTE K4 ##
    # use taylor expansion at center
    for i = 1:CenterMax
        du1_dt, du2_dt, dw1_dt, dw2_dt, dv1_dt, dv2_dt = CowlingTDEqs.CenterEqs.BDNK.compute_all(λ, cs, p[1], ε[1], ν[1], u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3, η, ζ, τε, τP, τQ, L, h)

        k4[1, i] = k * du1_dt
        k4[2, i] = k * du2_dt
        k4[3, i] = k * dw1_dt
        k4[4, i] = k * dw2_dt
        k4[5, i] = k * dv1_dt
        k4[6, i] = k * dv2_dt
    end

    @inbounds for i = (CenterMax+1):(nPoints-1)
        k4[1, i] = k * v1_k3[i]
        k4[2, i] = k * v2_k3[i]
        k4[3, i] = k * (v1_k3[i+1] - v1_k3[i-1])/(2*h)
        k4[4, i] = k * (v2_k3[i+1] - v2_k3[i-1])/(2*h)
        k4[5, i] = k * CowlingTDEqs.RK4.BDNK.dv1_dt(i, u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
        k4[6, i] = k * CowlingTDEqs.RK4.BDNK.dv2_dt(i, u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    end

    i = nPoints
    k4[1, i] = k * v1_k3[i]
    k4[2, i] = k * v2_k3[i]
    k4[3, i] = k * (v1_k3[-2+i]-4v1_k3[-1+i]+3v1_k3[i])/(2h)
    k4[4, i] = k * (v2_k3[-2+i]-4v2_k3[-1+i]+3v2_k3[i])/(2h)
    k4[5, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv1_dt(i, u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h)
    k4[6, i] = k * CowlingTDEqs.RK4Surface.BDNK.dv2_dt(i, u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h) 
end


@views function compute_IR!(sol_three_levels::AbstractArray, IR1::Vector{Float64}, IR2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, cs_prime_prime::Vector{Float64}, r::Vector{Float64}, u1_tm1::Vector{Float64}, u2_tm1::Vector{Float64}, w1_tm1::Vector{Float64}, w2_tm1::Vector{Float64}, v1_tm1::Vector{Float64}, v2_tm1::Vector{Float64}, u1_t0::Vector{Float64}, u2_t0::Vector{Float64}, w1_t0::Vector{Float64}, w2_t0::Vector{Float64}, v1_t0::Vector{Float64}, v2_t0::Vector{Float64}, u1_tp1::Vector{Float64}, u2_tp1::Vector{Float64}, w1_tp1::Vector{Float64}, w2_tp1::Vector{Float64}, v1_tp1::Vector{Float64}, v2_tp1::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64, k::Float64, nPoints::Int64)
    # fill arrays
    u1_tm1 .= sol_three_levels[1, 1, :];
    u2_tm1 .= sol_three_levels[1, 2, :];
    w1_tm1 .= sol_three_levels[1, 3, :];
    w2_tm1 .= sol_three_levels[1, 4, :];
    v1_tm1 .= sol_three_levels[1, 5, :];
    v2_tm1 .= sol_three_levels[1, 6, :];

    u1_t0 .= sol_three_levels[2, 1, :];
    u2_t0 .= sol_three_levels[2, 2, :];
    w1_t0 .= sol_three_levels[2, 3, :];
    w2_t0 .= sol_three_levels[2, 4, :];
    v1_t0 .= sol_three_levels[2, 5, :];
    v2_t0 .= sol_three_levels[2, 6, :];

    u1_tp1 .= sol_three_levels[3, 1, :];
    u2_tp1 .= sol_three_levels[3, 2, :];
    w1_tp1 .= sol_three_levels[3, 3, :];
    w2_tp1 .= sol_three_levels[3, 4, :];
    v1_tp1 .= sol_three_levels[3, 5, :];
    v2_tp1 .= sol_three_levels[3, 6, :];

    IR1[1] = 0.0
    IR2[1] = 0.0

    for i = 2:nPoints-1
        IR1[i] = CowlingTDEqs.Leapfrog.BDNK.Eq1(i, v1_tm1, v2_tm1, u1_t0, u2_t0, w1_t0, w2_t0, v1_t0, v2_t0, v1_tp1, v2_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
        IR2[i] = CowlingTDEqs.Leapfrog.BDNK.Eq2(i, v1_tm1, v2_tm1, u1_t0, u2_t0, w1_t0, w2_t0, v1_t0, v2_t0, v1_tp1, v2_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
    end
    i = nPoints
    IR1[i] = CowlingTDEqs.LeapfrogSurface.BDNK.Eq1(i, v1_tm1, v2_tm1, u1_t0, u2_t0, w1_t0, w2_t0, v1_t0, v2_t0, v1_tp1, v2_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
    IR2[i] = CowlingTDEqs.LeapfrogSurface.BDNK.Eq2(i, v1_tm1, v2_tm1, u1_t0, u2_t0, w1_t0, w2_t0, v1_t0, v2_t0, v1_tp1, v2_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
end

function kreiss_oliger(var::AbstractVector{Float64}, j::Int, coef::Float64, nPoints::Int64)
    if j == 1
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[1]
        u_plus_1 = var[2]
        u_plus_2 = var[3]
    elseif j == 2
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[2]
        u_plus_1 = var[3]
        u_plus_2 = var[4]
    elseif j == nPoints-1
        u_minus_2 = var[nPoints-3]
        u_minus_1 = var[nPoints-2]
        u = var[nPoints-1]
        u_plus_1 = var[nPoints]
        u_plus_2 = 0.0
    elseif j == nPoints
        u_minus_2 = var[nPoints-2]
        u_minus_1 = var[nPoints-1]
        u = var[nPoints]
        u_plus_1 = 0.0
        u_plus_2 = 0.0
    else
        u_minus_2 = var[j-2]
        u_minus_1 = var[j-1]
        u = var[j]
        u_plus_1 = var[j+1]
        u_plus_2 = var[j+2]
    end

    return coef * (u_plus_2 - 4.0 * u_plus_1 + 6.0 * u - 4.0 * u_minus_1 + u_minus_2) / 16.0
end

function solve(star::NeutronStarOscillations.Star, h::Float64; print_progress::Bool=true)
    fname = CowlingTimeDomain.td_fname(star, h)
    if abs(star.δu_ID(0.0)) > 1e-16
        error("|δu(r=0)| > 1e-16. δu must be zero at the center by regularity. δu(0.0) = $(star.δu_ID(0.0))")
    end

    # convert time from ms to km
    km_to_ms = 1e3 / (sec_to_km)
    target_dt_save = star.dt_save / km_to_ms
    total_time = star.T / km_to_ms

    k = star.CFL * h; # time step from CFL condition
    total_num_time_steps = Int(ceil(total_time / k));
    time_steps = range(start = k, step = k, stop = total_time) |> collect
    dt_save_idx = argmin(@. abs(time_steps - target_dt_save)); dt_save = time_steps[dt_save_idx]; # actual dt_save that will be used

    save_times = range(start = dt_save, step = dt_save, stop = total_time) |> collect
    save_time_idx = [argmin(@. abs(time_steps - save_times[i])) for i in eachindex(save_times)];

    ################ TOV BACKGROUND ################
    h_TOV = 1e-4;
    r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; TD=true, save_to_file=false);

    TOV_length = length(m);
    cs = [sqrt(star.dp_dε(ε[i])) for i in 1:TOV_length];
    cs_prime = [cs_prime_func(m[i], p[i], ε[i], cs[i], r[i], star.d2p_dε2) for i in 1:TOV_length];
    cs_prime[1] = 0.0
    λ = @. -log(1 - 2 * m / r);
    λ[1] = 0.0;

    cs_prime_prime = zero(cs_prime)
    FiniteDiffOrder4.compute_first_derivative(cs_prime_prime, cs_prime, diff(r)[1], length(r));

    ################ INITIAL DATA ################
    u1_0 = @. star.δu_ID(r);
    u2_0 = @. star.δε_ID(r);
    w1_0 = @. star.δu_dr_ID(r);
    w2_0 = @. star.δε_dr_ID(r);
    v1_0 = @. star.δu_dt_ID(r);
    v2_0 = @. star.δε_dt_ID(r);

    # downsample
    ds_fact = argmin(@. abs(r - h)) - 1;
    r = r[1:ds_fact:end];
    diff(r)[1] ≈ h ? nothing : error("Grid spacing does not match desired value of h after downsampling");
    m = m[1:ds_fact:end];
    p = p[1:ds_fact:end];
    ε = ε[1:ds_fact:end];
    ν = ν[1:ds_fact:end];
    λ = λ[1:ds_fact:end];
    cs = cs[1:ds_fact:end];
    cs_prime = cs_prime[1:ds_fact:end];
    cs_prime_prime = cs_prime_prime[1:ds_fact:end];
    u1_0 = u1_0[1:ds_fact:end];
    u2_0 = u2_0[1:ds_fact:end];
    w1_0 = w1_0[1:ds_fact:end];
    w2_0 = w2_0[1:ds_fact:end];
    v1_0 = v1_0[1:ds_fact:end];
    v2_0 = v2_0[1:ds_fact:end];
    TOV_length = length(m);

    # arrays for characteristic speed computation
    Λ0 = zeros(Float64, TOV_length);
    Λ1 = zeros(Float64, TOV_length);
    Λ2 = zeros(Float64, TOV_length);

    for i in eachindex(Λ0)
        ΛΛ0, ΛΛ1, ΛΛ2 = BDNKCharacteristicSpeeds.compute(p[i], ε[i], cs[i], star.η, star.ζ, star.τε, star.τP, star.τQ, star.L)
        Λ0[i] = ΛΛ0
        Λ1[i] = ΛΛ1
        Λ2[i] = ΛΛ2
    end

    ################ SET UP SOLUTION ARRAYS ################
    # 6 perturbation variables and 2 independent residuals (IR) we want to compute. To compute the IR, we must retain three time levels of the perturbations. 
    ds_fact = argmin(@. abs(r - star.h_save)) - 1; # downsampling factor for saving data
    nPointsSpace_save = 1 + (TOV_length-1) ÷ ds_fact; # number of points to save
    sol_save = zeros(Float64, star.save_every+1, 8, nPointsSpace_save); # 6 perturbation variables, 2 IRs (we do save_every + 1 to ensure that at the final save there is always sufficient space to tack on the solution at the final time step)
    spatial_N = TOV_length
    sol_three_levels = zeros(Float64, 3, 6, spatial_N); # 6 perturbation variables

    # set up arrays for IR computation
    u1_tm1 = zeros(Float64, spatial_N);
    u2_tm1 = zeros(Float64, spatial_N);
    w1_tm1 = zeros(Float64, spatial_N);
    w2_tm1 = zeros(Float64, spatial_N);
    v1_tm1 = zeros(Float64, spatial_N);
    v2_tm1 = zeros(Float64, spatial_N);

    u1_t0 = zeros(Float64, spatial_N);
    u2_t0 = zeros(Float64, spatial_N);
    w1_t0 = zeros(Float64, spatial_N);
    w2_t0 = zeros(Float64, spatial_N);
    v1_t0 = zeros(Float64, spatial_N);
    v2_t0 = zeros(Float64, spatial_N);

    u1_tp1 = zeros(Float64, spatial_N);
    u2_tp1 = zeros(Float64, spatial_N);
    w1_tp1 = zeros(Float64, spatial_N);
    w2_tp1 = zeros(Float64, spatial_N);
    v1_tp1 = zeros(Float64, spatial_N);
    v2_tp1 = zeros(Float64, spatial_N);

    ################ SET UP INITIAL DATA ################
    nPointsSpace = TOV_length;

    if isfile(fname)
        print_progress ? println("File $fname already exists and will be overwritten.") : nothing
        rm(fname)
    end

    num_saved_steps = length(save_time_idx);
    file = initialize_solution_file(fname, m, p, ε, ν, r, cs, cs_prime, cs_prime_prime, Λ0, Λ1, Λ2, u1_0, u2_0, w1_0, w2_0, v1_0, v2_0, ds_fact, nPointsSpace_save, num_saved_steps, star.save_every)

    ################ SET UP RK4 STEP ARRAYS ################
    k1 = zeros(Float64, 6, nPointsSpace);
    k2 = zeros(Float64, 6, nPointsSpace);
    k3 = zeros(Float64, 6, nPointsSpace);
    k4 = zeros(Float64, 6, nPointsSpace);

    u1_n = zeros(nPointsSpace);
    u2_n = zeros(nPointsSpace);
    w1_n = zeros(nPointsSpace);
    w2_n = zeros(nPointsSpace);
    v1_n = zeros(nPointsSpace);
    v2_n = zeros(nPointsSpace);

    u1_np1 = zeros(nPointsSpace);
    u2_np1 = zeros(nPointsSpace);
    w1_np1 = zeros(nPointsSpace);
    w2_np1 = zeros(nPointsSpace);
    v1_np1 = zeros(nPointsSpace);
    v2_np1 = zeros(nPointsSpace);

    u1_n[1:nPointsSpace] .= u1_0;
    u2_n[1:nPointsSpace] .= u2_0;
    w1_n[1:nPointsSpace] .= w1_0;
    w2_n[1:nPointsSpace] .= w2_0;
    v1_n[1:nPointsSpace] .= v1_0;
    v2_n[1:nPointsSpace] .= v2_0;

    u1_k1 = zero(u1_n);
    u2_k1 = zero(u2_n);
    w1_k1 = zero(w1_n);
    w2_k1 = zero(w2_n);
    v1_k1 = zero(v1_n);
    v2_k1 = zero(v2_n);

    u1_k2 = zero(u1_n);
    u2_k2 = zero(u2_n);
    w1_k2 = zero(w1_n);
    w2_k2 = zero(w2_n);
    v1_k2 = zero(v1_n);
    v2_k2 = zero(v2_n);

    u1_k3 = zero(u1_n);
    u2_k3 = zero(u2_n);
    w1_k3 = zero(w1_n);
    w2_k3 = zero(w2_n);
    v1_k3 = zero(v1_n);
    v2_k3 = zero(v2_n);

    IR1 = zero(r);
    IR2 = zero(r);

    # time integration loop
    saved_idx_counter = 0;
    save_idx0 = 2; # first index to save to in solution array (we skip the initial data which is already saved)
    saved_time = Float64[0.0];
    for i = 1:total_num_time_steps
        print_string = "Completion: $(round(100 * i/(total_num_time_steps); digits=5))%\r"  
        print_progress ? print(print_string) : nothing

        # Add KO dissipation. Utilize the advanced time step arrays to add the KO dissipation (this is because we do not want to do something like u1_n[j] += -kreiss_oliger(u1_n, j, coef) since we do not want to edit in place)
        u1_np1 .= u1_n;
        u2_np1 .= u2_n;
        w1_np1 .= w1_n;
        w2_np1 .= w2_n;
        v1_np1 .= v1_n;
        v2_np1 .= v2_n;
        @inbounds for j in 3:nPointsSpace-2 # don't apply at all points
            u1_np1[j] += -kreiss_oliger(u1_n, j, star.KO, nPointsSpace)
            u2_np1[j] += -kreiss_oliger(u2_n, j, star.KO, nPointsSpace)
            w1_np1[j] += -kreiss_oliger(w1_n, j, star.KO, nPointsSpace)
            w2_np1[j] += -kreiss_oliger(w2_n, j, star.KO, nPointsSpace)
            v1_np1[j] += -kreiss_oliger(v1_n, j, star.KO, nPointsSpace)
            v2_np1[j] += -kreiss_oliger(v2_n, j, star.KO, nPointsSpace)
        end

        u1_n .= u1_np1;
        u2_n .= u2_np1;
        w1_n .= w1_np1;
        w2_n .= w2_np1;
        v1_n .= v1_np1;
        v2_n .= v2_np1;

        compute_ks!(
            u1_n, u2_n, w1_n, w2_n, v1_n, v2_n,
            u1_k1, u2_k1, w1_k1, w2_k1, v1_k1, v2_k1, 
            u1_k2, u2_k2, w1_k2, w2_k2, v1_k2, v2_k2, 
            u1_k3, u2_k3, w1_k3, w2_k3, v1_k3, v2_k3,
            k1, k2, k3, k4,
            λ, m, p, ε, ν, cs, cs_prime, cs_prime_prime, r,
            star.η, star.ζ, star.τε, star.τP, star.τQ,
            star.L, h, k, nPointsSpace)

        for ii in 1:nPointsSpace
            u1_np1[ii] = u1_n[ii] + (k1[1, ii] + 2k2[1, ii] + 2k3[1, ii] + k4[1, ii]) / 6.0
            u2_np1[ii] = u2_n[ii] + (k1[2, ii] + 2k2[2, ii] + 2k3[2, ii] + k4[2, ii]) / 6.0
            w1_np1[ii] = w1_n[ii] + (k1[3, ii] + 2k2[3, ii] + 2k3[3, ii] + k4[3, ii]) / 6.0
            w2_np1[ii] = w2_n[ii] + (k1[4, ii] + 2k2[4, ii] + 2k3[4, ii] + k4[4, ii]) / 6.0
            v1_np1[ii] = v1_n[ii] + (k1[5, ii] + 2k2[5, ii] + 2k3[5, ii] + k4[5, ii]) / 6.0
            v2_np1[ii] = v2_n[ii] + (k1[6, ii] + 2k2[6, ii] + 2k3[6, ii] + k4[6, ii]) / 6.0
        end
    
        # update three-level solution array
        sol_three_levels[1, :, :] .= sol_three_levels[2, :, :];
        sol_three_levels[2, :, :] .= sol_three_levels[3, :, :];
        sol_three_levels[3, 1, :] .= u1_np1[:];
        sol_three_levels[3, 2, :] .= u2_np1[:];
        sol_three_levels[3, 3, :] .= w1_np1[:];
        sol_three_levels[3, 4, :] .= w2_np1[:];
        sol_three_levels[3, 5, :] .= v1_np1[:];
        sol_three_levels[3, 6, :] .= v2_np1[:];


        # compute IR and send solution to file if at a save time (note that we use i - 1 instead of i because we can only compute the IR with the advanced time step)
        if i - 1 ∈ save_time_idx
            append!(saved_time, time_steps[i - 1])
            compute_IR!(sol_three_levels, IR1, IR2, m, p, ε, ν, cs, cs_prime, cs_prime_prime, r, u1_tm1, u2_tm1, w1_tm1, w2_tm1, v1_tm1, v2_tm1, u1_t0, u2_t0, w1_t0, w2_t0, v1_t0, v2_t0, u1_tp1, u2_tp1, w1_tp1, w2_tp1, v1_tp1, v2_tp1, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, h, k, nPointsSpace)
            saved_idx_counter += 1;
            update_solution_array!(sol_save, sol_three_levels, IR1, IR2, saved_idx_counter, ds_fact, nPointsSpace)
            if saved_idx_counter == star.save_every
                HDF5Helper.Cowling.save_to_file_BDNK!(file, sol_save, save_idx0, save_idx0 + star.save_every - 1)
                saved_idx_counter = 0;
                save_idx0 += star.save_every;
            end
        end

        # initialization for next time step
        u1_n .= u1_np1;
        u2_n .= u2_np1;
        w1_n .= w1_np1;
        w2_n .= w2_np1;
        v1_n .= v1_np1;
        v2_n .= v2_np1;

        if i == total_num_time_steps
            # do not compute IR at final time step so save solution here manually
            if i ∈ save_time_idx
                saved_idx_counter += 1;
                append!(saved_time, time_steps[save_time_idx[end]])
                IR1 .= NaN;
                IR2 .= NaN;
                update_solution_array!(sol_save, u1_np1, u2_np1, w1_np1, w2_np1, v1_np1, v2_np1, IR1, IR2, saved_idx_counter, ds_fact, nPointsSpace)
            end
            HDF5Helper.Cowling.save_final_chunk_to_file_BDNK!(file, sol_save, save_idx0, save_idx0 + saved_idx_counter - 1, saved_idx_counter, saved_time .* km_to_ms)
        end
    end
end

end

end
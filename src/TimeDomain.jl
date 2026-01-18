#=

    Module comprising of functions which solve the perfect fluid, Eckart and BDNK equations of motion governing linear radial perturbations of polytropic neutron stars in the
    time domain. (Note we don't provide separate perfect fluid functions; one takes η = ζ = 0 in the Eckart functions). The Eckart equations of motion are integrated in time
    using RK4 while the BDNK system is integrated using an implicit Cranck-Nicholson discretization. Both numerical schemes are second-second order accurate since we use
    second-order finite differences for spatial derivatives. Functions to plot the results and compute convergence factors are also included.

=#

module TimeDomain
using LaTeXStrings
using CairoMakie
using ..PlotSettings
using ..QuickPlots
using ..Animations
using NeutronStarOscillations
using Printf
using HDF5

one_norm(x::AbstractArray) = sum(abs, x) / length(x)
td_fname(star::NeutronStarOscillations.Star, h::Float64)::String = star.data_path * @sprintf("Time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h);
plot_fname(star::NeutronStarOscillations.Star, h::Float64, var::String, time::Float64)::String = star.fig_path * @sprintf("Time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e_%s_time_%.2f.png", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h, var, time);
animation_fname(star::NeutronStarOscillations.Star, h::Float64, var::String)::String = star.fig_path * @sprintf("Time_domain_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s_h_%.1e_%s.mp4", star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T, h, var);
convergence_plot_fname(star::NeutronStarOscillations.Star, var::String)::String = star.fig_path * @sprintf("Time_domain_%s_convergence_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_η_%.1e_ζ_%.1e_τε_%.1e_τP_%.1e_τQ_%.1e_L_%s_KO_%s_CFL_%s_T_%s.png", var, star.kappa, star.n, star.εc_SI, star.ptol_TD, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, star.KO, star.CFL, star.T);

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
const static_vars::Vector{String} = [ "cs", "cs_prime", "m", "p", "ε", "ν"];
const eckart_IRs::Vector{String} = ["IR"]; # independent residuals for Eckart stars
const BDNK_IRs::Vector{String} = ["Constrained_IR1", "Constained_IR2", "Wave_IR1", "Wave_IR2", "Wave_IR3"]; # independent residuals for BDNK stars

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
        linestyles, alphas, linewidths, labels; legend = legend, position = position, fix_ylims = fix_ylims, ymin = lim_y_min, ymax = lim_y_max)
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

function plot_initial_data_convergence(star::NeutronStarOscillations.Star, h::Float64)
    if star.τε == star.τP == star.τQ == 0.0
        throw(ArgumentError("Characteristic speeds can only be computed for BDNK stars!"))
    end
    
    sol = load_solution(star, h);
    

    x = [sol["solution/r"][:], sol["solution/r"][:]];
    y = [sol["solution/Q_u2"][:], sol["solution/Q_u3"][:]];
    labels = [L"\delta{\epsilon}", L"\delta\lambda"]
    ylabel = L"Q_{N}";
    xlabel = L"r\,[\mathrm{km}]";

    close(sol)

    lim_y_min = 0.0; lim_y_max = 16.0;
    NeutronStarOscillations.QuickPlots.plot11(
        x, y;
        labels = labels,
        xlabel = xlabel,
        ylabel = ylabel,
        legend = true,
        position = :rt,
        framevisible = true,
        lim_y_min = lim_y_min,
        lim_y_max = lim_y_max,
        hlines = [4.0],
        )
end

module Eckart
using Printf
using HDF5
using ..TimeDomain
using ...HDF5Helper
using ...TOV
using NeutronStarOscillations
const c = 299792458;
const sec_to_km = c * 1e-3;

# non-trivial time integration equation. At surface use backwards differences
RK4Eq(i::Int64, ξ::Vector{Float64}, Ξ::Vector{Float64}, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, h::Float64)::Float64 = -0.08333333333333333*(3*exp(ν)*r*(r - 2*m)*(4*r*(r - 2*m)*(2*ξ[i] - ξ[i-1] - ξ[i+1])*cs^2 - h*(ξ[i-1] - ξ[i+1])*(2*m + 8*π*r^3*p - 2*(2*r - 5*m + 4*π*r^3*ε)*cs^2 - 4*r*(r - 2*m)*cs*cs_prime)) - 6*exp(ν)*h^2*ξ[i]*(r^2 - (r - 2*m)^2 + 64*π^2*r^6*p^2 + (r^2 - 5*(r - 2*m)^2 + 8*π*r^3*(r - 4*m)*ε)*cs^2 + 4*r*(2*r - 5*m)*(r - 2*m)*cs*cs_prime + 8*π*r^3*p*(4*m + r*(1 - 8*π*r^2*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime)) + 4*exp(ν/2.)*η*(-(h^2*Ξ[i]*(r^2 - (r - 2*m)^2 - 4*(5*m^2 - 8*m*(r + π*r^3*ε) + 2*(r^2 + 4*π*r^4*ε))*cs^2 + 64*π^2*r^6*p^2*(1 + cs^2) - 8*r*(r - 2*m)*(r - m)*cs*cs_prime + 16*π*r^3*p*(r - 4*(r - 3*m + 2*π*r^3*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime))) + r*(r - 2*m)*(4*r*(r - 2*m)*(2*Ξ[i] - Ξ[i-1] - Ξ[i+1])*cs^2 - 2*h*(Ξ[i-1] - Ξ[i+1])*(m - 2*(r - 3*m + 2*π*r^3*ε)*cs^2 + 4*π*r^3*p*(1 + cs^2) - 2*r*(r - 2*m)*cs*cs_prime))) + 3*exp(ν/2.)*ζ*(-(h^2*Ξ[i]*(r^2 - 6*r*(r - 2*m) + 5*(r - 2*m)^2 + 4*(-2*r^2 - 5*m^2 + 4*π*r^4*ε + m*(8*r - 16*π*r^3*ε))*cs^2 + 64*π^2*r^6*p^2*(1 + cs^2) + 8*r*(2*r - 5*m)*(r - 2*m)*cs*cs_prime + 16*π*r^3*p*(-2*r + 6*m - (r - 6*m + 8*π*r^3*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime))) + r*(r - 2*m)*(4*r*(r - 2*m)*(2*Ξ[i] - Ξ[i-1] - Ξ[i+1])*cs^2 - 2*h*(Ξ[i-1] - Ξ[i+1])*(m - 2*(r - 3*m + 2*π*r^3*ε)*cs^2 + 4*π*r^3*p*(1 + cs^2) - 2*r*(r - 2*m)*cs*cs_prime))))/(h^2*r^3*(r - 2*m))
RK4EqSurface(i::Int64, ξ::Vector{Float64}, Ξ::Vector{Float64}, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, h::Float64)::Float64 = (3*exp(ν)*r*(r - 2*m)*(4*r*(r - 2*m)*(2*ξ[i] - ξ[i-3] + 4*ξ[i-2] - 5*ξ[i-1])*cs^2 - h*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*(2*m + 8*π*r^3*p - 2*(2*r - 5*m + 4*π*r^3*ε)*cs^2 - 4*r*(r - 2*m)*cs*cs_prime)) + 6*exp(ν)*h^2*ξ[i]*(r^2 - (r - 2*m)^2 + 64*π^2*r^6*p^2 + (r^2 - 5*(r - 2*m)^2 + 8*π*r^3*(r - 4*m)*ε)*cs^2 + 4*r*(2*r - 5*m)*(r - 2*m)*cs*cs_prime + 8*π*r^3*p*(4*m + r*(1 - 8*π*r^2*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime)) + 4*exp(ν/2.)*η*(h^2*Ξ[i]*(r^2 - (r - 2*m)^2 - 4*(5*m^2 - 8*m*(r + π*r^3*ε) + 2*(r^2 + 4*π*r^4*ε))*cs^2 + 64*π^2*r^6*p^2*(1 + cs^2) - 8*r*(r - 2*m)*(r - m)*cs*cs_prime + 16*π*r^3*p*(r - 4*(r - 3*m + 2*π*r^3*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime)) + r*(r - 2*m)*(4*r*(r - 2*m)*(2*Ξ[i] - Ξ[i-3] + 4*Ξ[i-2] - 5*Ξ[i-1])*cs^2 - 2*h*(3*Ξ[i] + Ξ[i-2] - 4*Ξ[i-1])*(m - 2*(r - 3*m + 2*π*r^3*ε)*cs^2 + 4*π*r^3*p*(1 + cs^2) - 2*r*(r - 2*m)*cs*cs_prime))) + 3*exp(ν/2.)*ζ*(h^2*Ξ[i]*(r^2 - 6*r*(r - 2*m) + 5*(r - 2*m)^2 + 4*(-2*r^2 - 5*m^2 + 4*π*r^4*ε + m*(8*r - 16*π*r^3*ε))*cs^2 + 64*π^2*r^6*p^2*(1 + cs^2) + 8*r*(2*r - 5*m)*(r - 2*m)*cs*cs_prime + 16*π*r^3*p*(-2*r + 6*m - (r - 6*m + 8*π*r^3*ε)*cs^2 - 2*r*(r - 2*m)*cs*cs_prime)) + r*(r - 2*m)*(4*r*(r - 2*m)*(2*Ξ[i] - Ξ[i-3] + 4*Ξ[i-2] - 5*Ξ[i-1])*cs^2 - 2*h*(3*Ξ[i] + Ξ[i-2] - 4*Ξ[i-1])*(m - 2*(r - 3*m + 2*π*r^3*ε)*cs^2 + 4*π*r^3*p*(1 + cs^2) - 2*r*(r - 2*m)*cs*cs_prime))))/(12.0*h^2*r^3*(r - 2*m))

# leapfrog discretization of the time integration equation used to monitor convergence
LeapFrogEq(i::Int64, ξ::Vector{Float64}, Ξ_tm1::Vector{Float64}, Ξ_t0::Vector{Float64}, Ξ_tp1::Vector{Float64}, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, h::Float64, k::Float64)::Float64 =((r - 2*m)*(-6*exp(ν)*p*ξ[i] - 6*exp(ν)*ε*ξ[i] + (3*exp(ν)*r*p*(-ξ[i-1] + ξ[i+1]))/h + (3*exp(ν)*r*ε*(-ξ[i-1] + ξ[i+1]))/h - 30*exp(ν)*p*ξ[i]*cs^2 - 30*exp(ν)*ε*ξ[i]*cs^2 + (15*exp(ν)*r*p*(-ξ[i-1] + ξ[i+1])*cs^2)/h + (15*exp(ν)*r*ε*(-ξ[i-1] + ξ[i+1])*cs^2)/h + (12*exp(ν)*r^2*p*(-2*ξ[i] + ξ[i-1] + ξ[i+1])*cs^2)/h^2 + (12*exp(ν)*r^2*ε*(-2*ξ[i] + ξ[i-1] + ξ[i+1])*cs^2)/h^2 - 30*exp(ν/2.)*ζ*(p + ε)*Ξ_t0[i]*cs^2 - 16*exp(ν/2.)*η*(p + ε)*Ξ_t0[i]*cs^2 + (15*exp(ν/2.)*r*ζ*(p + ε)*(-Ξ_t0[i-1] + Ξ_t0[i+1])*cs^2)/h + (20*exp(ν/2.)*r*η*(p + ε)*(-Ξ_t0[i-1] + Ξ_t0[i+1])*cs^2)/h + (12*exp(ν/2.)*r^2*ζ*(p + ε)*(-2*Ξ_t0[i] + Ξ_t0[i-1] + Ξ_t0[i+1])*cs^2)/h^2 + (16*exp(ν/2.)*r^2*η*(p + ε)*(-2*Ξ_t0[i] + Ξ_t0[i-1] + Ξ_t0[i+1])*cs^2)/h^2 + (2*exp(ν/2.)*r^2*(1 + 8*π*r^2*p)*(p + ε)*(-((3*ζ + 4*η)*(-1 + 8*π*r^2*ε)*Ξ_t0[i]*cs^2) + 3*exp(ν/2.)*ξ[i]*(1 + 8*π*r^2*p + (1 - 8*π*r^2*ε)*cs^2)))/(r - 2*m)^2 + 60*exp(ν)*r*p*ξ[i]*cs*cs_prime + 60*exp(ν)*r*ε*ξ[i]*cs*cs_prime + (12*exp(ν)*r^2*p*(-ξ[i-1] + ξ[i+1])*cs*cs_prime)/h + (12*exp(ν)*r^2*ε*(-ξ[i-1] + ξ[i+1])*cs*cs_prime)/h - (30*exp(ν/2.)*ζ*(p + ε)*Ξ_t0[i]*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(r - 2*m) + (8*exp(ν/2.)*η*(p + ε)*Ξ_t0[i]*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(r - 2*m) + (6*exp(ν/2.)*r*ζ*(p + ε)*(Ξ_t0[i-1] - Ξ_t0[i+1])*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(h*(r - 2*m)) + (8*exp(ν/2.)*r*η*(p + ε)*(Ξ_t0[i-1] - Ξ_t0[i+1])*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(h*(r - 2*m)) + (r*(p + ε)*((6*r^2*(Ξ_tm1[i] - Ξ_tp1[i]))/k + (3*exp(ν)*r*((ξ[i-1] - ξ[i+1])*(1 + 8*π*r^2*p + (1 - 8*π*r^2*ε)*cs^2) - 4*h*ξ[i]*(8*π*r*(p - ε*cs^2) + (1 + 8*π*r^2*p)*cs*cs_prime)))/h + exp(ν/2.)*(-((r*(3*ζ + 4*η)*(-1 + 8*π*r^2*ε)*(Ξ_t0[i-1] - Ξ_t0[i+1])*cs^2)/h) + 2*Ξ_t0[i]*(-12*η*(1 + 16*π*r^2*p)*cs^2 + 16*π*r^2*(3*ζ - 2*η)*ε*cs^2 + ((3*ζ + 4*η)*(1 + 8*π*r^2*p)*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(r - 2*m)))))/(r - 2*m)))/(12.0*r^3*(p + ε))
LeapFrogEqSurface(i::Int64, ξ::Vector{Float64}, Ξ_tm1::Vector{Float64}, Ξ_t0::Vector{Float64}, Ξ_tp1::Vector{Float64}, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, h::Float64, k::Float64)::Float64 = ((r - 2*m)*(-6*exp(ν)*p*ξ[i] - 6*exp(ν)*ε*ξ[i] + (3*exp(ν)*r*p*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1]))/h + (3*exp(ν)*r*ε*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1]))/h - 30*exp(ν)*p*ξ[i]*cs^2 - 30*exp(ν)*ε*ξ[i]*cs^2 + (12*exp(ν)*r^2*p*(2*ξ[i] - ξ[i-3] + 4*ξ[i-2] - 5*ξ[i-1])*cs^2)/h^2 + (12*exp(ν)*r^2*ε*(2*ξ[i] - ξ[i-3] + 4*ξ[i-2] - 5*ξ[i-1])*cs^2)/h^2 + (15*exp(ν)*r*p*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*cs^2)/h + (15*exp(ν)*r*ε*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*cs^2)/h - 30*exp(ν/2.)*ζ*(p + ε)*Ξ_t0[i]*cs^2 - 16*exp(ν/2.)*η*(p + ε)*Ξ_t0[i]*cs^2 + (12*exp(ν/2.)*r^2*ζ*(p + ε)*(2*Ξ_t0[i] - Ξ_t0[i-3] + 4*Ξ_t0[i-2] - 5*Ξ_t0[i-1])*cs^2)/h^2 + (16*exp(ν/2.)*r^2*η*(p + ε)*(2*Ξ_t0[i] - Ξ_t0[i-3] + 4*Ξ_t0[i-2] - 5*Ξ_t0[i-1])*cs^2)/h^2 + (15*exp(ν/2.)*r*ζ*(p + ε)*(3*Ξ_t0[i] + Ξ_t0[i-2] - 4*Ξ_t0[i-1])*cs^2)/h + (20*exp(ν/2.)*r*η*(p + ε)*(3*Ξ_t0[i] + Ξ_t0[i-2] - 4*Ξ_t0[i-1])*cs^2)/h + (2*exp(ν/2.)*r^2*(1 + 8*π*r^2*p)*(p + ε)*(-((3*ζ + 4*η)*(-1 + 8*π*r^2*ε)*Ξ_t0[i]*cs^2) + 3*exp(ν/2.)*ξ[i]*(1 + 8*π*r^2*p + (1 - 8*π*r^2*ε)*cs^2)))/(r - 2*m)^2 + 60*exp(ν)*r*p*ξ[i]*cs*cs_prime + 60*exp(ν)*r*ε*ξ[i]*cs*cs_prime + (12*exp(ν)*r^2*p*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*cs*cs_prime)/h + (12*exp(ν)*r^2*ε*(3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*cs*cs_prime)/h - (30*exp(ν/2.)*ζ*(p + ε)*Ξ_t0[i]*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(r - 2*m) + (8*exp(ν/2.)*η*(p + ε)*Ξ_t0[i]*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(r - 2*m) - (6*exp(ν/2.)*r*ζ*(p + ε)*(3*Ξ_t0[i] + Ξ_t0[i-2] - 4*Ξ_t0[i-1])*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(h*(r - 2*m)) - (8*exp(ν/2.)*r*η*(p + ε)*(3*Ξ_t0[i] + Ξ_t0[i-2] - 4*Ξ_t0[i-1])*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))/(h*(r - 2*m)) + (r*(p + ε)*((6*r^2*(Ξ_tm1[i] - Ξ_tp1[i]))/k - (3*exp(ν)*r*((3*ξ[i] + ξ[i-2] - 4*ξ[i-1])*(1 + 8*π*r^2*p + (1 - 8*π*r^2*ε)*cs^2) + 4*h*ξ[i]*(8*π*r*(p - ε*cs^2) + (1 + 8*π*r^2*p)*cs*cs_prime)))/h - (exp(ν/2.)*(r*(r - 2*m)*(3*ζ + 4*η)*(1 - 8*π*r^2*ε)*(3*Ξ_t0[i] + Ξ_t0[i-2] - 4*Ξ_t0[i-1])*cs^2 + 2*h*Ξ_t0[i]*(12*(r - 2*m)*η*(1 + 16*π*r^2*p)*cs^2 - 16*π*r^2*(r - 2*m)*(3*ζ - 2*η)*ε*cs^2 - (3*ζ + 4*η)*(1 + 8*π*r^2*p)*(4*π*r^3*p*(1 + cs^2) - 2*r^2*cs*cs_prime + m*(1 + cs^2 + 4*r*cs*cs_prime)))))/(h*(r - 2*m))))/(r - 2*m)))/(12.0*r^3*(p + ε))

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function initialize_solution_file(filename::String, m::AbstractVector{Float64}, p::AbstractVector{Float64}, ε::AbstractVector{Float64}, ν::AbstractVector{Float64}, r::AbstractVector{Float64}, cs::AbstractVector{Float64}, cs_prime::AbstractVector{Float64}, u1::AbstractVector{Float64}, v1::AbstractVector{Float64}, ds::Int64, N::Int, num_time_steps::Int, chunk_size::Int)
    fmode = "w"
    file = h5open(filename, fmode); # filename should include the path

    # create group for solution data
    solution_group_name = "solution"
    HDF5Helper.create_file_group!(file, solution_group_name);

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
    HDF5Helper.create_dataset!(file, solution_group_name, "xi", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "xi_dot", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "IR", Float64, dataspace, chunk);


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

# compute the k1, k2, k3, k4 terms for the RK4 time integration
@views function compute_ks!(u1::Vector{Float64}, v1::Vector{Float64}, u1_k1::Vector{Float64}, v1_k1::Vector{Float64}, u1_k2::Vector{Float64}, v1_k2::Vector{Float64}, u1_k3::Vector{Float64}, v1_k3::Vector{Float64}, k1::Matrix{Float64}, k2::Matrix{Float64}, k3::Matrix{Float64}, k4::Matrix{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, h::Float64, k::Float64, nPoints::Int64)
    ### BCs at r=0 for regularity: ξ[t,0] = u1[t, 0] = 0, ∂_{t}ξ[t,0] = v1[t, 0] = 0, so the equations are already solved at r = 0 ###

    # compute k1
    @inbounds for i = 2:nPoints-1
        k1[1, i] = k * v1[i]
        k1[2, i] = k * RK4Eq(i, u1, v1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

        u1_k1[i] = u1[i] + 0.5 * k1[1, i]
        v1_k1[i] = v1[i] + 0.5 * k1[2, i]
    end

    i = nPoints
    k1[1, i] = k * v1[i]
    k1[2, i] = k * RK4EqSurface(i, u1, v1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

    u1_k1[i] = u1[i] + 0.5 * k1[1, i]
    v1_k1[i] = v1[i] + 0.5 * k1[2, i]

    # compute k2
    @inbounds for i = 2:nPoints-1
        k2[1, i] = k * v1_k1[i]
        k2[2, i] = k * RK4Eq(i, u1_k1, v1_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

        u1_k2[i] = u1[i] + 0.5 * k2[1, i]
        v1_k2[i] = v1[i] + 0.5 * k2[2, i]
    end

    i = nPoints
    k2[1, i] = k * v1_k1[i]
    k2[2, i] = k * RK4EqSurface(i, u1_k1, v1_k1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

    u1_k2[i] = u1[i] + 0.5 * k2[1, i]
    v1_k2[i] = v1[i] + 0.5 * k2[2, i]    

    # compute k3
    @inbounds for i = 2:nPoints-1
        k3[1, i] = k * v1_k2[i]
        k3[2, i] = k * RK4Eq(i, u1_k2, v1_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

        u1_k3[i] = u1[i] + k3[1, i]
        v1_k3[i] = v1[i] + k3[2, i]
    end

    i = nPoints
    k3[1, i] = k * v1_k2[i]
    k3[2, i] = k * RK4EqSurface(i, u1_k2, v1_k2, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)

    u1_k3[i] = u1[i] + k3[1, i]
    v1_k3[i] = v1[i] + k3[2, i]

    # compute k4
    @inbounds for i = 2:nPoints-1
        k4[1, i] = k * v1_k3[i]
        k4[2, i] = k * RK4Eq(i, u1_k3, v1_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)
    end

    i = nPoints
    k4[1, i] = k * v1_k3[i]
    k4[2, i] = k * RK4EqSurface(i, u1_k3, v1_k3, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h)
end

# compute independent residual
@views function compute_IR!(sol_three_levels::AbstractArray, IR::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, v1_tm1::Vector{Float64}, u1_t0::Vector{Float64}, v1_t0::Vector{Float64}, v1_tp1::Vector{Float64}, η::Float64, ζ::Float64, h::Float64, k::Float64, nPoints::Int64)
    # fill arrays
    v1_tm1 .= sol_three_levels[1, 2, :];

    u1_t0 .= sol_three_levels[2, 1, :];
    v1_t0 .= sol_three_levels[2, 2, :];

    v1_tp1 .= sol_three_levels[3, 2, :];

    for i = 2:nPoints-1
        IR[i] = LeapFrogEq(i, u1_t0, v1_tm1, v1_t0, v1_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h, k)
    end

    i = nPoints
    IR[i] = LeapFrogEqSurface(i, u1_t0, v1_tm1, v1_t0, v1_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], r[i], η, ζ, h, k)
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

# note that we assume stationary initial data both here and in all time domain functions
function solve(star::NeutronStarOscillations.Star, h::Float64; print_progress::Bool = true)
    fname = TimeDomain.td_fname(star, h)
    if abs(star.ξ_ID(0.0)) > 1e-16
        error("|ξ(r=0)| > 1e-16. ξ must be zero at the center by regularity. ξ(0.0) = $(star.ξ_ID(0.0))")
    end

    # length scale not explicitly in functions for Eckart fluid; the inputs η and ζ are actually η * L and ζ * L
    ηTimesL = star.η * star.L
    ζTimesL = star.ζ * star.L

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
    sol_save = zeros(Float64, star.save_every+1, 3, nPointsSpace_save); # 2 perturbation variables, 1 IR (we do save_every + 1 to ensure that at the final save there is always sufficient space to tack on the solution at the final time step)
    spatial_N = TOV_length
    sol_three_levels = zeros(Float64, 3, 2, spatial_N); # 2 perturbation variables

    # set up arrays for IR computation
    v1_tm1 = zeros(Float64, spatial_N);
    u1_t0 = zeros(Float64, spatial_N);
    v1_t0 = zeros(Float64, spatial_N);
    v1_tp1 = zeros(Float64, spatial_N);

    ################ SET UP INITIAL DATA ################
    nPointsSpace = TOV_length;
    u1_0 = @. star.ξ_ID(r) # ξ
    v1_0 = zero(r); # ∂ξ/∂t — # assume stationary initial data

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

        compute_ks!(u1_n, v1_n, u1_k1, v1_k1, u1_k2, v1_k2, u1_k3, v1_k3, k1, k2, k3, k4, m, p, ε, ν, cs, cs_prime, r, ηTimesL, ζTimesL, h, k, nPointsSpace)

        for i in 1:nPointsSpace
            u1_np1[i] = u1_n[i] + (k1[1, i] + 2k2[1, i] + 2k3[1, i] + k4[1, i]) / 6.0
            v1_np1[i] = v1_n[i] + (k1[2, i] + 2k2[2, i] + 2k3[2, i] + k4[2, i]) / 6.0
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
            compute_IR!(sol_three_levels, IR, m, p, ε, ν, cs, cs_prime, r, v1_tm1, u1_t0, v1_t0, v1_tp1, ηTimesL, ζTimesL, h, k, nPointsSpace)
            saved_idx_counter += 1;
            update_solution_array!(sol_save, sol_three_levels, IR, saved_idx_counter, ds_fact, nPointsSpace)
  
            if saved_idx_counter == star.save_every
                HDF5Helper.save_to_file_ECKART!(file, sol_save, save_idx0, save_idx0 + star.save_every- 1)
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
            HDF5Helper.save_final_chunk_to_file_ECKART!(file, sol_save, save_idx0, save_idx0 + saved_idx_counter - 1, saved_idx_counter, saved_time .* km_to_ms)
        end
    end
end

end

module BDNK
using Printf
using HDF5
using Dierckx
using SparseArrays
using LinearAlgebra
using NeutronStarOscillations
using ..TimeDomain
using ...HDF5Helper
using ...TOV
using ...BDNKDiffEqCoefficients
using ...BDNKLinearSystem
using ...BDNKCharacteristicSpeeds
using ...BDNKIndependentResidual
using ...BDNKInitialData
using ...FiniteDiffOrder4

const c = 299792458;
const sec_to_km = c * 1e-3;
one_norm(x::AbstractArray) = sum(abs, x) / length(x)

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function initialize_solution_file(filename::String, m::AbstractVector{Float64}, p::AbstractVector{Float64}, ε::AbstractVector{Float64}, ν::AbstractVector{Float64}, r::AbstractVector{Float64}, cs::AbstractVector{Float64}, cs_prime::AbstractVector{Float64}, Λ0::AbstractVector{Float64}, Λ1::AbstractVector{Float64}, Λ2::AbstractVector{Float64}, u1::AbstractVector{Float64}, u2::AbstractVector{Float64}, u3::AbstractVector{Float64}, v1::AbstractVector{Float64}, v2::AbstractVector{Float64}, v3::AbstractVector{Float64}, Q_u2::AbstractVector{Float64}, Q_u3::AbstractVector{Float64},  ds::Int64, N::Int, num_time_steps::Int, chunk_size::Int)
    fmode = "w"
    file = h5open(filename, fmode); # filename should include the path

    # create group for solution data
    solution_group_name = "solution"
    HDF5Helper.create_file_group!(file, solution_group_name);

    # first store simulation parameters
    param_group_name = "parameters"
    # HDF5Helper.create_file_group!(file, param_group_name);
    # for (key, value) in param_dict
    #     file["parameters"][key] = value
    # end

    # save data fixed at run time
    file["solution"]["r"] = r[1:ds:end];
    file["solution"]["m"] = m[1:ds:end];
    file["solution"]["p"] = p[1:ds:end];
    file["solution"]["ε"] = ε[1:ds:end];
    file["solution"]["ν"] = ν[1:ds:end];
    file["solution"]["cs"] = cs[1:ds:end];
    file["solution"]["cs_prime"] = cs_prime[1:ds:end];
    file["solution"]["Q_u2"] = Q_u2[1:ds:end];
    file["solution"]["Q_u3"] = Q_u3[1:ds:end];
    file["solution"]["Λ0"] = Λ0[1:ds:end];
    file["solution"]["Λ1"] = Λ1[1:ds:end];
    file["solution"]["Λ2"] = Λ2[1:ds:end];


    # create datasets for the solution
    dataspace = ((num_time_steps+1, N), (num_time_steps+1, N))
    chunk = (chunk_size, N)
    HDF5Helper.create_dataset!(file, solution_group_name, "du", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "deps", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "dlambda", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "du_dot", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "deps_dot", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "dlambda_dot", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "Constrained_IR1", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "Constrained_IR2", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "Wave_IR1", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "Wave_IR2", Float64, dataspace, chunk);
    HDF5Helper.create_dataset!(file, solution_group_name, "Wave_IR3", Float64, dataspace, chunk);



    # manually save initial data
    file[solution_group_name]["du"][1, :] = u1[1:ds:end];
    file[solution_group_name]["deps"][1, :] = u2[1:ds:end];
    file[solution_group_name]["dlambda"][1, :] = u3[1:ds:end];
    file[solution_group_name]["du_dot"][1, :] = v1[1:ds:end];
    file[solution_group_name]["deps_dot"][1, :] = v2[1:ds:end];
    file[solution_group_name]["dlambda_dot"][1, :] = v3[1:ds:end];
    file[solution_group_name]["Constrained_IR1"][1, :] = zero(v1[1:ds:end]) * NaN; # don't compute IR at initial time
    file[solution_group_name]["Constrained_IR2"][1, :] = zero(v1[1:ds:end]) * NaN;
    file[solution_group_name]["Wave_IR1"][1, :] = zero(v1[1:ds:end]) * NaN;
    file[solution_group_name]["Wave_IR2"][1, :] = zero(v1[1:ds:end]) * NaN;
    file[solution_group_name]["Wave_IR3"][1, :] = zero(v1[1:ds:end]) * NaN;
    return file
end

@views function update_solution_array!(sol::AbstractArray{Float64}, u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64}, IR1::Vector{Float64}, IR2::Vector{Float64}, IR3::Vector{Float64}, IR4::Vector{Float64}, IR5::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = u1[1:ds_fact:nMax];
    sol[time_idx, 2, :] = u2[1:ds_fact:nMax];
    sol[time_idx, 3, :] = u3[1:ds_fact:nMax];
    sol[time_idx, 4, :] = v1[1:ds_fact:nMax];
    sol[time_idx, 5, :] = v2[1:ds_fact:nMax];
    sol[time_idx, 6, :] = v3[1:ds_fact:nMax];
    sol[time_idx, 7, :] = IR1[1:ds_fact:nMax];
    sol[time_idx, 8, :] = IR2[1:ds_fact:nMax];
    sol[time_idx, 9, :] = IR3[1:ds_fact:nMax];
    sol[time_idx, 10, :] = IR4[1:ds_fact:nMax];
    sol[time_idx, 11, :] = IR5[1:ds_fact:nMax];
end

@views function update_solution_array!(sol::AbstractArray{Float64}, sol_three_levels::AbstractArray{Float64}, IR1::Vector{Float64}, IR2::Vector{Float64}, IR3::Vector{Float64}, IR4::Vector{Float64}, IR5::Vector{Float64}, time_idx::Int64, ds_fact::Int64, nMax::Int64)
    sol[time_idx, 1, :] = sol_three_levels[2, 1, 1:ds_fact:nMax];
    sol[time_idx, 2, :] = sol_three_levels[2, 2, 1:ds_fact:nMax];
    sol[time_idx, 3, :] = sol_three_levels[2, 3, 1:ds_fact:nMax];
    sol[time_idx, 4, :] = sol_three_levels[2, 4, 1:ds_fact:nMax];
    sol[time_idx, 5, :] = sol_three_levels[2, 5, 1:ds_fact:nMax];
    sol[time_idx, 6, :] = sol_three_levels[2, 6, 1:ds_fact:nMax];
    sol[time_idx, 7, :] = IR1[1:ds_fact:nMax];
    sol[time_idx, 8, :] = IR2[1:ds_fact:nMax];
    sol[time_idx, 9, :] = IR3[1:ds_fact:nMax];
    sol[time_idx, 10, :] = IR4[1:ds_fact:nMax];
    sol[time_idx, 11, :] = IR5[1:ds_fact:nMax];
end


# compute independent residuals
@views function compute_IR!(sol_three_levels::AbstractArray, IR1::Vector{Float64}, IR2::Vector{Float64}, IR3::Vector{Float64}, IR4::Vector{Float64}, IR5::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, cs_prime_prime::Vector{Float64}, r::Vector{Float64}, u1_tm1::Vector{Float64}, u2_tm1::Vector{Float64}, u3_tm1::Vector{Float64}, v1_tm1::Vector{Float64}, v2_tm1::Vector{Float64}, v3_tm1::Vector{Float64}, u1_t0::Vector{Float64}, u2_t0::Vector{Float64}, u3_t0::Vector{Float64}, w1_t0::Vector{Float64}, w2_t0::Vector{Float64}, w3_t0::Vector{Float64}, v1_t0::Vector{Float64}, v2_t0::Vector{Float64}, v3_t0::Vector{Float64}, u1_tp1::Vector{Float64}, u2_tp1::Vector{Float64}, u3_tp1::Vector{Float64}, v1_tp1::Vector{Float64}, v2_tp1::Vector{Float64}, v3_tp1::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64, k::Float64, nPoints::Int64, CNSystem)
    # fill arrays
    u1_tm1 .= sol_three_levels[1, 1, :];
    u2_tm1 .= sol_three_levels[1, 2, :];
    u3_tm1 .= sol_three_levels[1, 3, :];
    v1_tm1 .= sol_three_levels[1, 4, :];
    v2_tm1 .= sol_three_levels[1, 5, :];
    v3_tm1 .= sol_three_levels[1, 6, :];

    u1_t0 .= sol_three_levels[2, 1, :];
    u2_t0 .= sol_three_levels[2, 2, :];
    u3_t0 .= sol_three_levels[2, 3, :];
    v1_t0 .= sol_three_levels[2, 4, :];
    v2_t0 .= sol_three_levels[2, 5, :];
    v3_t0 .= sol_three_levels[2, 6, :];

    u1_tp1 .= sol_three_levels[3, 1, :];
    u2_tp1 .= sol_three_levels[3, 2, :];
    u3_tp1 .= sol_three_levels[3, 3, :];
    v1_tp1 .= sol_three_levels[3, 4, :];
    v2_tp1 .= sol_three_levels[3, 5, :];
    v3_tp1 .= sol_three_levels[3, 6, :];

    # compute w's
    FiniteDiffOrder4.compute_first_derivative(w1_t0, u1_t0, h, nPoints);
    FiniteDiffOrder4.compute_first_derivative(w2_t0, u2_t0, h, nPoints);
    FiniteDiffOrder4.compute_first_derivative(w3_t0, u3_t0, h, nPoints);
    
    # independent residual at origin only defined for the one equation we solve there
    IR1[1] = BDNKIndependentResidual.ConstrainedSystem.Center.Eq1(v2_tp1, v2_tm1, u1_t0, u2_t0, u3_t0, v1_t0, v2_t0, CNSystem, h, k)
    IR2[1] = 0.0
    IR3[1] = 0.0
    IR4[1] = 0.0
    IR5[1] = 0.0


    for i = 2:nPoints-1
        IR1[i] = BDNKIndependentResidual.ConstrainedSystem.Interior.Eq1(v1_tp1, v2_tp1, v1_tm1, v2_tm1, u1_t0, u2_t0, u3_t0, v1_t0, v2_t0, CNSystem, i, h, k)
        IR2[i] = BDNKIndependentResidual.ConstrainedSystem.Interior.Eq2(v1_tp1, v2_tp1, v1_tm1, v2_tm1, u1_t0, u2_t0, u3_t0, v1_t0, v2_t0, CNSystem, i, h, k)
        IR3[i] = BDNKIndependentResidual.WaveSystem.Interior.Eq1(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
        IR4[i] = BDNKIndependentResidual.WaveSystem.Interior.Eq2(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
        IR5[i] = BDNKIndependentResidual.WaveSystem.Interior.Eq3(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
    end

    i = nPoints
    IR1[i] = BDNKIndependentResidual.ConstrainedSystem.Surface.Eq1(v1_tp1, v2_tp1, v1_tm1, v2_tm1, u1_t0, u2_t0, u3_t0, v1_t0, v2_t0, CNSystem, i, h, k)
    IR2[i] = BDNKIndependentResidual.ConstrainedSystem.Surface.Eq2(v1_tp1, v2_tp1, v1_tm1, v2_tm1, u1_t0, u2_t0, u3_t0, v1_t0, v2_t0, CNSystem, i, h, k)
    IR3[i] = BDNKIndependentResidual.WaveSystem.Surface.Eq1(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
    IR4[i] = BDNKIndependentResidual.WaveSystem.Surface.Eq2(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
    IR5[i] = BDNKIndependentResidual.WaveSystem.Surface.Eq3(i, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, v1_tp1, v2_tp1, v3_tp1, m[i], p[i], ε[i], ν[i], cs[i], cs_prime[i], cs_prime_prime[i], r[i], η, ζ, τε, τP, τQ, L, h, k)
end

# Kreiss-Oliger dissipation
function kreiss_oliger(var::AbstractVector{Float64}, j::Int, coef::Float64, nPoints::Int64, spacing::Int64)
    if j == 1
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[1]
        u_plus_1 = var[2]
        u_plus_2 = var[3]
        error("Kreiss-Oliger dissipation not defined at first grid point")
    elseif j == 2
        u_minus_2 = var[1]
        u_minus_1 = var[1]
        u = var[2]
        u_plus_1 = var[3]
        u_plus_2 = var[4]
        error("Kreiss-Oliger dissipation not defined at second grid point")
    elseif j == nPoints-1
        u_minus_2 = var[nPoints-3]
        u_minus_1 = var[nPoints-2]
        u = var[nPoints-1]
        u_plus_1 = var[nPoints]
        u_plus_2 = 0.0
        error("Kreiss-Oliger dissipation not defined at second to last grid point")
    elseif j == nPoints
        u_minus_2 = var[nPoints-2]
        u_minus_1 = var[nPoints-1]
        u = var[nPoints]
        u_plus_1 = 0.0
        u_plus_2 = 0.0
        error("Kreiss-Oliger dissipation not defined at last grid point")
    else
        u_minus_2 = var[j-2*spacing]
        u_minus_1 = var[j-1*spacing]
        u = var[j]
        u_plus_1 = var[j+1*spacing]
        u_plus_2 = var[j+2*spacing]
    end

    return coef * (u_plus_2 - 4.0 * u_plus_1 + 6.0 * u - 4.0 * u_minus_1 + u_minus_2) / 16.0
end

# note that we assume stationary initial data both here and in all time domain functions
function solve(star::NeutronStarOscillations.Star, h::Float64; δε_center::Float64 = 1000.0, print_progress::Bool = true)
    fname = TimeDomain.td_fname(star, h)
    # ensure boundary condition that u1(r=0) = 0 is enforced
    if abs(star.δu_ID(0.0)) > 1e-16
        error("|δu(r=0)| > 1e-16. δu must be zero at the center by regularity. δu(0.0) = $(star.δu_ID(0.0))")
    end

    ################ PRELIMINARIES ################
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

    ################ SOLVE FOR INITIAL DATA AND TOV BACKGROUND ################
    # run at high resolutions
    h_ID = 5e-4
    r_1, u2_lev1, u3_lev1 = BDNKInitialData.compute_initial_data(star, 4.0h_ID; u2_0=δε_center);
    r_2, u2_lev2, u3_lev2 = BDNKInitialData.compute_initial_data(star, 2.0h_ID; u2_0=δε_center);
    r_3, u2_lev3, u3_lev3 = BDNKInitialData.compute_initial_data(star, 1.0h_ID; u2_0=δε_center);

    # interpolate to compute convergence factor
    spline_order = 5; s = 0.0;
    δε_spline_1 = Spline1D(r_1, u2_lev1; k=spline_order, s=s);
    δε_spline_2 = Spline1D(r_2, u2_lev2; k=spline_order, s=s);
    δε_spline_3 = Spline1D(r_3, u2_lev3; k=spline_order, s=s);

    δλ_spline_1 = Spline1D(r_1, u3_lev1; k=spline_order, s=s);
    δλ_spline_2 = Spline1D(r_2, u3_lev2; k=spline_order, s=s);
    δλ_spline_3 = Spline1D(r_3, u3_lev3; k=spline_order, s=s);

    # compute high resolution TOV and desample to user input
    h_TOV = 1e-4;
    r, u2_0, u3_0, u1_0, w1_0, m, p, ε, ν, cs, cs_prime = BDNKInitialData.compute_initial_data(star, h_TOV / 2.0; u2_0=δε_center, return_all=true);
    cs_prime[1] = 0.0
    cs_prime_prime = zero(cs_prime)
    FiniteDiffOrder4.compute_first_derivative(cs_prime_prime, cs_prime, diff(r)[1], length(r));

    # normalize initial data so that max(abs.(δε)) = 1.0. Note that we normalize u2 and u3 since we solve a linear ODE for them at t=0, so we can freely rescale them
    max_abs_δε = maximum(abs.(u2_0));
    u2_0 .= u2_0 ./ max_abs_δε;
    u3_0 .= u3_0 ./ max_abs_δε;

    w1_0 = zero(u1_0);
    w2_0 = zero(u2_0);
    w3_0 = zero(u3_0);
    FiniteDiffOrder4.compute_first_derivative(w1_0, u1_0, diff(r)[1], length(r));
    FiniteDiffOrder4.compute_first_derivative(w2_0, u2_0, diff(r)[1], length(r));
    FiniteDiffOrder4.compute_first_derivative(w3_0, u3_0, diff(r)[1], length(r));

    # downsample
    ds_fact = argmin(@. abs(r - h)) - 1;
    r = r[1:ds_fact:end];
    diff(r)[1] ≈ h ? nothing : error("Grid spacing does not match desired value of h after downsampling");
    m = m[1:ds_fact:end];
    p = p[1:ds_fact:end];
    ε = ε[1:ds_fact:end];
    ν = ν[1:ds_fact:end];
    cs = cs[1:ds_fact:end];
    cs_prime = cs_prime[1:ds_fact:end];
    cs_prime_prime = cs_prime_prime[1:ds_fact:end];
    u1_0 = u1_0[1:ds_fact:end];
    u2_0 = u2_0[1:ds_fact:end];
    u3_0 = u3_0[1:ds_fact:end];
    w1_0 = w1_0[1:ds_fact:end];
    w2_0 = w2_0[1:ds_fact:end];
    w3_0 = w3_0[1:ds_fact:end];
    TOV_length = length(m);

    # intial data convergence factor
    δλ_1_ds = @. δλ_spline_1(r);
    δλ_2_ds = @. δλ_spline_2(r);
    δλ_3_ds = @. δλ_spline_3(r);

    δε_1_ds = @. δε_spline_1(r);
    δε_2_ds = @. δε_spline_2(r);
    δε_3_ds = @. δε_spline_3(r);

    Q_δλ = [abs(δλ_2_ds[i] - δλ_1_ds[i]) / abs(δλ_3_ds[i] - δλ_2_ds[i]) for i in eachindex(δλ_1_ds)];
    Q_δε = [abs(δε_2_ds[i] - δε_1_ds[i]) / abs(δε_3_ds[i] - δε_2_ds[i]) for i in eachindex(δε_1_ds)];

    # compute background PDE coefficients 
    CNSystem = BDNKDiffEqCoefficients.CNSystem(m, p, ε, ν, cs, cs_prime, cs_prime_prime, r, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, h, TOV_length);
    
    # compute characteristic speeds——we compute the three terms in the expression for the characteristic speed separately in case one has chosen a frame in which one at least of the c^2 are imaginary
    Λ0 = zeros(Float64, TOV_length);
    Λ1 = zeros(Float64, TOV_length);
    Λ2 = zeros(Float64, TOV_length);

    for i in eachindex(Λ0)
        ΛΛ0, ΛΛ1, ΛΛ2 = BDNKCharacteristicSpeeds.compute(p[i], ε[i], cs[i], star.η, star.ζ, star.τε, star.τP, star.τQ, star.L)
        Λ0[i] = ΛΛ0
        Λ1[i] = ΛΛ1
        Λ2[i] = ΛΛ2
    end

    # we always assume stationary initial data
    v1_0 = zero(u1_0)
    v2_0 = zero(u1_0)

    ################ SET UP ARRAYS FOR TIME EVOLUTION ################
    # we have 5 TOV variables, 5 perturbation variables (velocity, energy and metric perturbations and the time derivatives of the former 2) and 5 independent residuals (IRs——2 from the wave-equations we solve and 3 from the original wave-like equations that come out directly from the Einstein-BDNK system). To compute the IR, we must retain three time levels of the perturbations. We separately store these 3 time levels from another (desampled) solution array that will be saved to file
    ds_fact = argmin(@. abs(r - star.h_save)) - 1; # downsampling factor for saving data
    nPointsSpace_save = 1 + (TOV_length-1) ÷ ds_fact; # number of points to save
    sol_save = zeros(Float64, star.save_every+1, 11, nPointsSpace_save); # 5 perturbation variables, 5 IRs (we do save_every + 1 to ensure that at the final save there is always sufficient space to tack on the solution at the final time step)
    spatial_N = TOV_length
    sol_three_levels = zeros(Float64, 3, 6, spatial_N); # 5 perturbation variables but we separately compute and save the time derivative of δλ since this is needed for computing the independent residual

    # set up arrays for IR computation which couples three time levels
    u1_tm1 = zeros(Float64, spatial_N);
    u2_tm1 = zeros(Float64, spatial_N);
    u3_tm1 = zeros(Float64, spatial_N);

    v1_tm1 = zeros(Float64, spatial_N);
    v2_tm1 = zeros(Float64, spatial_N);
    v3_tm1 = zeros(Float64, spatial_N);

    u1_t0 = zeros(Float64, spatial_N);
    u2_t0 = zeros(Float64, spatial_N);
    u3_t0 = zeros(Float64, spatial_N);
    w1_t0 = zeros(Float64, spatial_N);
    w2_t0 = zeros(Float64, spatial_N);
    w3_t0 = zeros(Float64, spatial_N);
    v1_t0 = zeros(Float64, spatial_N);
    v2_t0 = zeros(Float64, spatial_N);
    v3_t0 = zeros(Float64, spatial_N);

    u1_tp1 = zeros(Float64, spatial_N);
    u2_tp1 = zeros(Float64, spatial_N);
    u3_tp1 = zeros(Float64, spatial_N);
    v1_tp1 = zeros(Float64, spatial_N);
    v2_tp1 = zeros(Float64, spatial_N);
    v3_tp1 = zeros(Float64, spatial_N);


    ################ SET UP AND SAVE INITIAL DATA ################
    nPointsSpace = TOV_length;

    if isfile(fname)
        print_progress ? println("File $fname already exists and will be overwritten.") : nothing
        rm(fname)
    end

    num_saved_steps = length(save_time_idx);
    v3_0 = zero(u3_0);
    BDNKLinearSystem.compute_v3!(u1_0, u2_0, u3_0, w1_0, w2_0, w3_0, v1_0, v2_0, v3_0, m, p, ε, ν, cs, r, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L)
    file = initialize_solution_file(fname, m, p, ε, ν, r, cs, cs_prime, Λ0, Λ1, Λ2, u1_0, u2_0, u3_0, v1_0, v2_0, v3_0, Q_δε, Q_δλ, ds_fact, nPointsSpace_save, num_saved_steps, star.save_every)

    ################ SET UP STEP ARRAYS FOR CURRENT AND ADVANCED TIME STEP ################
    u1_n = zeros(nPointsSpace);
    u2_n = zeros(nPointsSpace);
    u3_n = zeros(nPointsSpace);
    v1_n = zeros(nPointsSpace);
    v2_n = zeros(nPointsSpace);
    v3_n = zeros(nPointsSpace);

    u1_np1 = zeros(nPointsSpace);
    u2_np1 = zeros(nPointsSpace);
    u3_np1 = zeros(nPointsSpace);
    v1_np1 = zeros(nPointsSpace);
    v2_np1 = zeros(nPointsSpace);
    v3_np1 = zeros(nPointsSpace);

    u1_n[1:nPointsSpace] .= u1_0;
    u2_n[1:nPointsSpace] .= u2_0;
    u3_n[1:nPointsSpace] .= u3_0;
    v1_n[1:nPointsSpace] .= v1_0;
    v2_n[1:nPointsSpace] .= v2_0;
    v3_n[1:nPointsSpace] .= v3_0;

    ################ SET UP LINEAR SYSTEM WHICH IS SOLVED AT EACH TIME STEP FOR THE SOLUTION AT THE ADVANCED TIME STEP ################
    # set up linear system arrays — note that we do not solve for (u1, v1, u3) at j = 1 due to regularity conditions, hence the -3 in the sizes below
    linear_system_RHS = zeros(5 * nPointsSpace - 3); # five equations at each spatial point j = 2:n (three non-trivial BDNK equations and two trivial time reduction equations)
    discretized_eqs = zeros(5 * nPointsSpace - 3);
    jacobian = zeros(5 * nPointsSpace - 3, 5 * nPointsSpace - 3);

    total_u1_vars = nPointsSpace-1
    total_u2_vars = nPointsSpace
    total_u3_vars = nPointsSpace-1
    total_u_vars = total_u1_vars + total_u2_vars + total_u3_vars

    total_v1_vars = nPointsSpace-1
    total_v2_vars = nPointsSpace
    total_v_vars = total_v1_vars + total_v2_vars

    total_vars = total_u_vars + total_v_vars

    # dummy arrays to compute the RHS of the linear system (i.e., the terms in the discretized equations which depend on the current time step)
    u1_np1_dummy = zeros(nPointsSpace);
    u2_np1_dummy = zeros(nPointsSpace);
    u3_np1_dummy = zeros(nPointsSpace);
    v1_np1_dummy = zeros(nPointsSpace);
    v2_np1_dummy = zeros(nPointsSpace);

    # independent residual arrays
    IR1 = zero(r);
    IR2 = zero(r);
    IR3 = zero(r);
    IR4 = zero(r);
    IR5 = zero(r);

    # parameters for time integration loop
    saved_idx_counter = 0;
    save_idx0 = 2; # first index to save to in solution array (we skip the initial data which is already saved)
    saved_time = Float64[0.0];
    eq_resid = -2.0;
    eq_resids = Float64[NaN];

    # compute jacobian, turn it into a sparse matrix and factorize it
    BDNKLinearSystem.Jacobian.compute_jacobian!(jacobian, CNSystem, h, k, nPointsSpace);
    sparse_jac = sparse(jacobian);
    jac_fact = factorize(sparse_jac);
    sol = zeros(size(jacobian, 1));

    ################ TIME INTEGRATION ################
    for i = 1:total_num_time_steps
        print_string = "Completion: $(round(100 * i/(total_num_time_steps); digits=5))%\r"  
        print_progress ? print(print_string) : nothing

        # Add KO dissipation. Utilize the advanced time step arrays to add the KO dissipation (this is because we do not want to do something like u1_n[j] += -kreiss_oliger(u1_n, j, coef) since we do not want to edit in place)
        u1_np1 .= u1_n;
        u2_np1 .= u2_n;
        u3_np1 .= u3_n;
        v1_np1 .= v1_n;
        v2_np1 .= v2_n;

        KO_edge = 3 # start applying KO dissipation at this grid point to avoid boundary issues
        @inbounds for j in KO_edge:nPointsSpace-KO_edge+1 # don't apply at all points
            u1_np1[j] += -kreiss_oliger(u1_n, j, star.KO, nPointsSpace, 1)
            u2_np1[j] += -kreiss_oliger(u2_n, j, star.KO, nPointsSpace, 1)
            u3_np1[j] += -kreiss_oliger(u3_n, j, star.KO, nPointsSpace, 1)
            v1_np1[j] += -kreiss_oliger(v1_n, j, star.KO, nPointsSpace, 1)
            v2_np1[j] += -kreiss_oliger(v2_n, j, star.KO, nPointsSpace, 1)
        end

        u1_n .= u1_np1;
        u2_n .= u2_np1;
        u3_n .= u3_np1;
        v1_n .= v1_np1;
        v2_n .= v2_np1;

        # check regularity conditions at center satisfies
        if abs(u1_np1[1]) > 1e-16 || abs(v1_np1[1]) > 1e-16 || abs(u3_np1[1]) > 1e-16
            error("u1 or v1 at center is non-zero, violating regularity condition. u1_np1[1] = $(u1_np1[1]), v1_np1[1] = $(v1_np1[1]), u1_np1_dummy[1] = $(u1_np1_dummy[1]), v1_np1_dummy[1] = $(v1_np1_dummy[1])")
        end

        # compute linear system
        BDNKLinearSystem.RHS.compute_RHS!(linear_system_RHS, u1_np1_dummy, u2_np1_dummy, u3_np1_dummy, v1_np1_dummy, v2_np1_dummy, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, h, k, nPointsSpace);

        # solve linear system and update solution at advanced time step
        ldiv!(sol, jac_fact, linear_system_RHS)
        u1_np1[2:end] .= sol[1:total_u1_vars];
        u2_np1 .= sol[total_u1_vars + 1:total_u1_vars+total_u2_vars];
        u3_np1[2:end] .= sol[total_u1_vars + total_u2_vars + 1:total_u_vars];
        v1_np1[2:end] .= sol[total_u_vars + 1:total_u_vars + total_v1_vars];
        v2_np1 .= sol[total_u_vars + total_v1_vars + 1:total_vars];

        # eliminated v3 from the equations but compute it separately here since it's necessary to evaluate the independent residual for the wave-like (i.e., not constrained) system. We note that v3 contains some errors at the boundary due to computation of radial derivatives of u2, u3, but the wave-like IRs tend to converge there with resolution 
        FiniteDiffOrder4.compute_first_derivative(w1_0, u1_np1, h, nPointsSpace);
        FiniteDiffOrder4.compute_first_derivative(w2_0, u2_np1, h, nPointsSpace);
        FiniteDiffOrder4.compute_first_derivative(w3_0, u3_np1, h, nPointsSpace);
        BDNKLinearSystem.compute_v3!(u1_np1, u2_np1, u3_np1, w1_0, w2_0, w3_0, v1_np1, v2_np1, v3_np1, m, p, ε, ν, cs, r, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L);
        
        # compute residual from solved linear system
        BDNKLinearSystem.RHS.compute_RHS!(discretized_eqs, u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, h, k, nPointsSpace);
        eq_resid = one_norm(discretized_eqs);

        # update three-level solution array
        sol_three_levels[1, :, :] .= sol_three_levels[2, :, :];
        sol_three_levels[2, :, :] .= sol_three_levels[3, :, :];
        sol_three_levels[3, 1, :] .= u1_np1[:];
        sol_three_levels[3, 2, :] .= u2_np1[:];
        sol_three_levels[3, 3, :] .= u3_np1[:];
        sol_three_levels[3, 4, :] .= v1_np1[:];
        sol_three_levels[3, 5, :] .= v2_np1[:];
        sol_three_levels[3, 6, :] .= v3_np1[:];

        # compute IR and send solution to file if at a save time (note that we use i - 1 instead of i because we can only compute the IR with the advanced time step)
        if i - 1 ∈ save_time_idx
            append!(saved_time, time_steps[i - 1])
            append!(eq_resids, eq_resid)
            compute_IR!(sol_three_levels, IR1, IR2, IR3, IR4, IR5, m, p, ε, ν, cs, cs_prime, cs_prime_prime, r, u1_tm1, u2_tm1, u3_tm1, v1_tm1, v2_tm1, v3_tm1, u1_t0, u2_t0, u3_t0, w1_t0, w2_t0, w3_t0, v1_t0, v2_t0, v3_t0, u1_tp1, u2_tp1, u3_tp1, v1_tp1, v2_tp1, v3_tp1, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L, h, k, nPointsSpace, CNSystem)
            saved_idx_counter += 1;
            update_solution_array!(sol_save, sol_three_levels, IR1, IR2, IR3, IR4, IR5, saved_idx_counter, ds_fact, nPointsSpace)
            if saved_idx_counter == star.save_every
                HDF5Helper.save_to_file_BDNK!(file, sol_save, save_idx0, save_idx0 + star.save_every- 1)
                saved_idx_counter = 0;
                save_idx0 += star.save_every;
            end
        end

        # initialization for next time step
        u1_n .= u1_np1;
        u2_n .= u2_np1;
        u3_n .= u3_np1;
        v1_n .= v1_np1;
        v2_n .= v2_np1;
        v3_n .= v3_np1;

        # at final time step, save solution manually if not already saved
        if i == total_num_time_steps
            # do not compute IR at final time step so save solution here manually
            if i ∈ save_time_idx
                saved_idx_counter += 1;
                append!(saved_time, time_steps[save_time_idx[end]])
                append!(eq_resids, eq_resid)
                IR1 .= NaN;
                IR2 .= NaN;
                IR3 .= NaN;
                IR4 .= NaN;
                IR5 .= NaN;
                update_solution_array!(sol_save, u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, v3_np1, IR1, IR2, IR3, IR4, IR5, saved_idx_counter, ds_fact, nPointsSpace)
            end
            HDF5Helper.save_final_chunk_to_file_BDNK!(file, sol_save, save_idx0, save_idx0 + saved_idx_counter - 1, saved_idx_counter, saved_time .* km_to_ms, eq_resids)
        end
    end
end

end


end
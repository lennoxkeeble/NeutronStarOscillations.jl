#=

    Module comprising of helper functions which solve the Tolman-Oppenheimer-Volkoff (TOV) equations. Two solvers are implemented: an implicit Crank-Nicholson solver
    (second-order accurate) and an explicit RK4 solver. Functions to plot the results and compute convergence factors are also included.

=#

module TOV
using HDF5
using Printf
using ..PlotSettings
using LaTeXStrings
using CairoMakie
using NeutronStarOscillations

function compute_convergence_factor(star::NeutronStarOscillations.Star, hmin::Float64, type::String)
    if isequal(type, "Explicit")
        r_4h, mass_4h, pressure_4h, eps_4h, nu_4h = Explicit.load(star, 4.0hmin);
        r_2h, mass_2h, pressure_2h, eps_2h, nu_2h = Explicit.load(star, 2.0hmin);
        r_h, mass_h, pressure_h, eps_h, nu_h = Explicit.load(star, hmin);
    elseif type == "Implicit"
        r_4h, mass_4h, pressure_4h, eps_4h, nu_4h = Implicit.load(star, 4.0hmin);
        r_2h, mass_2h, pressure_2h, eps_2h, nu_2h = Implicit.load(star, 2.0hmin);
        r_h, mass_h, pressure_h, eps_h, nu_h = Implicit.load(star, hmin);
    else
        throw(ArgumentError("type must be either 'Explicit' or 'Implicit'"))
    end
    
    max_r = minimum([r_4h[end], r_2h[end], r_h[end]])

    r_max_idx_4h = length(r_4h[r_4h .<= max_r])
    r_max_idx_2h = 1 + 2 * (r_max_idx_4h - 1)
    r_max_idx_h = 1 + 4 * (r_max_idx_4h - 1)

    # check radii are equal
    if !isapprox(r_4h[r_max_idx_4h], r_2h[r_max_idx_2h], atol = 1e-8) || !isapprox(r_4h[r_max_idx_4h], r_h[r_max_idx_h], atol = 1e-8)
        println("radius (4h) = $(r_4h[r_max_idx_4h]); radius (2h) = $(r_2h[r_max_idx_2h]); radius (h) = $(r_h[r_max_idx_h]); diff_2h = $(abs(r_4h[r_max_idx_4h] - r_2h[r_max_idx_2h])); diff_h = $(abs(r_4h[r_max_idx_4h] - r_h[r_max_idx_h]))")
        throw(DomainError("radii values are not equal for (h, 2h, 4h) pair"))
    end

    # mask radii
    r_4h = r_4h[1:r_max_idx_4h]
    r_2h = r_2h[1:r_max_idx_2h]
    r_h = r_h[1:r_max_idx_h]

    Q_mass = zeros(length(r_4h))
    Q_pressure = zeros(length(r_4h))
    Q_eps = zeros(length(r_4h))
    Q_nu = zeros(length(r_4h))

    for i in eachindex(r_4h)
        radius_idx_2h = 1 + 2 * (i - 1)
        radius_idx_h = 1 + 4 * (i - 1)    
        if !isapprox(r_2h[radius_idx_2h], r_h[radius_idx_h], atol = 1e-8)
            throw(DomainError("radius values are not equal for (h, 2h) pair. 
                radius (2h) = $(r_2h[radius_idx_2h]); radius (h) = $(r_h[radius_idx_h]); diff = $(abs(r_2h[radius_idx_2h] - r_h[radius_idx_h]))"))
        elseif !isapprox(r_4h[i], r_h[radius_idx_h], atol = 1e-8)
            throw(DomainError("radius values are not equal for (h, 4h) pair. 
                radius (4h) = $(r_4h[i]); radius (h) = $(r_h[radius_idx_h]); diff = $(abs(r_4h[i] - r_h[radius_idx_h]))"))
        end

        abs_2h_4h_mass = abs(mass_4h[i] - mass_2h[radius_idx_2h])
        abs_2h_4h_pressure = abs(pressure_4h[i] - pressure_2h[radius_idx_2h])
        abs_2h_4h_eps = abs(eps_4h[i] - eps_2h[radius_idx_2h])
        abs_2h_4h_nu = abs(nu_4h[i] - nu_2h[radius_idx_2h])

        abs_h_2h_mass = abs(mass_2h[radius_idx_2h] - mass_h[radius_idx_h])
        abs_h_2h_pressure = abs(pressure_2h[radius_idx_2h] - pressure_h[radius_idx_h])
        abs_h_2h_eps = abs(eps_2h[radius_idx_2h] - eps_h[radius_idx_h])
        abs_h_2h_nu = abs(nu_2h[radius_idx_2h] - nu_h[radius_idx_h])

        Q_mass[i] = abs_2h_4h_mass / abs_h_2h_mass
        Q_pressure[i] = abs_2h_4h_pressure / abs_h_2h_pressure
        Q_eps[i] = abs_2h_4h_eps / abs_h_2h_eps
        Q_nu[i] = abs_2h_4h_nu / abs_h_2h_nu
    end

    return r_4h, Q_mass, Q_pressure, Q_eps, Q_nu
end


var_plot_fname(var::String, star::NeutronStarOscillations.Star, h::Float64)::String = star.fig_path * "TOV_" * var * @sprintf("_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_h_%.1e.png", star.kappa, star.n, star.εc_SI, star.ptol, h);
function plot_var(var::String, star::NeutronStarOscillations.Star, h::Float64, type::String; 
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

    if isequal(type, "Explicit")
        r, mass, pressure, eps, nu = Explicit.load(star, h);
    elseif type == "Implicit"
        r, mass, pressure, eps, nu = Implicit.load(star, h);
    else
        throw(ArgumentError("type must be either 'Explicit' or 'Implicit'"))
    end

    if isequal(var, "mass")
        y = [mass[1:desample_factor:end]] / NeutronStarOscillations.Msun_to_km;
    elseif isequal(var, "pressure")
        y = [pressure[1:desample_factor:end]] / star.pc_SI;
    elseif isequal(var, "energy_density")
        y = [eps[1:desample_factor:end]] / star.εc_SI;
    elseif isequal(var, "nu")
        y = [nu[1:desample_factor:end]] / nu[end];
    else
        throw(ArgumentError("var must be one of 'mass', 'pressure', 'energy_density', or 'nu'"))
    end

    x = [r[1:desample_factor:end]];
    fname = var_plot_fname(var, star, h);

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


# implicit Crank Nicholson solver (second-order accurate)# explicit RK4 solver (fourth-order accurate)
module Explicit
using HDF5
using DifferentialEquations
using StaticArrays
using Printf
using NeutronStarOscillations

fname(star::NeutronStarOscillations.Star, h::Float64; TD::Bool=false)::String = star.data_path * @sprintf("Explicit_TOV_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_h_%.1e.h5", star.kappa, star.n, star.εc_SI, TD ? star.ptol_TD : star.ptol, h);

function load(star::NeutronStarOscillations.Star, h::Float64; TD::Bool=false)
    h5f = h5open(fname(star, h; TD=TD), "r")
        r = h5f["radius"][:];
        m = h5f["mass"][:];
        p = h5f["pressure"][:];
        ε = h5f["energy_density"][:];
        ν = h5f["nu"][:];
    close(h5f)
    return r, m, p, ε, ν
end

@inline m_prime(r::Float64, p::Float64, ε::Float64, m::Float64)::Float64 = 4*π*r^2*ε
@inline nu0_prime(r::Float64, p::Float64, ε::Float64, m::Float64)::Float64 = (2*m + 8*π*r^3*p)/(r^2 - 2*r*m)
@inline p_prime(r::Float64, p::Float64, ε::Float64, m::Float64)::Float64 = -(((m + 4*π*r^3*p)*(p + ε))/(r*(r - 2*m)))
@inline n0_prime(r::Float64, p::Float64, ε::Float64, m::Float64, ε_prime::Float64, βn::Float64, βε::Float64)::Float64 = -((βε*ε_prime)/βn)

# EOS::Function is an equation of state ε = ε(p)
function solve(star::NeutronStarOscillations.Star, h::Float64; pressure_floor_fix::Bool=false, fixed_mass_dp::Bool=false, mass_tol::Float64=0.001, TD::Bool=false)
    if TD
        solve(star.EOS_ε, star.εc_SI, star.pc_SI, h, star.TOV_max_steps, star.ptol_TD; save_to_file=true, fname=Explicit.fname(star, h; TD=TD), pressure_floor_fix=pressure_floor_fix, fixed_mass_dp=fixed_mass_dp, mass_tol=mass_tol)
    else
        solve(star.EOS_ε, star.εc_SI, star.pc_SI, h, star.TOV_max_steps, star.ptol; save_to_file=true, fname=Explicit.fname(star, h; TD=TD), pressure_floor_fix=pressure_floor_fix, fixed_mass_dp=fixed_mass_dp, mass_tol=mass_tol)
    end
end

function solve(EOS::Function, εc::Float64, pc::Float64, h::Float64, max_steps::Int64, ptol::Float64; save_to_file::Bool=false, fname::String="", pressure_floor_fix::Bool=false, fixed_mass_dp::Bool=false, mass_tol::Float64=0.001)
    # initial conditions
    p = [pc];
    nu = [0.0];
    m = [0.0];
    ε = [εc];
    r = [0.0];

    # [p, m, nu]
    k1 = zeros(4);
    k2 = zeros(4);
    k3 = zeros(4);
    k4 = zeros(4);

    # parameters enforce regularity at the center at first time step. Their effect is removed after the first step.
    regularity = 0.0;
    shift = 1.0;

    # stop conditions 
    # non_zero_pressure = true;
    # num_steps = true;
    step_number = 1;

    try
        while p[end] > ptol && length(p) < max_steps
            r0 = r[end];
            p0 = p[end];
            m0 = m[end];
            ε0 = ε[end];
            nu0 = nu[end];

            # compute k1
            k1[1] = p_prime(r0 + shift, p0, ε0, m0) * regularity;
            k1[2] = m_prime(r0 + shift, p0, ε0, m0) * regularity;
            k1[3] = nu0_prime(r0 + shift, p0, ε0, m0) * regularity;

            # compute k2
            r1 = r0 + h/2;
            p1 = p0 + h/2*k1[1];
            m1 = m0 + h/2*k1[2];
            ε1 = EOS(p1);

            k2[1] = p_prime(r1, p1, ε1, m1);
            k2[2] = m_prime(r1, p1, ε1, m1);
            k2[3] = nu0_prime(r1, p1, ε1, m1);

            # compute k3
            r2 = r0 + h/2;
            p2 = p0 + h/2*k2[1];
            m2 = m0 + h/2*k2[2];
            ε2 = EOS(p2);

            k3[1] = p_prime(r2, p2, ε2, m2);
            k3[2] = m_prime(r2, p2, ε2, m2);
            k3[3] = nu0_prime(r2, p2, ε2, m2);

            # compute k4
            r3 = r0 + h;
            p3 = p0 + h*k3[1];
            m3 = m0 + h*k3[2];
            ε3 = EOS(p3);


            k4[1] = p_prime(r3, p3, ε3, m3);
            k4[2] = m_prime(r3, p3, ε3, m3);
            k4[3] = nu0_prime(r3, p3, ε3, m3);

            # update values
            push!(r, step_number * h);
            push!(p, p0 + h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]));
            push!(m, m0 + h/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]));
            push!(nu, nu0 + h/6*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]));
            push!(ε, EOS(p[end]));

            step_number += 1

            # remove regularity control after first time step
            regularity = 1.0;
            shift = 0.0;
        end
    catch e
        println("An error occurred during the integration. Stopping at p = $(p[end]) and r = $(r[end]).")
        rethrow(e)
    end

    if length(p) >= max_steps
        @warn "Warning: Maximum number of steps reached without reaching surface of the star."
    end

    # clean up final time step if pressure goes negative
    if p[end] < 0
        pop!(p)
        pop!(r)
        pop!(m)
        pop!(nu)
        pop!(ε)
    end

    if length(p) == max_steps
        @warn "Warning: Maximum number of steps reached without reaching surface of the star. Don't trust saved solution."
        if pressure_floor_fix
            @warn "Warning: Taking into account possible pressure floor in the solution where pressure cannot be driven arbitrarily close to zero. To disable, change pressure_floor_fix to false."
            pressure_floor = @. !isapprox(p[1:end], p[end], rtol=1e-8)
            p = p[pressure_floor];
            r = r[pressure_floor];
            m = m[pressure_floor];
            nu = nu[pressure_floor];
            ε = ε[pressure_floor];
        end
    end

    if fixed_mass_dp
        mask = @. (abs(m .- m[end]) > mass_tol)
        p = p[mask];
        r = r[mask];
        m = m[mask];
        nu = nu[mask];
        ε = ε[mask];
    end

    nu[end] # nu at surface of star
    nu_true = log(1 - 2 * m[end] / r[end]) # true value of nu at surface of star
    @. nu += -nu[end] + nu_true # adjust nu to match Schwarzschild solution at surface

    if save_to_file
        h5open(fname, "w") do file
            file["mass"] = m;
            file["pressure"] = p;
            file["energy_density"] = ε;
            file["nu"] = nu;
            file["radius"] = r;
            file["central_density"] = εc;
            file["ptol"] = ptol;
            file["h"] = h;
        end

    else
        return r, m, p, ε, nu
    end
end

end

# implicit Crank Nicholson solver (second-order accurate)
module Implicit
using HDF5
using LinearAlgebra
using Printf
using NeutronStarOscillations
# u0 = solution at current step
# u1 = solution at next step
# r0 = current radius
# h = step size
# f::Function = EOS in the form p = f(energy_density)
# u[1] = mass
# u[2] = pressure
# u[3] = energy density
# u[4] = nu

F1(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64)::Float64 = (-m0 + m1)/h - (π*(h + 2*r0)^2*(ε0 + ε1))/2.

F2(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64)::Float64 = (-2*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3))/((h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0)) + (-ν0 + ν1)/h

F3(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64)::Float64 = (-p0 + p1)/h + ((2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1))/(2.0*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0))

F4(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, f::Function)::Float64 = (p0 + p1)/2. + (-f(ε0) - f(ε1))/2.

A11(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (h*(h - 2*m0 - 2*m1 + 2*r0)*(h*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + (-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1)))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A12(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = 0

A13(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (h^2*π*(h + 2*r0)^3*(h - 2*m0 - 2*m1 + 2*r0)^2)/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A14(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (-2*h^2*π*(h + 2*r0)^3*(h - 2*m0 - 2*m1 + 2*r0)^2*(1/h + (2*m0 + 2*m1 + π*(h + 2*r0)^3*(2*p0 + 2*p1 + ε0 + ε1))/(2.0*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0))))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A21(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (-2*h^2*(h + 2*r0)*(1 + (p0 + p1)*π*(h + 2*r0)^2)*(p0 + p1 + ε0 + ε1)*df_dε(ε1))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A22(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = 0

A23(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (2*h*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0)^2*df_dε(ε1))/(h*((h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1)) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A24(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (2*h*((h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1)))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A31(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (-2*h^2*(h + 2*r0)*(1 + (p0 + p1)*π*(h + 2*r0)^2)*(p0 + p1 + ε0 + ε1))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A32(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = 0

A33(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (2*h*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0)^2)/(h*((h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1)) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A34(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (-4*h*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0)^2*(1/h + (2*m0 + 2*m1 + π*(h + 2*r0)^3*(2*p0 + 2*p1 + ε0 + ε1))/(2.0*(h + 2*r0)*(h - 2*m0 - 2*m1 + 2*r0))))/(h*((h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1)) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A41(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (4*h^2*(1 + (p0 + p1)*π*(h + 2*r0)^2)*(2*m0*(h - (h + 4*r0)*df_dε(ε1)) + 2*m1*(h - (h + 4*r0)*df_dε(ε1)) + (h + 2*r0)^2*(2*df_dε(ε1) + h*(p0 + p1)*π*(h + 2*r0)*(1 + df_dε(ε1)))))/((h - 2*m0 - 2*m1 + 2*r0)*(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1)))

A42(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = h

A43(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (4*h^2*π*(h + 2*r0)^3*(h + h*(p0 + p1)*π*(h + 2*r0)^2 + (h - 2*m0 - 2*m1 + 2*r0)*df_dε(ε1)))/(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1))

A44(m0::Float64, p0::Float64, ε0::Float64, ν0::Float64, m1::Float64, p1::Float64, ε1::Float64, ν1::Float64, r0::Float64, h::Float64, df_dε::Function)::Float64 = (-4*h^2*π*(h + 2*r0)^2*(4*m0^2 + 4*m1^2 - 4*m1*(h + 3*r0 + (p0 + p1)*π*r0*(h + 2*r0)^2) - 4*m0*(h - 2*m1 + 3*r0 + (p0 + p1)*π*r0*(h + 2*r0)^2) + (h + 2*r0)^2*(2 + (p0 + p1)*π*(h + 2*r0)*(2*(h + r0) + h*(p0 + p1)*π*(h + 2*r0)^2))))/((h - 2*m0 - 2*m1 + 2*r0)*(h*(h - 2*m0 - 2*m1 + 2*r0)*(2*m0 + 2*m1 + (p0 + p1)*π*(h + 2*r0)^3) + h^2*π*(h + 2*r0)^2*(h + 2*r0 + (p0 + p1)*π*(h + 2*r0)^3)*(p0 + p1 + ε0 + ε1) + (h - 2*m0 - 2*m1 + 2*r0)*(-2*m0*(h + 4*r0) - 2*m1*(h + 4*r0) + (h + 2*r0)^2*(2 + h*π*(h + 2*r0)*(2*p0 + 2*p1 + ε0 + ε1)))*df_dε(ε1)))

function compute_function_vector!(F::Vector{Float64}, u0::Vector{Float64}, u1::Vector{Float64}, r0::Float64, h::Float64, f::Function)
    F[1] = F1(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h);
    F[2] = F2(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h);
    F[3] = F3(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h);
    F[4] = F4(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, f);
end

function newton_iterate!(Δu::Vector{Float64}, F::Vector{Float64}, u0::Vector{Float64}, u1::Vector{Float64}, r0::Float64, h::Float64, df_dε::Function)
    Δu[1] = -(A11(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[1] +
                A12(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[2] +
                A13(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[3] +
                A14(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[4])
    Δu[2] = -(A21(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[1] +
                A22(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[2] +
                A23(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0 , h, df_dε) * F[3] +
                A24(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[4])
    Δu[3] = -(A31(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[1] +
                A32(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[2] +
                A33(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[3] +
                A34(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[4])
    Δu[4] = -(A41(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[1] +
                A42(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[2] +
                A43(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[3] +
                A44(u0[1], u0[2], u0[3], u0[4], u1[1], u1[2], u1[3], u1[4], r0, h, df_dε) * F[4])
end

fname(star::NeutronStarOscillations.Star, h::Float64)::String =
star.data_path * @sprintf("Implicit_TOV_κ_%d_n_%.1f_εc_%.2e_ptol_%.2e_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.ptol, h);

function load(star::NeutronStarOscillations.Star, h::Float64)
    h5f = h5open(fname(star, h), "r")
        r = h5f["radius"][:];
        m = h5f["mass"][:];
        p = h5f["pressure"][:];
        ε = h5f["energy_density"][:];
        ν = h5f["nu"][:];
    close(h5f)
    return r, m, p, ε, ν
end

function solve(star::NeutronStarOscillations.Star, h::Float64)
    solve(star.EOS_p, star.dp_dε, star.pc_SI, star.εc_SI, h, star.TOV_max_steps, star.ptol, star.TOV_iter_tol, star.TOV_max_iter; initial_r=star.TOV_initial_r, save_to_file=true, fname=Implicit.fname(star, h))
end

function solve(f::Function, df_dε::Function, pc::Float64, εc::Float64, h::Float64, max_steps::Int64, ptol::Float64, iter_tol::Float64, max_iter::Int64; initial_r::Float64=1e-15, save_to_file::Bool=false, fname::String = "")

    u0 = zeros(4); # initial conditions
    u1 = zeros(4); # next step solution
    Δu = zeros(4); # change in solution
    F = zeros(4);

    p = [pc];
    nu = [0.0];
    ε = [εc];
    r = [initial_r];
    m = [0.0]; # initial mass

    while p[end] > ptol && length(p) < max_steps
        n_iter = 0;
        r0 = r[end];
        u0[:] = [m[end], p[end], ε[end], nu[end]];
        u1[:] = u0;

        compute_function_vector!(F, u0, u1, r0, h, f)
        resid = norm(F);

        while abs(resid) > iter_tol && n_iter < max_iter
            newton_iterate!(Δu, F, u0, u1, r0, h, df_dε)
            u1 += Δu;
            compute_function_vector!(F, u0, u1, r0, h, f);
            resid = norm(F);
            n_iter += 1;
        end

        push!(m, u1[1]);
        push!(p, u1[2]);
        push!(ε, u1[3]);
        push!(nu, u1[4]);
        push!(r, r0 + h);
    end

    if length(p) == max_steps
        @warn "Warning: Maximum number of steps reached without reaching surface of the star. Don't trust saved solution."
    end


    nu[end] # nu at surface of star
    nu_true = log(1 - 2 * m[end] / r[end]) # true value of nu at surface of star
    @. nu += -nu[end] + nu_true # adjust nu to match Schwarzschild solution at surface

    if save_to_file
        h5open(fname, "w") do file
            file["mass"] = m;
            file["pressure"] = p;
            file["energy_density"] = ε;
            file["nu"] = nu;
            file["radius"] = r;
            file["central_density"] = εc;
            file["ptol"] = ptol;
            file["h"] = h;
        end

    else
        return r, m, p, ε, nu
    end
end

end

end
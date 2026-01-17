#=

    Module comprising of initial data solvers for the BDNK system. We assume stationary initial data for the energy, velocity and redshift perturbations. We then solve the
    lower-order EFEs, given a choice of the radial velocity perturbation, for the energy and redshift perturbations at time t=0.

=#

module BDNKInitialData
using ..QuickPlots
using ..TOV
using ..FiniteDiffOrder4
using LaTeXStrings
using Dierckx
using NeutronStarOscillations

# function plot_initial_data_convergence(r_lev1, var1_lev1, var2_lev1, r_lev2, var1_lev2, var2_lev2, r_lev3, var1_lev3, var2_lev3, var1_label, var2_label)
#     # check convergence
#     spline_order = 5; s = 0.0;
#     var1_spline_2 = Spline1D(r_lev2, var1_lev2; k=spline_order, s=s);
#     var2_spline_2 = Spline1D(r_lev2, var2_lev2; k=spline_order, s=s);
#     var1_spline_3 = Spline1D(r_lev3, var1_lev3; k=spline_order, s=s);
#     var2_spline_3 = Spline1D(r_lev3, var2_lev3; k=spline_order, s=s);

#     var2_lev2_ds = [var2_spline_2(r) for r in r_lev1];
#     var1_lev2_ds = [var1_spline_2(r) for r in r_lev1];
#     var2_lev3_ds = [var2_spline_3(r) for r in r_lev1];
#     var1_lev3_ds = [var1_spline_3(r) for r in r_lev1];

#     Q_var1 = [abs(var1_lev1[i] - var1_lev2_ds[i]) / abs(var1_lev2_ds[i] - var1_lev3_ds[i]) for i in 1:length(var1_lev1)];
#     Q_var2 = [abs(var2_lev1[i] - var2_lev2_ds[i]) / abs(var2_lev2_ds[i] - var2_lev3_ds[i]) for i in 1:length(var2_lev1)];

#     plot_ds = 100;
#     lim_x_min, lim_x_max, lim_y_min, lim_y_max = nothing, nothing, 0.0, 16.0;
#     colors = [:tomato, :aquamarine4, :turquoise]
#     linestyles = [:solid, :dash, :dot]
#     linewidths = [2.0, 2.0, 2.0]
#     alphas = [1.0, 1.0, 1.0]
#     xlabel = L"r\,[\mathrm{km}]"

#     labels = [var1_label, var2_label]
#     QuickPlots.plot11([r_lev1[1:plot_ds:end], r_lev1[1:plot_ds:end]], [Q_var1[1:plot_ds:end], Q_var2[1:plot_ds:end]], colors, labels, linestyles, linewidths, alphas, xlabel, L"Q_{N}", lim_x_min, lim_x_max, lim_y_min, lim_y_max, identity; legend = true)

#     lim_x_min, lim_x_max, lim_y_min, lim_y_max = nothing, nothing, nothing, nothing;
#     labels = ["Low", "Med", "High"]
#     QuickPlots.plot11([r_lev1[1:plot_ds:end], r_lev2[1:2*plot_ds:end], r_lev3[1:4*plot_ds:end]], [var1_lev1[1:plot_ds:end], var1_lev2[1:2*plot_ds:end], var1_lev3[1:4*plot_ds:end]], colors, labels, linestyles, linewidths, alphas, xlabel, var1_label, lim_x_min, lim_x_max, lim_y_min, lim_y_max, identity; legend = true, position = :lt)
#     QuickPlots.plot11([r_lev1[1:plot_ds:end], r_lev2[1:2*plot_ds:end], r_lev3[1:4*plot_ds:end]], [var2_lev1[1:plot_ds:end], var2_lev2[1:2*plot_ds:end], var2_lev3[1:4*plot_ds:end]], colors, labels, linestyles, linewidths, alphas, xlabel, var2_label, lim_x_min, lim_x_max, lim_y_min, lim_y_max, identity; legend = true, position = :lb)
# end

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)

# RK4 functions
du2_dr(u1::Float64, u2::Float64, u3::Float64, du1_dr::Float64, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (r*(-6*p*u1 - 6*u1*ε - (2*exp(ν/2.)*L*(3*ζ + 4*η)*τQ*u2*(2*r^2*(2*π*r*p*(1 + cs^2) + cs*cs_prime) + m*(1 + cs^2 - 4*r*cs*cs_prime)))/r^2 + (L*(3*ζ + 4*η)*τQ*(p + ε)*(r*(-3*exp(ν/2.)*((1 + 8*π*r^2*p)*u3 + 8*π*r^2*u2*cs^2) + 24*L*π*r*ζ*(p + ε)*(4*π*r^2*p*τP*u1 - (-1 + τP)*cs^2*(u1*(2 + 4*π*r^2*ε) + r*du1_dr)) + 32*L*π*r*η*(p + ε)*(4*π*r^2*p*τP*u1 - cs^2*(u1*(1 - 4*π*r^2*ε + τP*(2 + 4*π*r^2*ε)) + r*(-1 + τP)*du1_dr))) + 2*m*(3*exp(ν/2.)*((1 + 8*π*r^2*p)*u3 + 8*π*r^2*u2*cs^2) + 16*L*π*r*η*(p + ε)*(cs^2*(u1 - 2*r*du1_dr) + τP*(u1*(1 + 5*cs^2) + 2*r*cs^2*du1_dr)) + 12*L*π*r*ζ*(p + ε)*(-(cs^2*(5*u1 + 2*r*du1_dr)) + τP*(u1*(1 + 5*cs^2) + 2*r*cs^2*du1_dr)))))/(3.0*r*(r - 2*m))))/(2.0*exp(ν/2.)*L*(r - 2*m)*(3*ζ + 4*η)*τQ*cs^2)
du3_dr(u1::Float64, u2::Float64, u3::Float64, du1_dr::Float64, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (8*π*r^2*(r - 2*m)*u2 + (r - 2*m)*u3*(-1 + 8*π*r^2*ε) - (8*L*π*r*(3*ζ + 4*η)*τε*(p + ε)*(4*π*r^3*p*u1 - r*cs^2*(u1*(2 + 4*π*r^2*ε) + r*du1_dr) + m*(u1*(1 + 5*cs^2) + 2*r*cs^2*du1_dr)))/(3.0*exp(ν/2.)*cs^2))/(r - 2*m)^2

function compute_initial_data(star::NeutronStarOscillations.Star, h::Float64; u2_0::Float64=1e-6, return_all = false)
    NeutronStarOscillations.TOV.Explicit.solve(star, h; TD=true);
    r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.load(star, h; TD=true);

    TOV_length = length(m);
    cs = [sqrt(star.dp_dε(ε[i])) for i in 1:TOV_length];
    cs_prime = [cs_prime_func(m[i], p[i], ε[i], cs[i], r[i], star.d2p_dε2) for i in 1:TOV_length];
    cs_prime[1] = 0.0

    u1_arr = @. star.δu_ID(r);
    w1_arr = @. star.dδu_dr_ID(r);
    r_ds, u2, u3 = integrate(u1_arr, w1_arr, m, p, ε, ν, cs, cs_prime, r, u2_0, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L);
    
    if return_all
        ID_length = floor(Int, (TOV_length - 1) / 2) + 1;
        TOV_final_idx = 2 * (ID_length) - 1
        
        isapprox(r[1:2:TOV_final_idx], r_ds, atol=1e-12) || error("Radial grid points do not match between TOV solution and initial data solution.")
        
        return r_ds, u2, u3, u1_arr[1:2:TOV_final_idx], w1_arr[1:2:TOV_final_idx], m[1:2:TOV_final_idx], p[1:2:TOV_final_idx], ε[1:2:TOV_final_idx], ν[1:2:TOV_final_idx], cs[1:2:TOV_final_idx], cs_prime[1:2:TOV_final_idx]
    else
        return r_ds, u2, u3
    end
end

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function integrate(u1::Vector{Float64}, w1::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, u2_0::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)
    h = diff(r)[1];
    H = 2.0h;
    # initial conditions
    u3 = [0.0];
    u2 = [u2_0];
    rr = [r[1]];

    TOV_length = length(r); # minus four comes from removing the additional four steps in the background solution that allows us to evolve the pulsation equations until the surface of the star (since RK4 requires background solution at r + h/2 and r + h)
    ID_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [u3, u2]
    k1 = zeros(Float64, 2);
    k2 = zeros(Float64, 2);
    k3 = zeros(Float64, 2);
    k4 = zeros(Float64, 2);

    regularity = 1.0;

    for i in 1:ID_length-1
        TOV_sol_idx = 2 * i - 1
        
        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        cs_0 = cs[TOV_sol_idx];
        cs_0_prime = cs_prime[TOV_sol_idx];
        u1_0 = u1[TOV_sol_idx];
        w1_0 = w1[TOV_sol_idx];

        u3_0 = u3[end];
        u2_0 = u2[end];
 
        k1[1] = du3_dr(u1_0, u2_0, u3_0, w1_0, m0, p0, ε0, ν0, cs_0, cs_0_prime, r0 + regularity, η, ζ, τε, τP, τQ, L) * (1-regularity); # zero at center of star
        k1[2] = du2_dr(u1_0, u2_0, u3_0, w1_0, m0, p0, ε0, ν0, cs_0, cs_0_prime, r0 + regularity, η, ζ, τε, τP, τQ, L) * (1-regularity); # zero at center of star

        regularity = 0.0; # only nonzero at i = 1 to enforce regularity condition at r = 0

        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        cs_1 = cs[TOV_sol_idx+1];
        cs_1_prime = cs_prime[TOV_sol_idx+1];
        u1_1 = u1[TOV_sol_idx+1];
        w1_1 = w1[TOV_sol_idx+1];

        u3_1 = u3_0 + H/2*k1[1];
        u2_1 = u2_0 + H/2*k1[2];

        k2[1] = du3_dr(u1_1, u2_1, u3_1, w1_1, m1, p1, ε1, ν1, cs_1, cs_1_prime, r1, η, ζ, τε, τP, τQ, L);
        k2[2] = du2_dr(u1_1, u2_1, u3_1, w1_1, m1, p1, ε1, ν1, cs_1, cs_1_prime, r1, η, ζ, τε, τP, τQ, L);


        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ε2 = ε1;
        ν2 = ν1;
        cs_2 = cs_1;
        cs_2_prime = cs_1_prime;
        u1_2 = u1_1;
        w1_2 = w1_1;

        u3_2 = u3_0 + H/2*k2[1];
        u2_2 = u2_0 + H/2*k2[2];

        k3[1] = du3_dr(u1_2, u2_2, u3_2, w1_2, m2, p2, ε2, ν2, cs_2, cs_2_prime, r2, η, ζ, τε, τP, τQ, L);
        k3[2] = du2_dr(u1_2, u2_2, u3_2, w1_2, m2, p2, ε2, ν2, cs_2, cs_2_prime, r2, η, ζ, τε, τP, τQ, L);


        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        cs_3 = cs[TOV_sol_idx+2];
        cs_3_prime = cs_prime[TOV_sol_idx+2];
        u1_3 = u1[TOV_sol_idx+2];
        w1_3 = w1[TOV_sol_idx+2];

        u3_3 = u3_0 + H*k3[1];
        u2_3 = u2_0 + H*k3[2];

        k4[1] = du3_dr(u1_3, u2_3, u3_3, w1_3, m3, p3, ε3, ν3, cs_3, cs_3_prime, r3, η, ζ, τε, τP, τQ, L);
        k4[2] = du2_dr(u1_3, u2_3, u3_3, w1_3, m3, p3, ε3, ν3, cs_3, cs_3_prime, r3, η, ζ, τε, τP, τQ, L);


        # update values
        push!(rr, r0 + H);
        push!(u3, u3_0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]));
        push!(u2, u2_0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]));
    end
    # return rr, u3, u2, m[1:2:end], p[1:2:end], ε[1:2:end], ν[1:2:end], cs[1:2:end], cs_prime[1:2:end];
    return rr, u2, u3
end

end
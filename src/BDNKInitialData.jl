#=
    Module comprising of initial data solvers for the BDNK system. The Einstein-BDNK system consists of five equations: three second-order wave-like equations in the variables δu, δε, and δλ and two lower first-order constraint equations. The constraint equations reduce the number of variables for which we can freely specify initial data from six to four. The system we solve numerically is obtained by using one of the lower-order equations to eliminate time derivatives of δλ from all the equations, so in the final system we can sepcify initial data for only four of (δu, ∂_{t}δu, δε, ∂_{t}δε, δλ) since we have one remaining first-order constraint. We choose to specify intial data for δu, δε and their time derivatives, and use the constraint equation to solve for δλ.
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
du3_dr(u1::Float64, u2::Float64, u3::Float64, du1_dr::Float64, du2_dr::Float64, v1::Float64, v2::Float64, m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (L*π*r^4*(3*ζ + 4*η)*τε*(p + ε)*(-384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p*τQ*u2 - (48*exp(ν/2.)*L*π*(3*ζ + 4*η)*p*τQ*u3)/r - 384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p^2*τQ*u3 - (48*exp(ν/2.)*L*π*(3*ζ + 4*η)*τQ*u3*ε)/r - 384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p*τQ*u3*ε + (216*exp(ν/2.)*u2)/(L*r*(3*ζ + 4*η)*τε*(p + ε)) - (27*exp(ν/2.)*u3)/(L*π*r^3*(3*ζ + 4*η)*τε*(p + ε)) + (216*exp(ν/2.)*u3*ε)/(L*r*(3*ζ + 4*η)*τε*(p + ε)) - 768*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p*τQ*u2*cs^2 + (384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p*τP*τQ*u2*cs^2)/τε - (48*exp(ν/2.)*L*π*(3*ζ + 4*η)*p*τP*τQ*u3*cs^2)/(r*τε) - 384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*τQ*u2*ε*cs^2 + (384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*τP*τQ*u2*ε*cs^2)/τε - (48*exp(ν/2.)*L*π*(3*ζ + 4*η)*τP*τQ*u3*ε*cs^2)/(r*τε) + (384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*p*τP*τQ*u3*ε*cs^2)/τε + (384*exp(ν/2.)*L*π^2*r*(3*ζ + 4*η)*τP*τQ*u3*ε^2*cs^2)/τε - (1152*exp(ν/2.)*L*π^2*r*ζ*τQ*u2*(p + ε)*cs^2)/τε - (1536*exp(ν/2.)*L*π^2*r*η*τQ*u2*(p + ε)*cs^2)/τε + (144*exp(ν/2.)*L*π*ζ*τQ*u3*(p + ε)*cs^2)/(r*τε) + (192*exp(ν/2.)*L*π*η*τQ*u3*(p + ε)*cs^2)/(r*τε) - (1152*exp(ν/2.)*L*π^2*r*ζ*τQ*u3*ε*(p + ε)*cs^2)/τε - (1536*exp(ν/2.)*L*π^2*r*η*τQ*u3*ε*(p + ε)*cs^2)/τε + (24*u1*(6*r^2*cs^2 + (64*L^2*π^3*r^6*(3*ζ + 4*η)^2*p^2*τQ*(p + ε)*cs^2)/3. - 64*L^2*π^2*r^4*η*(3*ζ + 4*η)*τQ*ε*(p + ε)*cs^4 + 6*m^2*(1 + 5*cs^2) + 4*π*r^4*p*(-3 - 3*cs^2 + (16*L^2*π^2*r^2*(3*ζ + 4*η)^2*τQ*ε*(p + ε)*cs^2)/3. - 16*L^2*π*η*(3*ζ + 4*η)*τQ*(p + ε)*cs^4) + r*m*(-3*(1 + 9*cs^2) + (16*L^2*π^2*r^2*(3*ζ + 4*η)*τQ*ε*(p + ε)*cs^2*(3*ζ + 4*η*(1 + 6*cs^2)))/3. + 8*π*r^2*p*(2*L^2*π*ζ*(3*ζ + 4*η)*τQ*(p + ε)*cs^2 + 3*(1 + cs^2) + (8*L^2*π*η*(3*ζ + 4*η)*τQ*(p + ε)*cs^2*(1 + 6*cs^2))/3.))))/(r^3*(r - 2*m)*cs^2) - 192*exp(ν/2.)*L*π*(3*ζ + 4*η)*τQ*u2*cs*cs_prime + (72*du1_dr)/r - 96*exp(ν/2.)*L*π*(3*ζ + 4*η)*τQ*cs^2*du2_dr - (96*L*π*(3*ζ + 4*η)*p*τQ*v1)/exp(ν/2.) - (96*L*π*(3*ζ + 4*η)*τQ*ε*v1)/exp(ν/2.) + (72*v2)/(r*(p + ε)) - 384*L^2*π^2*r*ζ*(3*ζ + 4*η)*τQ*(p + ε)*cs^2*v2 - 512*L^2*π^2*r*η*(3*ζ + 4*η)*τQ*(p + ε)*cs^2*v2 - (18*m*((3*exp(ν/2.)*u3*(-1 + 8*π*r^2*ε))/(L*π*(3*ζ + 4*η)*τε*(p + ε)) + 8*exp(ν/2.)*r^2*u2*(3/(L*(3*ζ + 4*η)*τε*(p + ε)) + (2*L*π*(3*ζ + 4*η)*τQ*(1 + cs^2 - 4*r*cs*cs_prime))/3.) + 8*r^2*(du1_dr - (4*exp(ν/2.)*L*π*r*(3*ζ + 4*η)*τQ*cs^2*du2_dr)/3. + v2/(p + ε))))/r^4))/(3.0*exp(ν/2.)*(r - 2*m)*(9*r - 18*m + 16*L^2*π^2*r^3*(3*ζ + 4*η)^2*p^2*τP*τQ*cs^2 + 16*L^2*π^2*r^3*(3*ζ + 4*η)^2*τP*τQ*ε^2*cs^2 - 48*L^2*π^2*r^3*ζ*(3*ζ + 4*η)*τQ*ε*(p + ε)*cs^2 - 64*L^2*π^2*r^3*η*(3*ζ + 4*η)*τQ*ε*(p + ε)*cs^2 - 16*L^2*π^2*r^3*(3*ζ + 4*η)^2*p*τQ*(p + ε - 2*τP*ε)*cs^2))

function compute_initial_data(star::NeutronStarOscillations.Star, h::Float64; return_all = false)
    # compute TOV at high resolution
    h_TOV = 1e-4;
    r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV / 2.0; TD=true, save_to_file=false);

    TOV_length = length(m);
    cs = [sqrt(star.dp_dε(ε[i])) for i in 1:TOV_length];
    cs_prime = [cs_prime_func(m[i], p[i], ε[i], cs[i], r[i], star.d2p_dε2) for i in 1:TOV_length];
    cs_prime[1] = 0.0
    cs_prime_prime = zero(cs_prime)
    FiniteDiffOrder4.compute_first_derivative(cs_prime_prime, cs_prime, diff(r)[1], length(r));

    # downsample to user specified h divided by two (this is because the RK4 integration will solve at steps of 2h_TOV)
    ds_fact = argmin(@. abs(r - h/2)) - 1;
    r = r[1:ds_fact:end];
    diff(r)[1] ≈ h/2 ? nothing : error("Grid spacing does not match desired value of h after downsampling");
    m = m[1:ds_fact:end];
    p = p[1:ds_fact:end];
    ε = ε[1:ds_fact:end];
    ν = ν[1:ds_fact:end];
    cs = cs[1:ds_fact:end];
    cs_prime = cs_prime[1:ds_fact:end];
    cs_prime_prime = cs_prime_prime[1:ds_fact:end];
    TOV_length = length(m);

    # evaluate initial data functions
    u1_arr = @. star.δu_ID(r);
    u2_arr = @. star.δε_ID(r);
    w1_arr = @. star.δu_dr_ID(r);
    w2_arr = @. star.δε_dr_ID(r);
    v1_arr = @. star.δu_dt_ID(r);
    v2_arr = @. star.δε_dt_ID(r);
    r_ds, u3 = integrate(u1_arr, u2_arr, w1_arr, w2_arr, v1_arr, v2_arr, m, p, ε, ν, cs, cs_prime, r, star.η, star.ζ, star.τε, star.τP, star.τQ, star.L);
    
    if return_all
        ID_length = floor(Int, (TOV_length - 1) / 2) + 1;
        TOV_final_idx = 2 * (ID_length) - 1
        
        isapprox(r[1:2:TOV_final_idx], r_ds, atol=1e-12) || error("Radial grid points do not match between TOV solution and initial data solution.")
        
        return r_ds, u1_arr[1:2:TOV_final_idx], u2_arr[1:2:TOV_final_idx], u3, w1_arr[1:2:TOV_final_idx], w2_arr[1:2:TOV_final_idx], v1_arr[1:2:TOV_final_idx], v2_arr[1:2:TOV_final_idx], m[1:2:TOV_final_idx], p[1:2:TOV_final_idx], ε[1:2:TOV_final_idx], ν[1:2:TOV_final_idx], cs[1:2:TOV_final_idx], cs_prime[1:2:TOV_final_idx], cs_prime_prime[1:2:TOV_final_idx]
    else
        return r_ds, u3
    end
end

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function integrate(u1::Vector{Float64}, u2::Vector{Float64}, w1::Vector{Float64}, w2::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)
    h = diff(r)[1];
    H = 2.0h;

    # initial conditions
    u3 = [0.0];
    rr = [r[1]];

    TOV_length = length(r);
    ID_length = floor(Int, (TOV_length - 1) / 2) + 1;

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
        u2_0 = u2[TOV_sol_idx];
        w1_0 = w1[TOV_sol_idx];
        w2_0 = w2[TOV_sol_idx];
        v1_0 = v1[TOV_sol_idx];
        v2_0 = v2[TOV_sol_idx];

        u3_0 = u3[end];
 
        k1 = du3_dr(u1_0, u2_0, u3_0, w1_0, w2_0, v1_0, v2_0, m0, p0, ε0, ν0, cs_0, cs_0_prime, r0 + regularity, η, ζ, τε, τP, τQ, L) * (1-regularity); # zero at center of star
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
        u2_1 = u2[TOV_sol_idx+1];
        w1_1 = w1[TOV_sol_idx+1];
        w2_1 = w2[TOV_sol_idx+1];
        v1_1 = v1[TOV_sol_idx+1];
        v2_1 = v2[TOV_sol_idx+1];

        u3_1 = u3_0 + H/2*k1[1];
        k2 = du3_dr(u1_1, u2_1, u3_1, w1_1, w2_1, v1_1, v2_1, m1, p1, ε1, ν1, cs_1, cs_1_prime, r1, η, ζ, τε, τP, τQ, L);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ε2 = ε1;
        ν2 = ν1;
        cs_2 = cs_1;
        cs_2_prime = cs_1_prime;
        u1_2 = u1_1;
        u2_2 = u2_1;
        w1_2 = w1_1;
        w2_2 = w2_1;
        v1_2 = v1_1;
        v2_2 = v2_1;

        u3_2 = u3_0 + H/2*k2[1];

        k3 = du3_dr(u1_2, u2_2, u3_2, w1_2, w2_2, v1_2, v2_2, m2, p2, ε2, ν2, cs_2, cs_2_prime, r2, η, ζ, τε, τP, τQ, L);

        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        cs_3 = cs[TOV_sol_idx+2];
        cs_3_prime = cs_prime[TOV_sol_idx+2];
        u1_3 = u1[TOV_sol_idx+2];
        u2_3 = u2[TOV_sol_idx+2];
        w1_3 = w1[TOV_sol_idx+2];
        w2_3 = w2[TOV_sol_idx+2];
        v1_3 = v1[TOV_sol_idx+2];
        v2_3 = v2[TOV_sol_idx+2];

        u3_3 = u3_0 + H*k3[1];

        k4 = du3_dr(u1_3, u2_3, u3_3, w1_3, w2_3, v1_3, v2_3, m3, p3, ε3, ν3, cs_3, cs_3_prime, r3, η, ζ, τε, τP, τQ, L);

        # update values
        push!(rr, r0 + H);
        push!(u3, u3_0 + H/6*(k1 + 2*k2 + 2*k3 + k4));
    end
    # return rr, u3, u2, m[1:2:end], p[1:2:end], ε[1:2:end], ν[1:2:end], cs[1:2:end], cs_prime[1:2:end];
    return rr, u3
end

end
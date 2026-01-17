module NeutronStarOscillations
using LaTeXStrings
using Printf
using HDF5
using CairoMakie
using Dierckx
using PrecompileTools: @setup_workload, @compile_workload
export Star, get_star_type, solve_TOV, load_TOV, plot_TOV_var, compute_eigensystem, load_eigensystem, plot_fd_eigvecs, time_integrate, load_td_solution, plot_characteristic_speeds, plot_initial_data_convergence, plot_td_var, plot_IR_one_norms, animate_td_var, check_causal_frame

mutable struct Star
    ########## STELLAR PARAMETERS ##########
    εc_cgs::Float64 # central energy density [g/cm^3]
    εc_SI::Float64 # central energy density [km^(-2)]
    pc_SI::Float64 # central pressure [km^(-2)]
    # polytropic params for equation of state p = κ * ε^(1 + 1 / n)
    kappa::Float64 # polytropic prefactor [km^(-2 / n)]
    n::Float64 # polytropic index
    EOS_ε::Function # energy density ε(p) as a function of pressure
    EOS_p::Function # pressure p(ε) as a function of energy density
    dp_dε::Function # derivative of pressure with respect to energy density dp/dε
    d2p_dε2::Function # second derivative of pressure with respect to energy density
    # viscous parameters
    η::Float64 # dimensionless shear viscosity parameter 
    ζ::Float64 # dimensionless bulk viscosity parameter
    τε::Float64 # dimensionless relaxation time parameter
    τP::Float64 # dimensionless relaxation time parameter
    τQ::Float64 # dimensionless relaxation time parameter
    L::Float64 # length scale associated with viscous parameters [km]
    ########## NUMERICAL PARAMETERS ##########
    ptol::Float64 # pressure at which to terminate integration [km^(-2)] in TOV and frequency domain code
    ptol_TD::Float64 # pressure at which to terminate integration [km^(-2)] in time domain code
    data_path::String # path to save data files
    fig_path::String # path to save figure files
    # TOV-solver specific numerical parameters
    TOV_iter_tol::Float64 # tolerance for Newton iteration in implicit solver
    TOV_max_iter::Int64 # maximum number of iterations in implicit solver
    TOV_max_steps::Int64 # maximum number of steps to take in the solver
    TOV_initial_r::Float64 # initial radius for implicit solver [km]
    # frequency domain numerical parameters
    NL_reltol::Float64 # NonlinearSolve relative tolerance (used for iteration of shooting method integration to satisfy boundary conditions)
    NL_abstol::Float64 # NonlinearSolve absolute tolerance (used for iteration of shooting method integration to satisfy boundary conditions)
    NL_maxiter::Int64 # maximum number of NonlinearSolve iterations
    N_eigvals::Int64 # number of eigenvalue-eigenvector pairs to compute
    # time domain initial data
    ξ_ID::Function # initial data function for Lagrangian displacement (ξ) (used by perfect fluid / Eckart)
    δu_ID::Function # initial data function for perturbation of radial component of four-velocity (δu) (used by full BDNK)
    dδu_dr_ID::Function # initial data function for radial derivative of δu (used by full BDNK)
    # δε_ID::Function # initial data function for perturbation of energy density (used by Cowling BDNK)
    # dδε_dr_ID::Function # initial data function for radial derivative of δε (used by Cowling BDNK)
    # time domain numerical parameters
    KO::Float64 # Kreiss-Oliger dissipation factor
    CFL::Float64 # Courant-Friedrichs-Lewy factor
    dt_save::Float64 # time interval between saved data points [ms]
    h_save::Float64 # spatial grid spacing for saved data points [km]
    save_every::Int64 # save to file after 'save_every' time steps have been stored in memory (i.e., after every Δt = save_every * dt_save)
    T::Float64 # total integration time [ms]
    # checks to make sure that the parameters are valid
    Star(εc_cgs, εc_SI, pc_SI, kappa, n, EOS_ε, EOS_p, dp_dε, d2p_dε2, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T) = begin
    if εc_cgs <= 0.0 || εc_SI <= 0.0
        error("Central energy density must be > 0.")
    elseif pc_SI <= 0.0
        error("Central pressure must be > 0.")
    elseif kappa <= 0.0
        error("Polytropic prefactor 'kappa' must be > 0.")
    elseif n <= 0.0
        error("Polytropic index 'n' must be > 0.")
    elseif η < 0.0
        error("Shear viscosity parameter 'η' must be >= 0.")
    elseif ζ < 0.0
        error("Bulk viscosity parameter 'ζ' must be >= 0.")
    elseif τε < 0.0
        error("Relaxation time parameter 'τε' must be >= 0.")
    elseif τP < 0.0
        error("Relaxation time parameter 'τP' must be >= 0.")
    elseif τQ < 0.0
        error("Relaxation time parameter 'τQ' must be >= 0.")
    elseif (τε != 0.0 || τP != 0.0 || τQ != 0.0) && η == 0.0 && ζ == 0.0
        error("Relaxation time parameters 'τε', 'τP', and 'τQ' can only be non-zero if both of the viscosity parameters 'η' or 'ζ' are
         non-zero.")
    elseif L <= 0.0
        error("Length scale 'L' must be > 0.")
    elseif ptol <= 0.0
        error("Pressure tolerance 'ptol' must be > 0.")
    elseif ptol_TD <= 0.0
        error("Pressure tolerance for time domain solver 'ptol_TD' must be > 0.")
    elseif TOV_iter_tol <= 0.0
        error("TOV iteration tolerance 'TOV_iter_tol' must be > 0.")
    elseif TOV_max_iter <= 0
        error("TOV maximum iterations 'TOV_max_iter' must be > 0.") 
    elseif TOV_max_steps <= 0
        error("TOV maximum steps 'TOV_max_steps' must be > 0.")
    elseif TOV_initial_r <= 0.0
        error("TOV initial radius 'TOV_initial_r' must be > 0.")
    elseif NL_reltol <= 0.0
        error("NonlinearSolve relative tolerance 'NL_reltol' must be > 0.")
    elseif NL_abstol <= 0.0
        error("NonlinearSolve absolute tolerance 'NL_abstol' must be > 0.")
    elseif NL_maxiter <= 0
        error("NonlinearSolve maximum iterations 'NL_maxiter' must be > 0.")
    elseif N_eigvals <= 0
        error("Number of eigenvalues 'N_eigvals' must be > 0.")
    elseif !(0.0 <= KO <= 1.0)
        error("Kreiss-Oliger dissipation factor 'KO' must be between 0 and 1.")
    elseif CFL <= 0.0 || CFL > 1.0
        error("Courant-Friedrichs-Lewy factor 'CFL' must be between 0 and 1.")
    elseif !(0.0 <= dt_save <= T)
        error("Time interval between saved data points 'dt_save' must be between 0 and total integration time 'T'.")
    elseif h_save <= 0.0
        error("Spatial grid spacing for saved data points 'h_save' must be > 0.")
    elseif save_every <= 0
        error("Parameter 'save_every' must be a positive integer.")
    elseif T <= 0.0
        error("Total integration time 'T' must be > 0.")
    else
        new(εc_cgs, εc_SI, pc_SI, kappa, n, EOS_ε, EOS_p, dp_dε, d2p_dε2, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T)
    end
    end
end

# convenience constructor that specifies polytrope EoS functions
function Star(
    εc_cgs::Float64,
    kappa::Float64,
    n::Float64,
    η::Float64,
    ζ::Float64,
    τε::Float64,
    τP::Float64,
    τQ::Float64,
    L::Float64,
    ptol::Float64,
    ptol_TD::Float64,
    data_path::String,
    fig_path::String,
    TOV_iter_tol::Float64,
    TOV_max_iter::Int64,
    TOV_max_steps::Int64,
    TOV_initial_r::Float64,
    NL_reltol::Float64,
    NL_abstol::Float64,
    NL_maxiter::Int64,
    N_eigvals::Int64,
    ξ_ID::Function,
    δu_ID::Function,
    dδu_dr_ID::Function,
    # δε_ID::Function,
    # dδε_dr_ID::Function,
    KO::Float64,
    CFL::Float64,
    dt_save::Float64,
    h_save::Float64,
    save_every::Int64,
    T::Float64)::Star

    γ = 1.0 + 1.0/n;
    EOS_ε(p::Float64)::Float64 = (p / kappa)^(1/γ) 
    EOS_p(ε::Float64)::Float64 = kappa * ε^γ
    dp_dε(ε::Float64)::Float64 = γ * kappa * ε^(γ - 1)
    d2p_dε2(ε::Float64)::Float64 = γ * (γ - 1) * kappa * ε^(γ - 2);

    εc_SI = εc_cgs * gram_per_cm3_to_km_minus2;
    pc_SI = EOS_p(εc_SI);

    if ptol < 0
        ptol = abs(ptol) * pc_SI
    end

    if ptol_TD < 0
        ptol_TD = abs(ptol_TD) * pc_SI
    end

    return Star(εc_cgs, εc_SI, pc_SI, kappa, n, EOS_ε, EOS_p, dp_dε, d2p_dε2, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T)
end

# convenience constructor that uses Gaussian initial data for time domain simulations
function Star(
    εc_cgs::Float64,
    kappa::Float64,
    n::Float64,
    η::Float64,
    ζ::Float64,
    τε::Float64,
    τP::Float64,
    τQ::Float64,
    L::Float64,
    ptol::Float64,
    ptol_TD::Float64,
    data_path::String,
    fig_path::String,
    TOV_iter_tol::Float64,
    TOV_max_iter::Int64,
    TOV_max_steps::Int64,
    TOV_initial_r::Float64,
    NL_reltol::Float64,
    NL_abstol::Float64,
    NL_maxiter::Int64,
    N_eigvals::Int64,
    Gaussian_center::Float64,
    Gaussian_width::Float64,
    Gaussian_amplitude::Float64,
    KO::Float64,
    CFL::Float64,
    dt_save::Float64,
    h_save::Float64,
    save_every::Int64,
    T::Float64)::Star

    gaussian(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = A/exp((r - r0)^2/w^2)
    gaussian_prime(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = (2*A*(-r + r0))/(exp((r - r0)^2/w^2)*w^2)

    ξ_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width)
    δu_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width)
    dδu_dr_ID(r::Float64)::Float64 = gaussian_prime(r, Gaussian_amplitude, Gaussian_center, Gaussian_width)
    # δε_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width)
    # dδε_dr_ID(r::Float64)::Float64 = gaussian_prime(r, Gaussian_amplitude, Gaussian_center, Gaussian_width);

    return Star(εc_cgs, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T)
end


# convenience constructor that uses the perfect fluid eigenvector corresponding to the nth overtone as initial data for time domain simulations

# file name convention. Functions will, where necessary, prepend or append additional information to this base file name. For example, data will add ".h5", figures might add "_ξ.png" for a figure of the lagrangian displacement for an Eckart or PF Star, frequency domain code will add "FD_" prefix, etc.
include("main.jl")

# wrapper functions for main frequency domain and time domain functions
function get_star_type(star::NeutronStarOscillations.Star)::String
    if star.η == 0 && star.ζ == 0 && star.τε == 0 && star.τP == 0 && star.τQ == 0
        star_type = "PF"
    elseif (star.η > 0 || star.ζ > 0) && star.τε == 0 && star.τP == 0 && star.τQ == 0
        star_type = "Eckart"
    elseif star.η > 0 && star.ζ > 0 && star.τε > 0 && star.τP > 0 && star.τQ > 0
        star_type = "BDNK"
    else
        error("Invalid set of viscous parameters: η = $(star.η), ζ = $(star.ζ), τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ)")
    end

    return star_type
end

### FREQUENCY DOMAIN ###
function compute_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64, h_TOV::Float64, N_eigvals::Int64, nPointsMatrix::Int64, cowling::Bool; print_progress::Bool = true)
    star_type = get_star_type(star)
    
    if star_type == "BDNK"
        error("Frequency domain methods only for perfect and Eckart fluids but relaxation times specified as τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ).")
    end
    
    if cowling
        if star_type == "PF"
            print_progress ? println("Computing perfect fluid Cowling eigensystem...") : nothing
            NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.get_eigensystem(star, h_shoot, h_TOV, N_eigvals; nPointsMatrix = nPointsMatrix, print_progress = print_progress);
        else
            error("Cowling frequency domain methods only for perfect fluid but shear and bulk viscosity specified as η = $(star.η), ζ = $(star.ζ).")
        end
    else
        if star_type == "PF"
            print_progress ? println("Computing perfect fluid eigensystem...") : nothing
            NeutronStarOscillations.FrequencyDomain.PerfectFluid.Shoot.get_eigensystem(star, h_shoot, h_TOV, N_eigvals; nPointsMatrix = nPointsMatrix, print_progress = print_progress);
        else
            print_progress ? println("Computing Eckart fluid eigensystem...") : nothing
            NeutronStarOscillations.FrequencyDomain.Eckart.Shoot.get_eigensystem(star, h_shoot, h_TOV, N_eigvals; nPointsMatrix = nPointsMatrix, print_progress = print_progress);
        end
    end
end

function load_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64, cowling::Bool)
    star_type = get_star_type(star)
    
    if star_type == "BDNK"
        error("Frequency domain methods only for perfect and Eckart fluids but relaxation times specified as τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ).")
    end
    
    if cowling
        if star_type == "PF"
            return NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.load_eigensystem(star, h_shoot);
        else
            error("Cowling frequency domain methods only for perfect fluid but shear and bulk viscosity specified as η = $(star.η), ζ = $(star.ζ).")
        end
    else
        if star_type == "PF"
            return NeutronStarOscillations.FrequencyDomain.PerfectFluid.Shoot.load_eigensystem(star, h_shoot);
        else
            return NeutronStarOscillations.FrequencyDomain.Eckart.Shoot.load_eigensystem(star, h_shoot);
        end
    end
end

function plot_fd_eigvecs(modes::Vector{Int64}, star::NeutronStarOscillations.Star, method::String; N::Int64=0, h::Float64=0.0,
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

    NeutronStarOscillations.FrequencyDomain.plot_eigvecs(modes, star, method; N=N, h=h,
    desample_factor=desample_factor,
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
    scatter_lines = scatter_lines)
end


### TIME DOMAIN ###
function time_integrate(star::NeutronStarOscillations.Star, h::Float64, cowling::Bool; print_progress::Bool = true)
    star_type = get_star_type(star)
    if cowling
        if star_type == "PF"
            print_progress ? println("Integrating the perfect fluid cowling equations of motion in time...") : nothing
            NeutronStarOscillations.CowlingTimeDomain.PerfectFluid.solve(star, h; print_progress = print_progress);
        elseif star_type == "Eckart"
            error("Cowling time domain integration methods only for perfect and BDNK fluids but an Eckart star was specified.")
        else
            print_progress ? println("Integrating the BDNK cowling equations of motion in time...") : nothing
            NeutronStarOscillations.CowlingTimeDomain.BDNK.solve(star, h; print_progress = print_progress);
        end
    else
        if star_type == "PF"
            print_progress ? println("Integrating the perfect fluid equations of motion in time...") : nothing
            NeutronStarOscillations.TimeDomain.Eckart.solve(star, h; print_progress = print_progress); # PF time domain solver is Eckart with η = ζ = 0
        elseif star_type == "Eckart"
            print_progress ? println("Integrating the Eckart fluid equations of motion in time...") : nothing
            NeutronStarOscillations.TimeDomain.Eckart.solve(star, h; print_progress = print_progress);
        else
            print_progress ? println("Integrating the BDNK fluid equations of motion in time...") : nothing
            NeutronStarOscillations.TimeDomain.BDNK.solve(star, h; print_progress = print_progress);
        end
    end
end

function load_td_solution(star::NeutronStarOscillations.Star, h::Float64, cowling::Bool)
    if cowling
        return NeutronStarOscillations.CowlingTimeDomain.load_solution(star, h);
    else
        return NeutronStarOscillations.TimeDomain.load_solution(star, h);
    end
end

function plot_characteristic_speeds(star::NeutronStarOscillations.Star, h::Float64, cowling::Bool)
    star_type = get_star_type(star)
    if star_type == "BDNK"
        nothing
    else
        error("Characteristic speeds computed only for BDNK star but viscous parameters specified as η = $(star.η), ζ = $(star.ζ), τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ).")
    end

    if cowling
        NeutronStarOscillations.CowlingTimeDomain.plot_characteristic_speeds(star, h);
    else
        NeutronStarOscillations.TimeDomain.plot_characteristic_speeds(star, h);
    end 

end

function plot_initial_data_convergence(star::NeutronStarOscillations.Star, h::Float64, cowling::Bool)
    star_type = get_star_type(star)
    if star_type == "BDNK"
        nothing
    else
        error("Initial data solver only for BDNK star but viscous parameters specified as η = $(star.η), ζ = $(star.ζ), τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ).")
    end

    if cowling
        NeutronStarOscillations.CowlingTimeDomain.plot_initial_data_convergence(star, h);
    else
        NeutronStarOscillations.TimeDomain.plot_initial_data_convergence(star, h);
    end 

end

# plot single solution
function plot_td_var(star::NeutronStarOscillations.Star, h::Float64, var::String, time::Float64, cowling::Bool;
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

    if cowling
        NeutronStarOscillations.CowlingTimeDomain.plot_var(star, h, var, time;
        desample_factor=desample_factor,
        annotate_time=annotate_time,
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
        scatter_lines = scatter_lines)
    else
        NeutronStarOscillations.TimeDomain.plot_var(star, h, var, time;
        desample_factor=desample_factor,
        annotate_time=annotate_time,
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
        scatter_lines = scatter_lines)
    end
    

end

# plot one norm of independent residuals
function plot_IR_one_norms(star::NeutronStarOscillations.Star, h::Vector{Float64}, var::String, time::Float64, cowling::Bool;
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

    if cowling
        NeutronStarOscillations.CowlingTimeDomain.plot_IR_one_norms(star, h, var, time;
        desample_factor=desample_factor,
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
        scatter_lines = scatter_lines)
    else
        NeutronStarOscillations.TimeDomain.plot_IR_one_norms(star, h, var, time;
        desample_factor=desample_factor,
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
        scatter_lines = scatter_lines)
    end


end

# animate single solution
function animate_td_var(star::NeutronStarOscillations.Star, h::Float64, var::String, cowling::Bool;
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

    if cowling
        NeutronStarOscillations.CowlingTimeDomain.animate_var(star, h, var;
        desample_factor=desample_factor,
        stop_time=stop_time,
        animation_length=animation_length,
        framerate=framerate,
        colors = colors,
        labels = labels,
        linestyles = linestyles,
        linewidths = linewidths,
        alphas = alphas,
        xlabel = xlabel,
        ylabel = ylabel,
        lim_y_min = lim_y_min,
        lim_y_max = lim_y_max,
        fix_ylims = fix_ylims,
        yscale = yscale,
        text_xloc = text_xloc,
        text_yloc = text_yloc,
        text = text,
        xalign = xalign,
        yalign = yalign,
        legend = legend,
        position = position,
        width = width,
        height = height)
    else
        NeutronStarOscillations.TimeDomain.animate_var(star, h, var;
        desample_factor=desample_factor,
        stop_time=stop_time,
        animation_length=animation_length,
        framerate=framerate,
        colors = colors,
        labels = labels,
        linestyles = linestyles,
        linewidths = linewidths,
        alphas = alphas,
        xlabel = xlabel,
        ylabel = ylabel,
        lim_y_min = lim_y_min,
        lim_y_max = lim_y_max,
        fix_ylims = fix_ylims,
        yscale = yscale,
        text_xloc = text_xloc,
        text_yloc = text_yloc,
        text = text,
        xalign = xalign,
        yalign = yalign,
        legend = legend,
        position = position,
        width = width,
        height = height)
    end
end


## TOV equations ##
function solve_TOV(star::NeutronStarOscillations.Star, h::Float64)
    NeutronStarOscillations.TOV.Explicit.solve(star, h);
end

function load_TOV(star::NeutronStarOscillations.Star, h::Float64)
    return NeutronStarOscillations.TOV.Explicit.load(star, h);
end

function plot_TOV_var(var::String, star::NeutronStarOscillations.Star, h::Float64; 
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

    type = "Explicit"
    NeutronStarOscillations.TOV.plot_var(var, star, h, type; 
        desample_factor=desample_factor,
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
        scatter_lines = scatter_lines)
end

# convenience constructor that uses frequency domain eigenvector corresponding to the nth overtone as initial data for time domain simulations. Warning: this function computes eigenvectors using the frequency domain code and then interpolates them to create initial data functions. If the time domain grid extends slightly beyond the frequency domain grid, an error will be thrown since the interpolation functions are not defined outside the frequency domain grid. To avoid this, make sure the frequency domain code is ran at a higher resolution than the time domain code with the same (or smaller) pressure tolerance.
function Star(
    εc_cgs::Float64,
    kappa::Float64,
    n::Float64,
    η::Float64,
    ζ::Float64,
    τε::Float64,
    τP::Float64,
    τQ::Float64,
    L::Float64,
    ptol::Float64,
    ptol_TD::Float64,
    data_path::String,
    fig_path::String,
    TOV_iter_tol::Float64,
    TOV_max_iter::Int64,
    TOV_max_steps::Int64,
    TOV_initial_r::Float64,
    NL_reltol::Float64,
    NL_abstol::Float64,
    NL_maxiter::Int64,
    N_eigvals::Int64,
    mode::Int64,
    cowling::Bool,
    KO::Float64,
    CFL::Float64,
    dt_save::Float64,
    h_save::Float64,
    save_every::Int64,
    T::Float64; nPointsMatrix::Int64 = 500, h_shoot::Float64 = 1e-2, h_TOV::Float64 = 1e-5)::Star

    gaussian(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = A/exp((r - r0)^2/w^2)
    gaussian_prime(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = (2*A*(-r + r0))/(exp((r - r0)^2/w^2)*w^2)

    # define interim values to create a star object
    ξ_ID(r::Float64)::Float64 = 0.0
    δu_ID(r::Float64)::Float64 = 0.0
    dδu_dr_ID(r::Float64)::Float64 = 0.0
    star = Star(εc_cgs, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T)

    # create a new star with relaxation times set to zero (since we only have frequency domain methods for PF and Eckart stars)
    star_FD = Star(εc_cgs, kappa, n, η, ζ, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, ξ_ID, δu_ID, dδu_dr_ID, KO, CFL, dt_save, h_save, save_every, T);
    star_FD.τε = 0.0
    star_FD.τP = 0.0
    star_FD.τQ = 0.0

    # compute and load eigensystem
    N_eigvals = mode + 1; # compute up to and including desired mode
    compute_eigensystem(star_FD, h_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = false)
    init_freqs, shooting_freqs, r, evecs, resids, ret = load_eigensystem(star_FD, h_shoot, cowling)

    # lagrangian displacement initial data
    y1 = real.(evecs[:, mode + 1][:])

    # velocity perturbation initial data, since δu = ∂_{t} ξ = -iω * ξ
    y2 = real.(-im * shooting_freqs[mode + 1] * evecs[:, mode + 1][:])

    # interpolate
    ξ_spline = Spline1D(r, y1; k=NeutronStarOscillations.FrequencyDomain.spline_order, s=NeutronStarOscillations.FrequencyDomain.s, bc=NeutronStarOscillations.FrequencyDomain.bc)
    δu_spline = Spline1D(r, y2; k=NeutronStarOscillations.FrequencyDomain.spline_order, s=NeutronStarOscillations.FrequencyDomain.s, bc=NeutronStarOscillations.FrequencyDomain.bc)
    
    ξ_func(r::Float64)::Float64 = ξ_spline(r);
    δu_func(r::Float64)::Float64 = δu_spline(r);
    dδu_dr_func(r::Float64)::Float64 = derivative(δu_spline, r);

    star.ξ_ID = ξ_func
    star.δu_ID = δu_func
    star.dδu_dr_ID = dδu_dr_func
    return star
end

# upper bound on sound speed squared for causality in BDNK theory with τP = 1
cs2_max1(star::NeutronStarOscillations.Star) = (star.τQ*star.τε)/(star.τQ + star.τε + star.τQ*star.τε)

# upper bound on sound speed squared for causality in BDNK theory with τP > 1
cs2_max2(star::NeutronStarOscillations.Star) = (star.τε + star.τQ*(star.τP + star.τε) - sqrt(-4*(-1 + star.τP)*star.τQ^2*star.τε + (star.τε + star.τQ*(star.τP + star.τε))^2))/(2.0*(-1 + star.τP)*star.τQ)

function check_causal_frame(star::NeutronStarOscillations.Star)
    star_type = get_star_type(star)
    if star_type != "BDNK"
        error("Causality check only for BDNK star but viscous parameters specified as η = $(star.η), ζ = $(star.ζ), τε = $(star.τε), τP = $(star.τP), τQ = $(star.τQ).")
    end

    if star.τP < 1
        @warn("Causality violated: τP ≥ 1 for causality but τP = $(star.τP) < 1.")
    elseif star.τP == 1.0
        star.dp_dε(star.εc_SI) <= cs2_max1(star) ? println("Causal choice of frame.") : @warn("Causality violated: cs²(0) = $(star.dp_dε(star.εc_SI)) > cs²_max = $(cs2_max1(star)).")
    else
        star.dp_dε(star.εc_SI) <= cs2_max2(star) ? println("Causal choice of frame.") : @warn("Causality violated: cs²(0) = $(star.dp_dε(star.εc_SI)) > cs²_max = $(cs2_max2(star)).")
    end
end

# precompilation
@setup_workload begin
    # --- START OF SILENCER ---
    # Define the "Lying" Null Display locally
    struct SilentDisplay <: AbstractDisplay end
    Base.displayable(d::SilentDisplay, ::MIME) = true
    Base.display(d::SilentDisplay, x) = nothing
    
    # Push the display and force Makie to look for it
    Base.pushdisplay(SilentDisplay())
    # Save current inline state to restore it later
    old_inline = Makie.inline!()
    Makie.inline!(true)
    # --- END OF SILENCER ---
    @compile_workload begin
        include("./Examples/params.jl")

        # adjust star parameters to make runs faster
        h_save = 1e-1;
        total_time_ms = 0.0005;
        dt_save_ms = 0.0005;
        save_every = 1;

        PF_star.T = total_time_ms;
        PF_star.dt_save = dt_save_ms;
        PF_star.save_every = save_every;
        PF_star.h_save = h_save;
        Eckart_star.T = total_time_ms;
        Eckart_star.dt_save = dt_save_ms;
        Eckart_star.save_every = save_every;
        Eckart_star.h_save = h_save;
        BDNK_star.T = total_time_ms;
        BDNK_star.dt_save = dt_save_ms;
        BDNK_star.save_every = save_every;
        BDNK_star.h_save = h_save;

        # TOV solver
        h = 1e-4
        NeutronStarOscillations.solve_TOV(BDNK_star, h); # explicit RK4 solver

        var = "mass";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"m / M_{\odot}";
        legend = false;
        plot_desample_factor = 100; # desampling factor for plots to reduce number of points plotted
        NeutronStarOscillations.plot_TOV_var(var, BDNK_star, h; 
            desample_factor = plot_desample_factor,
            xlabel = xlabel,
            ylabel = ylabel,
            legend = legend,
        );

        ## frequency domain ##
        nPointsMatrix = 500;
        h1_shoot = 1e-1;
        h_TOV = 1e-4;

        # PF cowling
        cowling = true;
        NeutronStarOscillations.compute_eigensystem(PF_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = false);
        modes = [0, 1, 2, 3];
        method = "Shoot Cowling";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        legend = true;
        labels = [L"n=%$mode" for mode in modes];
        framevisible = true;

        NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
            modes, PF_star, method; h=h1_shoot,
            xlabel = xlabel,
            ylabel = ylabel,
            legend = legend,
            labels = labels,
            position = :lb,
            framevisible = framevisible
        );

        # PF
        cowling = false;
        NeutronStarOscillations.compute_eigensystem(PF_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = false);
        modes = [0, 1, 2, 3];
        method = "Shoot";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        legend = true;
        labels = [L"n=%$mode" for mode in modes];
        framevisible = true;

        NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
            modes, PF_star, method; h=h1_shoot,
            xlabel = xlabel,
            ylabel = ylabel,
            legend = legend,
            labels = labels,
            position = :lb,
            framevisible = framevisible
        );

        # Eckart
        cowling = false;
        NeutronStarOscillations.compute_eigensystem(Eckart_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = false);
        modes = [0, 1, 2, 3];
        method = "Shoot";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        legend = true;
        labels = [L"n=%$mode" for mode in modes];
        framevisible = true;

        NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
            modes, Eckart_star, method; h=h1_shoot,
            xlabel = xlabel,
            ylabel = ylabel,
            legend = legend,
            labels = labels,
            position = :lb,
            framevisible = framevisible
        );

        ## time domain ##
        # PF Cowling time domain
        h1 = h_save;
        cowling = true;
        NeutronStarOscillations.time_integrate(PF_star, h1, cowling; print_progress = false);
        time = 0.0;
        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";

        NeutronStarOscillations.plot_td_var(
            PF_star, h1, var, 0.0, cowling;
            desample_factor=1,
            annotate_time=true,
            xlabel = xlabel,
            ylabel = ylabel
        );

        resolutions = [h1];
        var = "IR";
        xlabel = L"t\,[\mathrm{ms}]";
        ylabel = L"||\mathrm{IR}||_1";
        labels = ["Lev1", "Lev2", "Lev3"];
        legend = true;
        position = :lt;
        NeutronStarOscillations.plot_IR_one_norms(
            PF_star, resolutions, var, time, cowling;
            desample_factor=1,
            xlabel = xlabel,
            ylabel = ylabel,
            yscale = identity,
            labels = labels,
            legend = legend,
            position = position
        );

        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        NeutronStarOscillations.animate_td_var(
            PF_star, h1, var, cowling;
            desample_factor=1,
            stop_time=PF_star.T, # stop time of animation in ms
            animation_length=10.0, # length of animation in seconds
            framerate=10, # frames per second
            xlabel = xlabel,
            ylabel = ylabel,
            fix_ylims = false, # whether to fix y-limits of plot
        )

        # PF time domain
        cowling = false;
        NeutronStarOscillations.time_integrate(PF_star, h1, cowling; print_progress = false);
        time = 0.0;
        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";

        NeutronStarOscillations.plot_td_var(
            PF_star, h1, var, 0.0, cowling;
            desample_factor=1,
            annotate_time=true,
            xlabel = xlabel,
            ylabel = ylabel
        );

        resolutions = [h1];
        var = "IR";
        xlabel = L"t\,[\mathrm{ms}]";
        ylabel = L"||\mathrm{IR}||_1";
        labels = ["Lev1", "Lev2", "Lev3"];
        legend = true;
        position = :lt;
        NeutronStarOscillations.plot_IR_one_norms(
            PF_star, resolutions, var, time, cowling;
            desample_factor=1,
            xlabel = xlabel,
            ylabel = ylabel,
            yscale = identity,
            labels = labels,
            legend = legend,
            position = position
        );

        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        NeutronStarOscillations.animate_td_var(
            PF_star, h1, var, cowling;
            desample_factor=1,
            stop_time=PF_star.T, # stop time of animation in ms
            animation_length=10.0, # length of animation in seconds
            framerate=10, # frames per second
            xlabel = xlabel,
            ylabel = ylabel,
            fix_ylims = false, # whether to fix y-limits of plot
        )

        # Eckart time domain
        cowling = false;
        NeutronStarOscillations.time_integrate(Eckart_star, h1, cowling; print_progress = false);
        time = 0.0;
        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";

        NeutronStarOscillations.plot_td_var(
            Eckart_star, h1, var, 0.0, cowling;
            desample_factor=1,
            annotate_time=true,
            xlabel = xlabel,
            ylabel = ylabel
        );

        resolutions = [h1];
        var = "IR";
        xlabel = L"t\,[\mathrm{ms}]";
        ylabel = L"||\mathrm{IR}||_1";
        labels = ["Lev1", "Lev2", "Lev3"];
        legend = true;
        position = :lt;
        NeutronStarOscillations.plot_IR_one_norms(
            Eckart_star, resolutions, var, time, cowling;
            desample_factor=1,
            xlabel = xlabel,
            ylabel = ylabel,
            yscale = identity,
            labels = labels,
            legend = legend,
            position = position
        );

        var = "xi";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\xi\,[\mathrm{km}]";
        NeutronStarOscillations.animate_td_var(
            Eckart_star, h1, var, cowling;
            desample_factor=1,
            stop_time=Eckart_star.T, # stop time of animation in ms
            animation_length=10.0, # length of animation in seconds
            framerate=10, # frames per second
            xlabel = xlabel,
            ylabel = ylabel,
            fix_ylims = false, # whether to fix y-limits of plot
        )
        
        # BDNK cowling time domain
        cowling = true;
        NeutronStarOscillations.time_integrate(BDNK_star, h1, cowling; print_progress = false);
        time = 0.0;
        var = "du";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\delta{u}";

        NeutronStarOscillations.plot_td_var(
            BDNK_star, h1, var, 0.0, cowling;
            desample_factor=1,
            annotate_time=true,
            xlabel = xlabel,
            ylabel = ylabel
        );

        resolutions = [h1];
        var = "IR1";
        xlabel = L"t\,[\mathrm{ms}]";
        ylabel = L"||\mathrm{IR}||_1";
        labels = ["Lev1", "Lev2", "Lev3"];
        legend = true;
        position = :lt;
        NeutronStarOscillations.plot_IR_one_norms(
            BDNK_star, resolutions, var, time, cowling;
            desample_factor=1,
            xlabel = xlabel,
            ylabel = ylabel,
            yscale = identity,
            labels = labels,
            legend = legend,
            position = position
        );

        var = "du";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\delta{u}";
        NeutronStarOscillations.animate_td_var(
            BDNK_star, h1, var, cowling;
            desample_factor=1,
            stop_time=BDNK_star.T, # stop time of animation in ms
            animation_length=10.0, # length of animation in seconds
            framerate=10, # frames per second
            xlabel = xlabel,
            ylabel = ylabel,
            fix_ylims = false, # whether to fix y-limits of plot
        )

        # BDNK time domain
        cowling = false;
        NeutronStarOscillations.time_integrate(BDNK_star, h1, cowling; print_progress = false);
        time = 0.0;
        var = "du";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\delta{u}";

        NeutronStarOscillations.plot_td_var(
            BDNK_star, h1, var, 0.0, cowling;
            desample_factor=1,
            annotate_time=true,
            xlabel = xlabel,
            ylabel = ylabel
        );

        resolutions = [h1];
        var = "Constrained_IR1";
        xlabel = L"t\,[\mathrm{ms}]";
        ylabel = L"||\mathrm{IR}||_1";
        labels = ["Lev1", "Lev2", "Lev3"];
        legend = true;
        position = :lt;
        NeutronStarOscillations.plot_IR_one_norms(
            BDNK_star, resolutions, var, time, cowling;
            desample_factor=1,
            xlabel = xlabel,
            ylabel = ylabel,
            yscale = identity,
            labels = labels,
            legend = legend,
            position = position
        );

        var = "du";
        xlabel = L"r\,[\mathrm{km}]";
        ylabel = L"\delta{u}";
        NeutronStarOscillations.animate_td_var(
            BDNK_star, h1, var, cowling;
            desample_factor=1,
            stop_time=BDNK_star.T, # stop time of animation in ms
            animation_length=10.0, # length of animation in seconds
            framerate=10, # frames per second
            xlabel = xlabel,
            ylabel = ylabel,
            fix_ylims = false, # whether to fix y-limits of plot
        )        
        
    end
end


end

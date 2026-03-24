#=

    Parameter file for example code. Default parameters correspond to a viscous BDNK star, but the appropriate viscous parameters can be set to zero to obtain Eckart and perfect fluid stars.

=#

using Revise
using NeutronStarOscillations
using LaTeXStrings

#################### STELLAR PARAMETERS ####################
eps_central = 5.5e15; # central energy density [g/cm^3] — (Float64)

# polytropic parameter for equation of state p = κ * ε^(1 + 1 / n)
n = 1.0; # polytropic index — (Float64)
kappa = 100.0; # polytropic prefactor [km^(-2 / n)] — (Float64)

# viscous parameters
η = 1.0e-2; # dimensionless shear viscosity parameter — (Float64) 
ζ = 1.0e-2; # dimensionless bulk viscosity parameter — (Float64)
τε = 15.0; # dimensionless relaxation time parameter — (Float64)
τP = 1.5; # dimensionless relaxation time parameter — (Float64)
τQ = 20.0; # dimensionless relaxation time parameter — (Float64)
L = 1.0; # length scale associated with viscous parameters [km] — (Float64)

#################### NUMERICAL PARAMETERS ####################
ptol = -1e-6; # pressure at which to terminate integration in TOV and frequency domain code [km^(-2)]. Positive value: sets pressure tolerance directly. Negative value: sets pressure tolerance as ptol * p(r=0) — (Float64)
ptol_TD = -1e-3; # pressure at which to terminate integration for time domain solver [km^(-2)]. Positive value: sets pressure tolerance directly. Negative value: sets pressure tolerance as ptol_TD * p(r=0) — (Float64)

# TOV-solver specific numerical parameters for implicit solver
TOV_iter_tol = 1e-15; # tolerance for Newton iteration in implicit solver — (Float64)
TOV_max_iter = 10; # maximum number of iterations in implicit solver — (Int64)
TOV_max_steps = Int(1e11); # maximum number of steps to take in the solver — (Int64)
TOV_initial_r = 1e-15; # initial radius for implicit solver [km] — (Float64)

# frequency domain numerical parameters
NL_reltol = 1e-12; # NonlinearSolve relative tolerance (used for iteration of shooting method integration to satisfy boundary conditions) — (Float64)
NL_abstol = 1e-12; # NonlinearSolve absolute tolerance (used for iteration of shooting method integration to satisfy boundary conditions) — (Float64)
NL_maxiter = 200; # maximum number of NonlinearSolve iterations — (Int64)
N_eigvals = 5; # number of eigenvalue-eigenvector pairs to compute — (Int64)

#= 
    Time domain initial data functions. The example scripts show different options for specifying initial data. For the perfeclt fluid / Eckart equations, the lagrangian displacement ξ and its time derivative must be specified at the initial time. For the BDNK fluid, one must specidy the radial velocity perturbation δu, the energy density perturbation δε, and their radial and time derivatives. Inside the BDNK time domain solver, this initial data is used to solve an ODE for the initial data in the metric perturbation δλ. We additionally provide in built functions for a Gaussian in δu and δε (so one just specifies the amplitude, center, and width) as well as an option for using the perfect fluid / Eckart eigenvectors as initial data for ξ and/or δu.
=#

Gaussian_amplitude = 1.0; # amplitude of Gaussian initial data — (Float64)
Gaussian_center = 5.0; # center of Gaussian initial data [km] — (Float64)
Gaussian_width = 0.5; # width of Gaussian initial data [km] — (Float64)
gaussian(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = A/exp((r - r0)^2/w^2)
gaussian_prime(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = (2*A*(-r + r0))/(exp((r - r0)^2/w^2)*w^2)

xi_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for Lagrangian displacement (ξ) (used by perfect fluid / Eckart) — (Function)
xi_dt_ID(r::Float64)::Float64 = 0.0 # initial data function for time derivative of ξ (used by perfect fluid / Eckart) — (Function)

du_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for perturbation of radial component of four-velocity (δu) (used by BDNK) — (Function)
du_dr_ID(r::Float64)::Float64 = gaussian_prime(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for radial derivative of δu (used by BDNK) — (Function)
du_dt_ID(r::Float64)::Float64 = 0.0 # initial data function for time derivative of δu (used by BDNK) — (Function)

de_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for perturbation of energy density (δε) (used by BDNK) — (Function)
de_dr_ID(r::Float64)::Float64 = gaussian_prime(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for radial derivative of δε (used by BDNK) — (Function)
de_dt_ID(r::Float64)::Float64 = 0.0 # initial data function for time derivative of δε (used by BDNK) — (Function)

# time domain numerical parameters
KO = 0.2; # Kreiss-Oliger dissipation coefficient for BDNK evolution — (Float64)
CFL = 0.1; # Courant-Friedrichs-Lewy factor — (Float64)
h_save = 4e-2; # spatial grid spacing for saved data points [km] — (Float64)
total_time_ms = 0.05; # total integration time [ms] — (Float64)
dt_save_ms = total_time_ms / 100.0; # time interval between saved data points [ms] — (Float64)
save_every = 50; # save to file after 'save_every' time steps have been stored in memory (i.e., after every Δt = save_every * dt_save_ms) — (Int64)

# file name conventions
data_path = "./Results/Data/";
fig_path  = "./Results/Figures/";
mkpath(data_path)
mkpath(fig_path)

#################### CREATING STAR OBJECTS ####################
# set all viscous parameters to zero for perfect fluid star as well as KO since only used in time integration of BDNK equations
PF_star = NeutronStarOscillations.Star(eps_central, kappa, n, 0.0, 0.0, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, xi_dt_ID, du_ID, du_dr_ID, du_dt_ID, de_ID, de_dr_ID, de_dt_ID, 0.0, CFL, dt_save_ms, h_save, save_every, total_time_ms);

# set relaxation times to zero for Eckart star as well as KO since only used in time integration of BDNK equations
Eckart_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, xi_dt_ID, du_ID, du_dr_ID, du_dt_ID, de_ID, de_dr_ID, de_dt_ID, 0.0, CFL, dt_save_ms, h_save, save_every, total_time_ms);

# BDNK star
BDNK_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, xi_dt_ID, du_ID, du_dr_ID, du_dt_ID, de_ID, de_dr_ID, de_dt_ID, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

#################### DEFAULT INITIAL DATA OPTIONS ####################
#=
    There are two default options for the initial data in ξ and δu. 
        
        
        (1) Gaussian initial data, like that specified explicitly above. One must set 'Gaussian_center', 'Gaussian_width', 'Gaussian_amplitude' to desired values.

            gaussian_ID_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, Gaussian_center, Gaussian_width, Gaussian_amplitude, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

        
        (2) Use eigenvectors from the perfect or Eckart fluid frequency domain solver. One must set 'mode' to the desired mode number (0 for fundamental, 1 for first overtone, etc.) and 'cowling' to true or false depending on whether one wants to use the Cowling approximation or not. Note that 'cowling=true' will only work for the perfect fluid case since there are no functions for an Eckart fluid in the Cowling approximation.

            mode = 0; cowling = false;
            eigvec_ID_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, mode, cowling, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

=#
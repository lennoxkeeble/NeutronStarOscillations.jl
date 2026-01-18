#=

    Parameter file for all example code. Default parameters correspond to a viscous BDNK star, but at the end we create objects corresponding to perfect fluid, Eckart, and
    BDNK stars which are used throughout the example code. To make it easier to run convergence tests, we omit parameters in these objects that change with resolution (e.g.,
    the number of grid points in the matrix method frequency domain solver or the grid spacing in the shooting method and time domain solvers). These are instead specified
    in the example scripts.

=#

using Revise
using NeutronStarOscillations
using LaTeXStrings

#################### STELLAR PARAMETERS ####################
eps_central = 3.0e15; # central energy density [g/cm^3] — (Float64)

# polytropic parameter for equation of state p = κ * ε^(1 + 1 / n)
n = 0.8; # polytropic index — (Float64)
kappa = 700.0; # polytropic prefactor [km^(-2 / n)] — (Float64)

# viscous parameters
η = 1.0e-4; # dimensionless shear viscosity parameter — (Float64) 
ζ = 1.0e-4; # dimensionless bulk viscosity parameter — (Float64)
τε = 10.0; # dimensionless relaxation time parameter — (Float64)
τP = 1.0; # dimensionless relaxation time parameter — (Float64)
τQ = 5.0; # dimensionless relaxation time parameter — (Float64)
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
NL_reltol = 1e-6; # NonlinearSolve relative tolerance (used for iteration of shooting method integration to satisfy boundary conditions) — (Float64)
NL_abstol = 1e-6; # NonlinearSolve absolute tolerance (used for iteration of shooting method integration to satisfy boundary conditions) — (Float64)
NL_maxiter = 200; # maximum number of NonlinearSolve iterations — (Int64)
N_eigvals = 5; # number of eigenvalue-eigenvector pairs to compute — (Int64)

#= 
    time domain initial data functions. In the Eckart case, u1 is the Lagrangian displacement. In the BDNK case, u1 is the radial velocty perturbation.
    In both cases we assume stationary initial data (i.e., time derivatives of perturbations are zero at t=0). The example scripts show different options for
    specifying initial data. Below is the most basic version where one directly specifies δu and its radial derivative. We provide in built functions for
    the Gaussian (so one just specifies the amplitude, center, and width) as well as an option for using the perfect fluid eigenvectors as initial data.
=#

Gaussian_amplitude = 1.0; # amplitude of Gaussian initial data — (Float64)
Gaussian_center = 5.0; # center of Gaussian initial data [km] — (Float64)
Gaussian_width = 0.5; # width of Gaussian initial data [km] — (Float64)
gaussian(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = A/exp((r - r0)^2/w^2)
gaussian_prime(r::Float64, A::Float64, r0::Float64, w::Float64)::Float64 = (2*A*(-r + r0))/(exp((r - r0)^2/w^2)*w^2)

xi_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for Lagrangian displacement (ξ) (used by perfect fluid / Eckart) — (Function)
du_ID(r::Float64)::Float64 = gaussian(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for perturbation of radial component of four-velocity (δu) (used by full BDNK) — (Function)
ddu_dr_ID(r::Float64)::Float64 = gaussian_prime(r, Gaussian_amplitude, Gaussian_center, Gaussian_width) # initial data function for radial derivative of δu (used by full BDNK) — (Function)

# time domain numerical parameters
KO = 0.1; # Kreiss-Oliger dissipation coefficient for BDNK evolution — (Float64)
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
PF_star = NeutronStarOscillations.Star(eps_central, kappa, n, 0.0, 0.0, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, du_ID, ddu_dr_ID, 0.0, CFL, dt_save_ms, h_save, save_every, total_time_ms);

# set relaxation times to zero for Eckart star as well as KO since only used in time integration of BDNK equations
Eckart_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, du_ID, ddu_dr_ID, 0.0, CFL, dt_save_ms, h_save, save_every, total_time_ms);

# BDNK star
BDNK_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, du_ID, ddu_dr_ID, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

#################### DEFAULT INITIAL DATA OPTIONS ####################
#=
    There are two default options for the initial data in ξ and δu. 
        
        
        (1) Gaussian initial data, like that specified explicitly above. One must set 'Gaussian_center', 'Gaussian_width', 'Gaussian_amplitude' to desired values.

            gaussian_ID_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, Gaussian_center, Gaussian_width, Gaussian_amplitude, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

        
        (2) Use eigenvectors from the perfect or Eckart fluid frequency domain solver. One must set 'mode' to the desired mode number (0 for fundamental, 1 for first overtone, etc.) and 'cowling' to true or false depending on whether one wants to use the Cowling approximation or not. Note that 'cowling=true' will only work for the perfect fluid case since there are no functions for an Eckart fluid in the Cowling approximation.

            mode = 0; cowling = false;
            eigvec_ID_star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, mode, cowling, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

        Warning about (2): this function computes eigenvectors using the frequency domain code and then interpolates them to create initial data functions. If the time domain grid extends slightly beyond the frequency domain grid, an error will be thrown since the interpolation functions are not defined outside the frequency domain grid. To avoid this, make sure the frequency domain code is ran at a higher resolution than the time domain code with the same (or smaller) pressure tolerance.
=#
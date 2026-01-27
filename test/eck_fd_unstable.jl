include("params.jl");
nPointsMatrix = 500;
h1_shoot = 2e-2;
h2_shoot = h1_shoot / 2;
h3_shoot = h1_shoot / 4;
h_TOV = 1e-4;
N_eigvals = 1;

cowling = false;
@time compute_eigensystem(Eckart_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = true);

η, ζ = 1.0, 10.0
function compute_eig(central_eps_factor)
    # set relaxation times to zero for Eckart star as well as KO since only used in time integration of BDNK equations
    Eckart_star = NeutronStarOscillations.Star(central_eps_factor * 1e15, kappa, n, η, ζ, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, xi_dt_ID, du_ID, du_dr_ID, du_dt_ID, de_ID, de_dr_ID, de_dt_ID, 0.0, CFL, dt_save_ms, h_save, save_every, total_time_ms);
    compute_eigensystem(Eckart_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = false);
    _, shooting_freqs_1, _, _, _, _ = load_eigensystem(Eckart_star, h1_shoot, cowling);
    return shooting_freqs_1[1]
end

using NonlinearSolve

NL_solver = NonlinearSolve.RobustMultiNewton()
NL_termination_condition = NonlinearSolve.AbsTerminationMode()

compute_eig(5.7)

# ε_init = 5.6;
# ε = [ε_init];
# resids = zero(ε)
# retcodes = String[];

# function bivariate_func(dx, x, params)
#     εc = x[1]
#     ω = compute_eig(εc)
#     dx[1] = real(ω)
#     dx[2] = imag(ω)
#     return nothing
# end

# NL_abstol = 1e-6;
# u0 = [ε_init, 0.0]
# using NonlinearSolve
# prob = NonlinearSolve.NonlinearProblem(bivariate_func, u0)
# sol = NonlinearSolve.solve(prob, NL_solver, termination_condition = NL_termination_condition,
# abstol = NL_abstol, maxiters = NL_maxiter)

# resids = zero(ε)
# retcodes = String[];

function f(u, p)
    return imag(compute_eig(u))
end

# function f(u, p)
#     return compute_eig(u)
# end

NL_abstol = 1e-14;
NL_maxiter = 100;
uspan = (5.6, 5.7);
using NonlinearSolve
prob = NonlinearSolve.IntervalNonlinearProblem(f, uspan)
@time sol = NonlinearSolve.solve(prob, NL_solver, termination_condition = NL_termination_condition,
abstol = NL_abstol, maxiters = NL_maxiter)

sol.u
sol.resid
sol.retcode

compute_eig(sol.u)


sol

ω[i] = sol.u[1] + im * sol.u[2]

# now compute eigenvector and boundary condition residual
rr, ξ, X, ξ_prime_BC = integrate(m, p, ε, ν, cs, cs_prime, r, ηTimesL, ζTimesL, ω[i]; shift = 1e-18)
evecs[:, i] .= ξ
resids[i] = ξ_prime_BC
push!(retcodes, string(sol.retcode))
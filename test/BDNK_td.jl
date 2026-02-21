include("params.jl")
res = 4e-2;
mode = 0;
KO = 0.2;
CFL = 0.1; # Courant-Friedrichs-Lewy factor — (Float64)
total_time_ms = 10.0; # total integration time [ms] — (Float64)
dt_save_ms = total_time_ms / 1000.0; # time interval between saved data points [ms] — (Float64)
save_every = 200; # save to file after 'save_every' time steps have been stored in memory (i.e., after every Δt = save_every * dt_save_ms) — (Int64)

cowling = false;
star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, τε, τP, τQ, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, mode, cowling, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);
@time time_integrate(star, res, cowling; print_progress = true);

sol = load_td_solution(star, res, cowling);
max_idx = 1000;
x = [sol["solution/t"][1:max_idx]];
y = [sol["solution/du"][1:max_idx, end]];
labels = ["KO=0.1", "KO=0.5"];
xlabel = L"t\,[\mathrm{ms}]";
ylabel = L"\delta{u}(R_{S})";
legend = true;
lim_x_min = nothing;
lim_x_max = nothing;
lim_y_min = nothing;
lim_y_max = nothing;
NeutronStarOscillations.QuickPlots.plot11(x, y;
    labels = labels,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    legend = legend,
    position = :lb)

# plot_characteristic_speeds(star, res, cowling)
# plot_initial_data_convergence(star, res, cowling)

time = 0.5;
var = "du";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\delta{u}";

plot_td_var(
    star, res, var, time, cowling;
    desample_factor=1,
    annotate_time=true,
    xlabel = xlabel,
    ylabel = ylabel
);

var = "du";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\delta{u}";
animate_td_var(
    star, res, var, cowling;
    desample_factor=1,
    stop_time=star.T, # stop time of animation in ms
    animation_length=10.0, # length of animation in seconds
    framerate=20, # frames per second
    xlabel = xlabel,
    ylabel = ylabel,
    fix_ylims = false, # whether to fix y-limits of plot
)

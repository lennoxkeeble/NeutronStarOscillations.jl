include("../../params.jl");
h1 = 4e-4;
h2 = h1 / 2.0;
h3 = h1 / 4.0;

# explicit RK4 solver (fourth-order accurate)
@time NeutronStarOscillations.TOV.Explicit.solve(BDNK_star, h1);
@time NeutronStarOscillations.TOV.Explicit.solve(BDNK_star, h2);
@time NeutronStarOscillations.TOV.Explicit.solve(BDNK_star, h3);

# implicit Crank Nicholson solver (second-order accurate)
@time NeutronStarOscillations.TOV.Implicit.solve(BDNK_star, h1);
@time NeutronStarOscillations.TOV.Implicit.solve(BDNK_star, h2);
@time NeutronStarOscillations.TOV.Implicit.solve(BDNK_star, h3);

# compute field convergence factors
hmin = h3;
@time r_explicit, Q_mass_explicit, Q_pressure_explicit, Q_eps_explicit, Q_nu_explicit = NeutronStarOscillations.TOV.compute_convergence_factor(BDNK_star, hmin, "Explicit");
@time r_implicit, Q_mass_implicit, Q_pressure_implicit, Q_eps_implicit, Q_nu_implicit = NeutronStarOscillations.TOV.compute_convergence_factor(BDNK_star, hmin, "Implicit");

# plot convergence factors for mass function from explicit and implicit solvers
x_plot = [r_explicit, r_implicit];
y_plot = [Q_mass_explicit, Q_mass_implicit];
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"Q_{N}";
lim_y_min = 0.0;
lim_y_max = 8.0;
legend = true;
labels = ["RK4", "CN"];
framevisible = true;
NeutronStarOscillations.QuickPlots.plot11(x_plot, y_plot;
    xlabel = xlabel,
    ylabel = ylabel,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    legend = legend,
    labels = labels,
    framevisible = framevisible
);

# plot convergence factors for all fields from explicit solver
x_plot = [r_explicit, r_explicit, r_explicit, r_explicit];
y_plot = [Q_mass_explicit, Q_pressure_explicit, Q_eps_explicit, Q_nu_explicit];
legend = true;
labels = [L"m", L"p", L"\epsilon", L"\nu"];
NeutronStarOscillations.QuickPlots.plot11(x_plot, y_plot;
    xlabel = xlabel,
    ylabel = ylabel,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    legend = legend,
    labels = labels,
    framevisible = framevisible
);
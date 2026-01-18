include("params.jl");
using Dierckx

h_TOV = 1e-4;
r_1, δε_1, δλ_1 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, 4.0h_TOV; u2_0 = 0.0);
r_2, δε_2, δλ_2 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, 2.0h_TOV; u2_0 = 0.0);
r_3, δε_3, δλ_3 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, h_TOV; u2_0 = 0.0);

# interpolate higher resolutions to coarsest grid
spline_order = 5; s = 0.0;
δε_spline_2 = Spline1D(r_2, δε_2; k=spline_order, s=s);
δλ_spline_2 = Spline1D(r_2, δλ_2; k=spline_order, s=s);
δε_spline_3 = Spline1D(r_3, δε_3; k=spline_order, s=s);
δλ_spline_3 = Spline1D(r_3, δλ_3; k=spline_order, s=s);

δλ_2_ds = [δλ_spline_2(r) for r in r_1];
δε_2_ds = [δε_spline_2(r) for r in r_1];
δλ_3_ds = [δλ_spline_3(r) for r in r_1];
δε_3_ds = [δε_spline_3(r) for r in r_1];

Q_δε = [abs(δε_1[i] - δε_2_ds[i]) / abs(δε_2_ds[i] - δε_3_ds[i]) for i in 1:length(δε_1)];
Q_δλ = [abs(δλ_1[i] - δλ_2_ds[i]) / abs(δλ_2_ds[i] - δλ_3_ds[i]) for i in 1:length(δλ_1)];

plot_ds = 100;
lim_x_min, lim_x_max, lim_y_min, lim_y_max = nothing, nothing, 0.0, 16.0;
x = [r_1[1:plot_ds:end], r_1[1:plot_ds:end]];
y = [Q_δε[1:plot_ds:end], Q_δλ[1:plot_ds:end]];
labels = [L"\delta\epsilon", L"\delta\lambda"];
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"Q_{N}";
legend = true;
NeutronStarOscillations.QuickPlots.plot11(x, y;
    labels = labels,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    legend = legend)

lim_x_min, lim_x_max, lim_y_min, lim_y_max = nothing, nothing, nothing, nothing;
x = [r_1[1:plot_ds:end], r_2[1:2*plot_ds:end], r_3[1:4*plot_ds:end]];
y = [δε_1[1:plot_ds:end], δε_2[1:2*plot_ds:end], δε_3[1:4*plot_ds:end]];
labels = ["Low", "Med", "High"];
ylabel = L"\delta\epsilon";
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

y = [δλ_1[1:plot_ds:end], δλ_2[1:2*plot_ds:end], δλ_3[1:4*plot_ds:end]];
ylabel = L"\delta\lambda";
NeutronStarOscillations.QuickPlots.plot11(x, y;
    labels = labels,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    legend = legend,
    position = :lt)
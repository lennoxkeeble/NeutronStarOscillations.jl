include("params.jl");
using Dierckx

h_TOV = 1e-2;
@time r_1, δλ_1 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, 4.0h_TOV; return_all=false);
@time r_2, δλ_2 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, 2.0h_TOV; return_all=false);
@time r_3, δλ_3 = NeutronStarOscillations.BDNKInitialData.compute_initial_data(BDNK_star, h_TOV; return_all=false);

# interpolate higher resolutions to coarsest grid
spline_order = 5; s = 0.0;
δλ_spline_2 = Spline1D(r_2, δλ_2; k=spline_order, s=s);
δλ_spline_3 = Spline1D(r_3, δλ_3; k=spline_order, s=s);

δλ_2_ds = [δλ_spline_2(r) for r in r_1];
δλ_3_ds = [δλ_spline_3(r) for r in r_1];

Q_δλ = [abs(δλ_1[i] - δλ_2_ds[i]) / abs(δλ_2_ds[i] - δλ_3_ds[i]) for i in 1:length(δλ_1)];

plot_ds = 1;
lim_x_min, lim_x_max, lim_y_min, lim_y_max = nothing, nothing, 0.0, 24.0;
x = [r_1[1:plot_ds:end]];
y = [Q_δλ[1:plot_ds:end]];
labels = [L"\delta\lambda"];
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
x = [r_1, r_2, r_3];
y = [δλ_1, δλ_2, δλ_3];
labels = ["Low", "Med", "High"];
xlabel = L"r\,[\mathrm{km}]";
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
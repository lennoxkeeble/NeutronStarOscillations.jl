include("../../params.jl");

# τQ_min(α::Float64)::Float64 = 1 / (1/α - 1)
# τP_func(τε::Float64, τQ::Float64, α::Float64)::Float64 = τε * (1/α - 1/τQ - 1)
# α = 0.8;
# τQ_min(α) |> println;
# τε_test = 10.0;
# τQ_test = 5.0;
# τP_func(τε_test, τQ_test, α) |> println

cs2_max_star = BDNK_star.dp_dε(BDNK_star.εc_SI)
cs2_max_causality = BDNK_star.τQ * BDNK_star.τε / (BDNK_star.τε + BDNK_star.τQ * (BDNK_star.τP + BDNK_star.τε))
cs2_max_star < cs2_max_causality ? println("Causal choice of frame") : @warn("Acausal choice of frame")

h1 = h_save;
h2 = h_save / 2.0;
h3 = h_save / 4.0;

@time NeutronStarOscillations.TimeDomain.BDNK.solve(BDNK_star, h1);
@time NeutronStarOscillations.TimeDomain.BDNK.solve(BDNK_star, h2);
@time NeutronStarOscillations.TimeDomain.BDNK.solve(BDNK_star, h3);

# load solution and print keys
sol1 = NeutronStarOscillations.TimeDomain.load_solution(BDNK_star, h1);
keys(sol1["solution"]) |> println;
close(sol1)

# plot characteristic speeds throughout star
NeutronStarOscillations.TimeDomain.plot_characteristic_speeds(BDNK_star, h3);

# plot convergence factor for initial data
NeutronStarOscillations.TimeDomain.plot_initial_data_convergence(BDNK_star, h3);

# plot solution at fixed time
time = 0.0;
var = "deps";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\delta{\epsilon}";

NeutronStarOscillations.TimeDomain.plot_var(
    BDNK_star, h3, var, time;
    desample_factor=1,
    annotate_time=true,
    xlabel = xlabel,
    ylabel = ylabel
);

# plot one-norm of independent residual for three different resolutions
resolutions = [h1, h2, h3];
var = "Constrained_IR1";
xlabel = L"t\,[\mathrm{ms}]";
ylabel = L"||\mathrm{IR}||_1";
labels = ["Lev1", "Lev2", "Lev3"];
legend = true;
position = :rt;
NeutronStarOscillations.TimeDomain.plot_IR_one_norms(
    BDNK_star, resolutions, var, time;
    desample_factor=1,
    xlabel = xlabel,
    ylabel = ylabel,
    yscale = log10,
    labels = labels,
    legend = legend,
    position = position
);

# make animation — saved to star.fig_path
var = "du";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\delta{u}";
NeutronStarOscillations.TimeDomain.animate_var(
    BDNK_star, h3, var;
    desample_factor=1,
    stop_time=BDNK_star.T, # stop time of animation in ms
    animation_length=10.0, # length of animation in seconds
    framerate=10, # frames per second
    xlabel = xlabel,
    ylabel = ylabel,
    fix_ylims = false, # whether to fix y-limits of plot
)
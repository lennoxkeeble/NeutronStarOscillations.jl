include("../../params.jl");
h1 = h_save;
h2 = h_save / 2.0;
h3 = h_save / 4.0;

@time NeutronStarOscillations.TimeDomain.Eckart.solve(Eckart_star, h1);
@time NeutronStarOscillations.TimeDomain.Eckart.solve(Eckart_star, h2);
@time NeutronStarOscillations.TimeDomain.Eckart.solve(Eckart_star, h3);

# load solution and print keys
sol1 = NeutronStarOscillations.TimeDomain.load_solution(Eckart_star, h1);
keys(sol1["solution"]) |> println;
close(sol1)

# plot solution at fixed time
time = 0.05;
var = "xi";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";

NeutronStarOscillations.TimeDomain.plot_var(
    Eckart_star, h3, var, time;
    desample_factor=1,
    annotate_time=true,
    xlabel = xlabel,
    ylabel = ylabel
);

# plot one-norm of independent residual for three different resolutions
resolutions = [h1, h2, h3];
var = "IR";
xlabel = L"t\,[\mathrm{ms}]";
ylabel = L"||\mathrm{IR}||_1";
labels = ["Lev1", "Lev2", "Lev3"];
legend = true;
position = :lt;
NeutronStarOscillations.TimeDomain.plot_IR_one_norms(
    Eckart_star, resolutions, var, time;
    desample_factor=1,
    xlabel = xlabel,
    ylabel = ylabel,
    yscale = log10,
    labels = labels,
    legend = legend,
    position = position
);

# make animation — saved to star.fig_path
var = "xi";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";
NeutronStarOscillations.TimeDomain.animate_var(
    Eckart_star, h3, var;
    desample_factor=1,
    stop_time=Eckart_star.T, # stop time of animation in ms
    animation_length=10.0, # length of animation in seconds
    framerate=10, # frames per second
    xlabel = xlabel,
    ylabel = ylabel,
    fix_ylims = false, # whether to fix y-limits of plot
)
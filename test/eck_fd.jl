include("params.jl")
nPointsMatrix = 500;
h1_shoot = 2e-2;
h2_shoot = h1_shoot / 2;
h3_shoot = h1_shoot / 4;
h_TOV = 1e-5;

cowling = false;
@time compute_eigensystem(Eckart_star, h1_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = true);
@time compute_eigensystem(Eckart_star, h2_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = true);
@time compute_eigensystem(Eckart_star, h3_shoot, h_TOV, N_eigvals, nPointsMatrix, cowling; print_progress = true);

# load eigensystems. Eigenvalues are (complex-valued) linear frequencies, i.e., ν = ω / (2π) and [ν] = kHz
init_freqs_1, shooting_freqs_1, radii_1, evecs_1, resids_1, ret_1 = load_eigensystem(Eckart_star, h1_shoot, cowling);
init_freqs_2, shooting_freqs_2, radii_2, evecs_2, resids_2, ret_2 = load_eigensystem(Eckart_star, h2_shoot, cowling);
init_freqs_3, shooting_freqs_3, radii_3, evecs_3, resids_3, ret_3 = load_eigensystem(Eckart_star, h3_shoot, cowling);

# print return codes for the shooting method iterations
println("*** NonlinearSolve return codes for shooting method ***")
println("(Low res) Eckart shooting method return codes: ", ret_1)
println("(Med res) Eckart shooting method return codes: ", ret_2)
println("(High res) Eckart shooting method return codes: ", ret_3)

# print highest resolution frequencies
println("*** Linear frequencies ν = ω / 2π, [v] = κHz ***")
println("(Low res) Eckart fluid frequencies (kHz): ", shooting_freqs_1)
println("(Med res) Eckart fluid frequencies (kHz): ", shooting_freqs_2)
println("(High res) Eckart fluid frequencies (kHz): ", shooting_freqs_3)

# plot eigenvectors
modes = [0, 1, 2, 3];
method = "Shoot";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";
legend = true;
labels = [L"n=%$mode" for mode in modes];
framevisible = true;

plot_fd_eigvecs(
    modes, Eckart_star, method; h=h3_shoot,
    xlabel = xlabel,
    ylabel = ylabel,
    legend = legend,
    labels = labels,
    position = :lb,
    framevisible = framevisible
);
include("../../params.jl");

#=
    compute eigenvalues and eigenvectors via matrix + shooting method at three different resolutions. The shooting function runs a matrix method first to locate the
    eigenvalues and uses these as the initial guess in the shooting method which consists of a nonlinear rootfinding in ω such that Δp = 0 at the surface.
=#
  
N_matrix = 500;
h1_shoot = 4e-2;
h2_shoot = h1_shoot / 2;
h3_shoot = h1_shoot / 4;
h_TOV = 1e-5;

@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.get_eigensystem(PF_star, h1_shoot, h_TOV, N_eigvals; nPointsMatrix = N_matrix, print_progress = true);
@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.get_eigensystem(PF_star, h2_shoot, h_TOV, N_eigvals; nPointsMatrix = N_matrix, print_progress = true);
@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.get_eigensystem(PF_star, h3_shoot, h_TOV, N_eigvals; nPointsMatrix = N_matrix, print_progress = true);

# load eigensystems. Eigenvalues are linear frequencies squared, i.e., ν^2 where ν = ω / (2π) and [ν] = kHz
init_freqs_1, shooting_freqs_1, radii_1, evecs_1, resids_1, ret_1 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.load_eigensystem(PF_star, h1_shoot);
init_freqs_2, shooting_freqs_2, radii_2, evecs_2, resids_2, ret_2 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.load_eigensystem(PF_star, h2_shoot);
init_freqs_3, shooting_freqs_3, radii_3, evecs_3, resids_3, ret_3 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.load_eigensystem(PF_star, h3_shoot);

# check convergence
diff_12 = abs.(shooting_freqs_1 .- shooting_freqs_2)
diff_23 = abs.(shooting_freqs_2 .- shooting_freqs_3)
prod(@. ((diff_23 < diff_12) || diff_23 == diff_12 == 0.0)) ? println("FREQUENCIES CONVERGING") : println("FREQUENCIES NOT CONVERGING")

# linear frequencies in kHz (if the equilibrium configuration is unstable, then the squared frequencies will be negative and the code below will throw an error due to the sqrt)
freqs_1 = sqrt.(shooting_freqs_1);
freqs_2 = sqrt.(shooting_freqs_2);
freqs_3 = sqrt.(shooting_freqs_3);

# print return codes for the shooting method iterations
println("*** NonlinearSolve return codes for shooting method ***")
println("(Low res) PF shooting method return codes: ", ret_1)
println("(Med res) PF shooting method return codes: ", ret_2)
println("(High res) PF shooting method return codes: ", ret_3)

# print highest resolution frequencies
println("*** Linear frequencies ν = ω / 2π, [v] = κHz ***")
println("(Low res) PF fluid frequencies (kHz): ", freqs_1)
println("(Med res) PF fluid frequencies (kHz): ", freqs_2)
println("(High res) PF fluid frequencies (kHz): ", freqs_3)

# plot eigenvectors
modes = [0, 1, 2, 3];
method = "Shoot Cowling";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";
legend = true;
labels = [L"n=%$mode" for mode in modes];
framevisible = true;

NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
    modes, PF_star, method; h=h3_shoot,
    xlabel = xlabel,
    ylabel = ylabel,
    legend = legend,
    labels = labels,
    position = :lb,
    framevisible = framevisible
);
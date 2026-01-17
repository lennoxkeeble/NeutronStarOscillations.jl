include("../../params.jl");
h_TOV = 1e-5;
N_eigvals = 5;

N1 = 500;
N2 = 1000;
N3 = 2000;

# compute eigenvalues and eigenvectors at three different resolutions
@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.get_eigensystem(PF_star, N1, h_TOV, N_eigvals);
@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.get_eigensystem(PF_star, N2, h_TOV, N_eigvals);
@time NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.get_eigensystem(PF_star, N3, h_TOV, N_eigvals);

# load eigensystems. Eigenvalues are linear frequencies squared, i.e., ν^2 where ν = ω / (2π) and [ν] = kHz
omega_squared_1, evecs_1, r_1 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.load_eigensystem(PF_star, N1);
omega_squared_2, evecs_2, r_2 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.load_eigensystem(PF_star, N2);
omega_squared_3, evecs_3, r_3 = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.load_eigensystem(PF_star, N3);
isapprox(omega_squared_3, NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.get_eigenvalues(PF_star, N3, h_TOV, N_eigvals))

# check convergence
diff_12 = abs.(omega_squared_1 .- omega_squared_2);
diff_23 = abs.(omega_squared_2 .- omega_squared_3);
prod(@. ((diff_23 < diff_12) || diff_23 == diff_12 == 0.0)) ? println("FREQUENCIES CONVERGING") : println("FREQUENCIES NOT CONVERGING")

# linear frequencies in kHz (if the equilibrium configuration is unstable, then the squared frequencies will be negative and the code below will throw an error due to the sqrt)
freqs_1 = sqrt.(omega_squared_1);
freqs_2 = sqrt.(omega_squared_2);
freqs_3 = sqrt.(omega_squared_3);

# print frequencies
println("(Low res) frequencies (kHz): ", freqs_1)
println("(Med res) frequencies (kHz): ", freqs_2)
println("(High res) frequencies (kHz): ", freqs_3)

# plot eigenvectors
modes = [0, 1, 2, 3];
method = "Matrix Cowling";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";
legend = true;
labels = [L"n=%$mode" for mode in modes];
framevisible = true;

NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
    modes, PF_star, method; N=N3,
    xlabel = xlabel,
    ylabel = ylabel,
    legend = legend,
    labels = labels,
    position = :lb,
    framevisible = framevisible
);
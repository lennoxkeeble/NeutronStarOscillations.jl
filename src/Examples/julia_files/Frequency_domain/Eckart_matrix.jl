include("../../params.jl");
h_TOV = 1e-5;
N_eigvals = 5;

N1 = 500;
N2 = 750;
N3 = 1000;

# compute eigenvalues and eigenvectors at three different resolutions
@time NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.get_eigensystem(Eckart_star, N1, h_TOV, N_eigvals);
@time NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.get_eigensystem(Eckart_star, N2, h_TOV, N_eigvals);
@time NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.get_eigensystem(Eckart_star, N3, h_TOV, N_eigvals);

# load eigensystems. Eigenvalues are (complex-valued) linear frequencies, i.e., ν = ω / (2π) and [ν] = kHz
omega_1, evecs_1, r_1 = NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.load_eigensystem(Eckart_star, N1);
omega_2, evecs_2, r_2 = NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.load_eigensystem(Eckart_star, N2);
omega_3, evecs_3, r_3 = NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.load_eigensystem(Eckart_star, N3);
isapprox(omega_3, NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.get_eigenvalues(Eckart_star, N3, h_TOV, N_eigvals))

# check convergence
diff_12 = abs.(omega_1 .- omega_2);
diff_23 = abs.(omega_2 .- omega_3);
prod(@. ((diff_23 < diff_12) || diff_23 == diff_12 == 0.0)) ? println("FREQUENCIES CONVERGING") : println("FREQUENCIES NOT CONVERGING")

# print frequencies
println("(Low res) frequencies (kHz): ", omega_1)
println("(Med res) frequencies (kHz): ", omega_2)
println("(High res) frequencies (kHz): ", omega_3)

# plot eigenvectors
modes = [0, 1, 2, 3];
method = "Matrix";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";
legend = true;
labels = [L"n=%$mode" for mode in modes];
framevisible = true;

NeutronStarOscillations.FrequencyDomain.plot_eigvecs(
    modes, Eckart_star, method; N=N3,
    xlabel = xlabel,
    ylabel = ylabel,
    legend = legend,
    labels = labels,
    position = :lb,
    framevisible = framevisible
);
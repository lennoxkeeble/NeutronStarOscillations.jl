#=

    Module comprising of functions which solve the perfect fluid and Eckart equations of motion governing linear radial perturbations of polytropic neutron stars in the
    frequency domain. We provide functions which implement matrix and shooting methods. The functions which we intend for user use combine these methods.
    The matrix method is used to locate the eigenvalues (i.e., to provide a good initial guess). The matrix method output is then input into the shooting method
    which iterates over the value of the eigenvalue, driving the residual (BC that the lagrangian displacement of the pressure vanishes at the surface) to zero.

=#

module FrequencyDomain
using NeutronStarOscillations
using Printf
using NonlinearSolve

NL_solver = NonlinearSolve.RobustMultiNewton()
# NL_termination_condition = NonlinearSolve.RelTerminationMode()
NL_termination_condition = NonlinearSolveBase.AbsTerminationMode()

# Dierckx interpolation parameters
const spline_order::Int64 = 5;
const s::Float64 = 0.0;
const bc::String = "error"; # boundary condition for spline interpolation ("nearest", "zero", "extrapolate", "error")

shoot_fname(star::NeutronStarOscillations.Star, h::Float64)::String = star.data_path * @sprintf("Frequency_domain_matrix_κ_%d_n_%.1f_εc_%.2e_η_%.1e_ζ_%.1e_L_%s_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.η, star.ζ, star.L, h);
matrix_fname(star::NeutronStarOscillations.Star, N::Int64)::String = star.data_path * @sprintf("Frequency_domain_matrix_κ_%d_n_%.1f_εc_%.2e_η_%.1e_ζ_%.1e_L_%s_N_%s.h5", star.kappa, star.n, star.εc_SI, star.η, star.ζ, star.L, N);
cowling_shoot_fname(star::NeutronStarOscillations.Star, h::Float64)::String = star.data_path * @sprintf("Frequency_domain_cowling_matrix_κ_%d_n_%.1f_εc_%.2e_η_%.1e_ζ_%.1e_L_%s_h_%.1e.h5", star.kappa, star.n, star.εc_SI, star.η, star.ζ, star.L, h);
cowling_matrix_fname(star::NeutronStarOscillations.Star, N::Int64)::String = star.data_path * @sprintf("Frequency_domain_cowling_matrix_κ_%d_n_%.1f_εc_%.2e_η_%.1e_ζ_%.1e_L_%s_N_%s.h5", star.kappa, star.n, star.εc_SI, star.η, star.ζ, star.L, N);


plot_fname(star::NeutronStarOscillations.Star, N::Int64)::String = star.fig_path * @sprintf("Frequency_domain_eigvecs_κ_%d_n_%.1f_εc_%.2e_η_%.1e_ζ_%.1e_L_%s_N_%s.png", star.kappa, star.n, star.εc_SI, star.η, star.ζ, star.L, N);
function plot_eigvecs(modes::Vector{Int64}, star::NeutronStarOscillations.Star, method::String; N::Int64=0, h::Float64=0.0,
    desample_factor::Int64=1,
    colors = NeutronStarOscillations.QuickPlots.default_colors,
    labels = NeutronStarOscillations.QuickPlots.default_labels,
    linestyles = NeutronStarOscillations.QuickPlots.default_linestyles,
    linewidths = NeutronStarOscillations.QuickPlots.default_linewidths,
    alphas = NeutronStarOscillations.QuickPlots.default_alphas,
    xlabel = "",
    ylabel = "",
    lim_x_min = nothing,
    lim_x_max = nothing,
    lim_y_min = nothing,
    lim_y_max = nothing,
    yscale = identity,
    text_xloc = nothing, # loacation of text overlayed on plot
    text_yloc = nothing, # loacation of text overlayed on plot
    text = nothing, # text overlayed on plot
    xalign = :left, # alignment of text overlayed on plot
    yalign = :top, # alignment of text overlayed on plot
    legend = false,
    position = :rt, # position of legend
    labelsize = NeutronStarOscillations.QuickPlots.fontsize, # size of legend and axis labels
    framevisible = true, # whether to draw box around legend
    hlines = nothing, # y-values to draw horizontal lines at in the form [y1, y2, ...]
    vlines = nothing, # x-values to draw vertical lines at in the form [x1, x2, ...]
    save_plot = false, # whether to save the plot to file
    scatter_lines = false # whether to plot scatter lines instead of normal lines
    )

    
    if isequal(method, "Matrix")
        if star.η == 0.0 && star.ζ == 0.0
            eigvals, evecs, r = NeutronStarOscillations.FrequencyDomain.PerfectFluid.Matrix.load_eigensystem(star, N);
        else
            eigvals, evecs, r = NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.load_eigensystem(star, N);
        end
    elseif isequal(method, "Shoot")
        if star.η == 0.0 && star.ζ == 0.0
            init_freqs, eigvals, r, evecs, resids, ret = NeutronStarOscillations.FrequencyDomain.PerfectFluid.Shoot.load_eigensystem(star, h);
        else
            init_freqs, eigvals, r, evecs, resids, ret = NeutronStarOscillations.FrequencyDomain.Eckart.Shoot.load_eigensystem(star, h);
        end
    elseif isequal(method, "Matrix Cowling")
        if star.η == 0.0 && star.ζ == 0.0
            eigvals, evecs, r = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.load_eigensystem(star, N);
        else
            error("Frequency domain Cowling only implemented for perfect fluid stars");
        end
    elseif isequal(method, "Shoot Cowling")
        if star.η == 0.0 && star.ζ == 0.0
            init_freqs, eigvals, r, evecs, resids, ret = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Shoot.load_eigensystem(star, h);
        else
            error("Frequency domain Cowling only implemented for perfect fluid stars");
        end
    end

    x = Vector{Float64}[];
    y = Vector{Float64}[];

    for i in eachindex(modes)
        mode = modes[i];
        push!(x, r[1:desample_factor:end]);
        push!(y, real.(evecs[:, mode + 1][1:desample_factor:end]));
    end

    fname = plot_fname(star, N);

    NeutronStarOscillations.QuickPlots.plot11(x, y;
    colors = colors,
    labels = labels,
    linestyles = linestyles,
    linewidths = linewidths,
    alphas = alphas,
    xlabel = xlabel,
    ylabel = ylabel,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
    yscale = yscale,
    text_xloc = text_xloc,
    text_yloc = text_yloc,
    text = text,
    xalign = xalign,
    yalign = yalign,
    legend = legend,
    position = position,
    labelsize = labelsize,
    framevisible = framevisible,
    hlines = hlines,
    vlines = vlines,
    save_plot = save_plot,
    fname = fname, 
    scatter_lines = scatter_lines
    )
end

# normalize eigenvectors such that |real(ξ(R))| = 1 km
function normalize_eigvecs!(evecs::AbstractArray)
    for i in 1:size(evecs, 2)
        evecs[:, i] .= evecs[:, i] ./ abs(real(evecs[end, i]));
    end
end

module PerfectFluidCowling

module Matrix
using LinearAlgebra
using Dierckx
using HDF5
using NeutronStarOscillations
using ...FrequencyDomain

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;

# write ODE as A2(r) ξ''(r) + A1(r) ξ'(r) + A0(r) ξ(r) = ω^2 ξ(r)
A0(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = (2*exp(ν)*(m^2*(1 + 5*cs*(cs - 2*r*cs_prime)) + r*m*(-1 - 8*π*r^2*p + cs*((-5 + 8*π*r^2*ε)*cs + r*(9 + 8*π*r^2*ε)*cs_prime)) + r^2*(cs^2 + 2*r*(π*r*(p + ε + 8*π*r^2*p*ε - 8*π*r^2*ε^2*cs^2) - (1 + 2*π*r^2*ε)*cs*cs_prime))))/(r^3*(r - 2*m))
A1(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = (exp(ν)*(4*π*r^3*p - 2*r*cs*((1 + 2*π*r^2*ε)*cs + r*cs_prime) + m*(1 + 5*cs^2 + 4*r*cs*cs_prime)))/r^2
A2(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = -((exp(ν)*(r - 2*m)*cs^2)/r)

# discretized equation at interior point i: Ai_minus1*ξ_{i-1} + Aii*ξ_{i} + Ai_plus1*ξ_{i+1} = ω^2 * ξ_{i}
Ai_minus1(A1::Float64, A2::Float64, A3::Float64, h::Float64)::Float64 = (2*A3 - A2*h)/(2.0*h^2)
Aii(A1::Float64, A2::Float64, A3::Float64, h::Float64)::Float64 = A1 - (2*A3)/h^2
Ai_plus1(A1::Float64, A2::Float64, A3::Float64, h::Float64)::Float64 = (2*A3 + A2*h)/(2.0*h^2)

# in the final row, enforce Δp = 0 at surface. Get an expression like aξ + bξ' = 0. Discretize this with centered differences, solve for ξ[R+h] and eliminate from equation discretized
# at final grid point, which is also discretized with centered differences.
An_minus1(A1::Float64, A2::Float64, A3::Float64, h::Float64)::Float64 = (2*A3)/h^2
Ann(A1::Float64, A2::Float64, A3::Float64, m::Float64, p::Float64, ε::Float64, r::Float64, h::Float64)::Float64 = ((2*A3*(5*h + 2*r) + h^2*(5*A2 - 2*A1*r))*m - r*(2*A3*(2*h + r) + h^2*(2*A2 - A1*r) + 4*h*(2*A3 + A2*h)*π*r^2*ε))/(h^2*r*(r - 2*m))

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function get_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)
    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluidCowling.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    matrix_fname = FrequencyDomain.cowling_matrix_fname(star, nPoints); 
    FrequencyDomain.PerfectFluidCowling.Matrix.get_eigensystem(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, N_eigvals; fname=matrix_fname, save_to_file=true);
end

function get_eigensystem(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, N_eigvals::Int64; fname::String = "", save_to_file::Bool=false)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    matrix_size = length(m_arr) - 1;
    A = zeros(matrix_size, matrix_size); # create a square matrix of
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h);
    eigsyst = eigen(A);
    eig_vals = eigsyst.values;
    eig_vecs = eigsyst.vectors;

    perm = sortperm(eig_vals);
    omega_squared = eig_vals[perm][1:N_eigvals];
    vecs = eig_vecs[:, perm][:, 1:N_eigvals];
    # don't solve for ξ at center so prepend this to arrays
    evecs = vcat(zeros(size(vecs, 2))', vecs)

    # normalize eigenvectors
    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Filename must be provided to save eigensystem."))
        end
        save_eigensystem(fname, omega_squared / (2π * kHz_to_km)^2, evecs, r)
    else
        return omega_squared / (2π * kHz_to_km)^2, evecs, r
    end
end

function save_eigensystem(fname::String, eigvals::Vector{Float64}, eigvecs::AbstractArray{Float64}, radius::Vector{Float64})
    h5open(fname, "w") do file
        file["eigvals"] = eigvals;
        file["eigvecs"] = eigvecs;
        file["radius"] = radius;
    end
end

function load_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64)
    fname = FrequencyDomain.cowling_matrix_fname(star, nPoints)
    h5f = h5open(fname, "r")
        eigvals = h5f["eigvals"][:];
        eigvecs = h5f["eigvecs"][:,:];
        radius = h5f["radius"][:];
    close(h5f)
    return eigvals, eigvecs, radius
end

function get_eigenvalues(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)
    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluidCowling.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    return FrequencyDomain.PerfectFluidCowling.Matrix.get_eigenvalues(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, N_eigvals);
end

# intended for computing eigenvalues as input to shooting method
function get_eigenvalues(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, N_eigvals::Int64)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    matrix_size = length(m_arr) - 1;
    A = zeros(matrix_size, matrix_size); # create a square matrix of
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h);
    eig_vals = eigvals(A);

    perm = sortperm(eig_vals);
    omega_squared = eig_vals[perm][1:N_eigvals];

    return omega_squared / (2π * kHz_to_km)^2
end

# fill matrix used in matrix method
function fill_matrix!(A::AbstractArray{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, h::Float64)
    # A is an n x n matrix. Row i of matrix A corresponds to the pulsaiton equation discretized at radius r[i+1]. Note that this index shifting is because the equation is solved at r=0 since ξ(0) = 0 is known.
    n = length(m)-1 

    # first row (discretized eqn at r[2], using ξ(0) = 0)
    A0 = Matrix.A0(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2]);
    A1 = Matrix.A1(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2]);
    A2 = Matrix.A2(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2]);

    A[1, 1] = Aii(A0, A1, A2, h)
    A[1, 2] = Ai_plus1(A0, A1, A2, h)

    # last row
    A0 = Matrix.A0(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end]);
    A1 = Matrix.A1(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end]);
    A2 = Matrix.A2(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end]);
    
    # # ghost BCs
    # A[end, end-1] = Ai_minus1(A0, A1, A2, h)
    # A[end, end] = Aii(A0, A1, A2, h)

    # enforcing lagrangian pressure variation Δp = 0 at surface
    A[end, end-1] = An_minus1(A0, A1, A2, h)
    A[end, end] = Ann(A0, A1, A2, m[end], p[end], ε[end], r[end], h)

    # interior points
    for i in 2:(n - 1)
        mi = m[i+1]; pressure_i = p[i+1]; εi = ε[i+1]; νi = ν[i+1]; csi = cs[i+1]; cs_prime_i = cs_prime[i+1]; ri = r[i+1];
        A0 = Matrix.A0(mi, pressure_i, εi, νi, csi, cs_prime_i, ri);
        A1 = Matrix.A1(mi, pressure_i, εi, νi, csi, cs_prime_i, ri);
        A2 = Matrix.A2(mi, pressure_i, εi, νi, csi, cs_prime_i, ri);
        A[i, i-1] = Ai_minus1(A0, A1, A2, h);
        A[i, i] = Aii(A0, A1, A2, h);
        A[i, i+1] = Ai_plus1(A0, A1, A2, h);
    end
end

end

module Shoot
using LinearAlgebra
using Dierckx
using HDF5
using NonlinearSolve
using StaticArrays
using NeutronStarOscillations
using ...FrequencyDomain
using ...PerfectFluidCowling
using Printf

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;

function get_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64, h_TOV::Float64, N_eigvals::Int64; nPointsMatrix::Int64, print_progress::Bool=false)
    # extract frequencies from matrix method
    matrix_freqs = NeutronStarOscillations.FrequencyDomain.PerfectFluidCowling.Matrix.get_eigenvalues(star, nPointsMatrix, h_TOV, N_eigvals)
    ω_init = matrix_freqs * (2π * kHz_to_km)^2

    # compute high resolution TOV solution and interpolate to desired computational grid
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r = range(start = 0.0, stop = Rs, step = h_shoot) |> collect;

    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    m = @. m_spline(r);
    p = @. p_spline(r);
    ε = @. ε_spline(r);
    ν = @. ν_spline(r);
    cs = @. sqrt(star.dp_dε(ε_spline(r)));
    cs_prime = @. FrequencyDomain.PerfectFluidCowling.Matrix.cs_prime_func(m, p, ε, cs, r, star.d2p_dε2);
    cs_prime[1] = 0.0;

    # compute eigenvalues and eigenvectors via matrix + shooting method
    fname = FrequencyDomain.cowling_shoot_fname(star, h_shoot);
    PerfectFluidCowling.Shoot.iterate_frequencies(m, p, ε, ν, cs, cs_prime, r, ω_init, star.NL_reltol, star.NL_abstol, star.NL_maxiter; fname=fname, save_to_file=true, print_progress = print_progress);
end

function iterate_frequencies(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω_init::Vector{<:Number}, NL_reltol::Float64, NL_abstol::Float64, NL_maxiter::Int64; fname::String = "", save_to_file::Bool=false, print_progress::Bool=false)
    # now use root finding to solve for ω with initial condition given by (low resolution) matrix method
    N_eigvals = length(ω_init);
    ω = zero(ω_init)

    TOV_length = length(r); # minus four comes from removing the additional four steps in the background solution that allows us to evolve the pulsation equations until the surface of the star (since RK4 requires background solution at r + h/2 and r + h)
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    evecs = zeros(ComplexF64, pulsation_length, N_eigvals)
    resids = zero(ω_init)
    retcodes = String[];

    rr = [];

    function bivariate_func(dx, x, params)
        ω0 = x[1]
        ξ_prime_BC = fast_integrate(m, p, ε, ν, cs, cs_prime, r, ω0; shift=1e-18);
        dx[1] = ξ_prime_BC
        return nothing
    end

    for i in 1:N_eigvals
        print_progress ? println("Root finding eigenvalue $i / $N_eigvals") : nothing
        u0 = [real(ω_init[i])]
        prob = NonlinearProblem(bivariate_func, u0)
        sol = solve(prob, NeutronStarOscillations.FrequencyDomain.NL_solver, termination_condition = NeutronStarOscillations.FrequencyDomain.NL_termination_condition,
        reltol = NL_reltol, maxiters = NL_maxiter)
        ω[i] = sol.u[1]

        # now compute eigenvector and boundary condition residual
        rr, ξ, X, ξ_prime_BC = integrate(m, p, ε, ν, cs, cs_prime, r, ω[i]; shift = 1e-18)
        evecs[:, i] .= ξ
        resids[i] = ξ_prime_BC
        push!(retcodes, string(sol.retcode))
    end
    print_progress ? println("Computed all eigenvalues") : nothing

    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Must provide filename to save eigenvalues"))
        else
            save_eigensystem(fname, ω_init ./ (2π * kHz_to_km)^2, ω ./ (2π * kHz_to_km)^2, rr, evecs, resids, retcodes)
        end
    else
        return ω_init ./ (2π * kHz_to_km)^2, ω ./ (2π * kHz_to_km)^2, rr, evecs, resids, retcodes
    end
end

function save_eigensystem(fname::String, ω_init::Vector{<:Number}, shoot_eigvals::Vector{<:Number}, radius::Vector{Float64}, evecs::LinearAlgebra.Matrix{<:Number}, resids::Vector{<:Number}, retcodes::Vector{String})
    h5open(fname, "w") do file
        file["ω_init"] = ω_init;
        file["shoot_eigvals"] = shoot_eigvals;
        file["radius"] = radius;
        file["eigvecs"] = evecs;
        file["resids"] = resids;
        file["retcodes"] = retcodes;
    end;
end

function load_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64)
    fname = FrequencyDomain.cowling_shoot_fname(star, h_shoot);
    h5f = h5open(fname, "r")
        ω_init = h5f["ω_init"][:];
        shoot_eigvals = h5f["shoot_eigvals"][:];
        radius = h5f["radius"][:];
        evecs = h5f["eigvecs"][:, :];
        resids = h5f["resids"][:];
        retcodes = h5f["retcodes"][:];
    close(h5f)
    return ω_init, shoot_eigvals, radius, evecs, resids, retcodes
end

# RK4 equations. Integrate second-order ODE in space. To reduce to first-order, define X ≡ ξ'[r], ω2 = ω^2
function X_prime(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, ξ, X, ω2)
    A1 = PerfectFluidCowling.Matrix.A0(m, p, ε, ν, cs, cs_prime, r);
    A2 = PerfectFluidCowling.Matrix.A1(m, p, ε, ν, cs, cs_prime, r);
    A3 = PerfectFluidCowling.Matrix.A2(m, p, ε, ν, cs, cs_prime, r);
    return (-(A2*X) + (-A1 + ω2)*ξ)/A3
end
ξ_prime(X) = X

# homogeneous boundary condition to enforce Lagrangian variation of pressure Δp = 0 at surface. Nonlinear rootfinder finds ω^2 such that this quantity is zero.
ξ_prime_BC(m::Float64, ε::Float64, r::Float64, ξ) = ((5*m - 2*(r + 2*π*r^3*ε))*ξ)/(r*(r - 2*m))

# TOV equations
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function fast_integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω2; shift::Float64=1e-18)
    h = diff(r)[1]; # step size for the TOV solution
    H = 2.0h;

    # initial conditions
    ξi = 0.0
    Xi = 1e-10

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(Number, 2);
    k2 = zeros(Number, 2);
    k3 = zeros(Number, 2);
    k4 = zeros(Number, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξi;
        X0 = Xi;
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, ξ0, X0, ω2) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0
    
        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, ξ1, X1, ω2);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, ξ2, X2, ω2);


        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, ξ3, X3, ω2);

        # update values
        ξi = ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        Xi = X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    ε0 = ε[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξi;
    X0 = Xi;

    BC_residual = X0 - ξ_prime_BC(m0, ε0, r0, ξ0)

    return BC_residual
end

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω2; shift::Float64=1e-18)
    h = diff(r)[1];
    H = 2.0h;

    # initial conditions
    ξ = [0.0];
    X = [1e-10];
    rr = [r[1]];

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(ComplexF64, 2);
    k2 = zeros(ComplexF64, 2);
    k3 = zeros(ComplexF64, 2);
    k4 = zeros(ComplexF64, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξ[end];
        X0 = X[end];
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, ξ0, X0, ω2) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0

        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, ξ1, X1, ω2);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, ξ2, X2, ω2);

        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, ξ3, X3, ω2);

        # update values
        push!(rr, r0 + H);
        push!(ξ, ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]));
        push!(X, X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]));

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    ε0 = ε[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξ[end];
    X0 = X[end];

    BC_residual = X0 - ξ_prime_BC(m0, ε0, r0, ξ0)

    return rr, ξ, X, BC_residual
end

end


end

module PerfectFluid
module Matrix
using LinearAlgebra
using Dierckx
using HDF5
using NeutronStarOscillations
using ...FrequencyDomain

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;

# discretized equation at interior point i: Ai_minus1*ξ_{i-1} + Aii*ξ_{i} + Ai_plus1*ξ_{i+1} = ω^2 * ξ_{i}
Ai_minus1(m::Float64, p::Float64, ε0::Float64, ν0::Float64, cs::Float64, cs_prime::Float64, r::Float64, h::Float64)::Float64 = (exp(ν0)*(h*(m + 4*π*r^3*p) + cs^2*((5*h - 4*r)*m - 2*r*(h - r + 2*h*π*r^2*ε0)) - 2*h*r*cs*(r - 2*m)*cs_prime))/(2.0*h^2*r^2)
Aii(m::Float64, p::Float64, ε0::Float64, ν0::Float64, cs::Float64, cs_prime::Float64, r::Float64, h::Float64)::Float64 = (-2*exp(ν0)*(h^2*(m^2 - 16*π^2*r^6*p^2 - m*(r + 8*π*r^3*p)) + cs^2*(h^2*r^2 + r^4 + (5*h^2 + 4*r^2)*m^2 + m*(-5*h^2*r - 4*r^3 + 8*h^2*π*r^3*ε0) + 2*h^2*π*r^4*(-ε0 + p*(-1 + 8*π*r^2*ε0))) + h^2*r*cs*(r - 2*m)*(-2*r + 5*m + 4*π*r^3*p)*cs_prime))/(h^2*r^3*(r - 2*m))
Ai_plus1(m::Float64, p::Float64, ε0::Float64, ν0::Float64, cs::Float64, cs_prime::Float64, r::Float64, h::Float64)::Float64 = (exp(ν0)*(-(h*(m + 4*π*r^3*p)) + cs^2*(-((5*h + 4*r)*m) + 2*r*(h + r + 2*h*π*r^2*ε0)) + 2*h*r*cs*(r - 2*m)*cs_prime))/(2.0*h^2*r^2)

# in the final row, enforce Δp = 0 at surface. Get an expression like aξ + bξ' = 0. Discretize this with centered differences, solve for ξ[R+h] and eliminate from equation discretized
# at final grid point, which is also discretized with centered differences.
An_minus1(m::Float64, p::Float64, ε0::Float64, ν0::Float64, cs::Float64, cs_prime::Float64, r::Float64, h::Float64)::Float64 = (2*exp(ν0)*cs^2*(r - 2*m))/(h^2*r)
Ann(m::Float64, p::Float64, ε0::Float64, ν0::Float64, cs::Float64, cs_prime::Float64, r::Float64, h::Float64)::Float64 = -((exp(ν0)*(h^2*(7*m^2 - 8*π*r^4*p*(1 + 2*π*r^2*p) + m*(-4*r + 8*π*r^3*p)) + cs^2*((35*h^2 + 20*h*r + 8*r^2)*m^2 + 2*r*m*(-15*h^2 - 9*h*r - 4*r^2 + 2*h*π*r^2*((5*h + 4*r)*p - h*ε0)) + 2*r^2*(3*h^2 + 2*h*r + r^2 + 2*h*π*r^2*(h*ε0 + p*(-3*h - 2*r + 4*h*π*r^2*ε0))))))/(h^2*r^3*(r - 2*m)))

# TOV functions
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

function get_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)
    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    matrix_fname = FrequencyDomain.matrix_fname(star, nPoints); 
    FrequencyDomain.PerfectFluid.Matrix.get_eigensystem(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, N_eigvals; fname=matrix_fname, save_to_file=true);
end

function get_eigensystem(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, N_eigvals::Int64; fname::String = "", save_to_file::Bool=false)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    matrix_size = length(m_arr) - 1;
    A = zeros(matrix_size, matrix_size); # create a square matrix of
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h);
    eigsyst = eigen(A);
    eig_vals = -eigsyst.values;
    eig_vecs = eigsyst.vectors;

    perm = sortperm(eig_vals);
    omega_squared = eig_vals[perm][1:N_eigvals];
    vecs = eig_vecs[:, perm][:, 1:N_eigvals];
    # don't solve for ξ at center so prepend this to arrays
    evecs = vcat(zeros(size(vecs, 2))', vecs)
    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Filename must be provided to save eigensystem."))
        end
        save_eigensystem(fname, omega_squared / (2π * kHz_to_km)^2, evecs, r)
    else
        return omega_squared / (2π * kHz_to_km)^2, evecs, r
    end
end

function save_eigensystem(fname::String, eigvals::Vector{Float64}, eigvecs::AbstractArray{Float64}, radius::Vector{Float64})
    h5open(fname, "w") do file
        file["eigvals"] = eigvals;
        file["eigvecs"] = eigvecs;
        file["radius"] = radius;
    end
end

function load_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64)
    fname = FrequencyDomain.matrix_fname(star, nPoints)
    h5f = h5open(fname, "r")
        eigvals = h5f["eigvals"][:];
        eigvecs = h5f["eigvecs"][:,:];
        radius = h5f["radius"][:];
    close(h5f)
    return eigvals, eigvecs, radius
end

function get_eigenvalues(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)

    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    return FrequencyDomain.PerfectFluid.Matrix.get_eigenvalues(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, N_eigvals);
end

# intended for computing eigenvalues as input to shooting method
function get_eigenvalues(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, N_eigvals::Int64)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    matrix_size = length(m_arr) - 1;
    A = zeros(matrix_size, matrix_size); # create a square matrix of
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h);
    eig_vals = -eigvals(A);

    perm = sortperm(eig_vals);
    omega_squared = eig_vals[perm][1:N_eigvals];

    return omega_squared / (2π * kHz_to_km)^2
end

# fill matrix used in matrix method
function fill_matrix!(A::AbstractArray{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, h::Float64)
    # A is an n x n matrix. Row i of matrix A corresponds to the pulsaiton equation discretized at radius r[i+1]. Note that this index shifting is because the equation is solved at r=0 since ξ(0) = 0 is known.
    n = length(m)-1

    # first row (discretized eqn at r[2], using ξ(0) = 0)
    A[1, 1] = Aii(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2], h)
    A[1, 2] = Ai_plus1(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2], h)

    # last row
    # enforcing lagrangian pressure variation Δp = 0 at surface
    A[end, end-1] = An_minus1(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)
    A[end, end] = Ann(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)

    # # ghost BCs
    # A[end, end-1] = Ai_minus1(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)
    # A[end, end] = Aii(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)

    # interior points
    for i in 2:(n - 1)
        mi = m[i+1]; pressure_i = p[i+1]; εi = ε[i+1]; νi = ν[i+1]; csi = cs[i+1]; cs_prime_i = cs_prime[i+1]; ri = r[i+1];
        A[i, i-1] = Ai_minus1(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, h);
        A[i, i] = Aii(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, h);
        A[i, i+1] = Ai_plus1(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, h);
    end
end

# intended for non-polytopic EOS
function get_frequencies(m::Function, p::Function, ε::Function, ν::Function, r::Vector{Float64}, target_TOV_length::Int64)
    TOV_length = length(m);
    ds = floor(Int, (TOV_length - 1) / (target_TOV_length - 1)); # downsample the TOV solution to match the pulsation solution

    m_ds = m[1:ds:end];
    p_ds = p[1:ds:end];
    ε_ds = ε[1:ds:end];
    ν_ds = ν[1:ds:end];
    r_ds = r[1:ds:end];
    TOV_ds_length = length(m_ds);
    h_ds = diff(r_ds)[1]; # step size for the downsampled TOV solution
    isapprox(h_ds * ones(TOV_ds_length-1), diff(r_ds), atol = 1e-12) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant. Maximum difference: $(maximum(abs.(diff(r_ds) .- h_ds)))"))

    dp = zero(p_ds);
    dε = zero(ε_ds);

    FiniteDiffOrder2.compute_first_derivative(dp, p_ds, h_ds, TOV_ds_length);
    FiniteDiffOrder2.compute_first_derivative(dε, ε_ds, h_ds, TOV_ds_length);

    cs_ds = @. sqrt(dp / dε); # sound speed squared
    cs_prime_ds = zero(cs_ds);
    FiniteDiffOrder2.compute_first_derivative(cs_prime_ds, cs_ds, h_ds, TOV_ds_length);    

    A = zeros(TOV_ds_length-1, TOV_ds_length-1); # create a square matrix of
    fill_matrix!(A, m_ds, p_ds, ε_ds, ν_ds, cs_ds, cs_prime_ds, r_ds, h_ds);
    omega_squared = -eigvals(A) |> sort
    return omega_squared / (2π * kHz_to_km)^2
end

end



module Shoot
using LinearAlgebra
using Dierckx
using HDF5
using NonlinearSolve
using StaticArrays
using NeutronStarOscillations
using ...FrequencyDomain
using ...PerfectFluid
using Printf

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;

function get_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64, h_TOV::Float64, N_eigvals::Int64; nPointsMatrix::Int64, print_progress::Bool=false)
    # extract frequencies from matrix method
    matrix_freqs = NeutronStarOscillations.FrequencyDomain.PerfectFluid.Matrix.get_eigenvalues(star, nPointsMatrix, h_TOV, N_eigvals)
    ω_init = matrix_freqs * (2π * kHz_to_km)^2

    # compute high resolution TOV solution and interpolate to desired computational grid
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r = range(start = 0.0, stop = Rs, step = h_shoot) |> collect;

    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    m = @. m_spline(r);
    p = @. p_spline(r);
    ε = @. ε_spline(r);
    ν = @. ν_spline(r);
    cs = @. sqrt(star.dp_dε(ε_spline(r)));
    cs_prime = @. FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(m, p, ε, cs, r, star.d2p_dε2);
    cs_prime[1] = 0.0;

    # compute eigenvalues and eigenvectors via matrix + shooting method
    fname = FrequencyDomain.shoot_fname(star, h_shoot);
    PerfectFluid.Shoot.iterate_frequencies(m, p, ε, ν, cs, cs_prime, r, ω_init, star.NL_reltol, star.NL_abstol, star.NL_maxiter; fname=fname, save_to_file=true, print_progress = print_progress);
end

function iterate_frequencies(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω_init::Vector{<:Number}, NL_reltol::Float64, NL_abstol::Float64, NL_maxiter::Int64; fname::String = "", save_to_file::Bool=false, print_progress::Bool=false)
    # now use root finding to solve for ω with initial condition given by (low resolution) matrix method
    N_eigvals = length(ω_init);
    ω = zero(ω_init)

    TOV_length = length(r); # minus four comes from removing the additional four steps in the background solution that allows us to evolve the pulsation equations until the surface of the star (since RK4 requires background solution at r + h/2 and r + h)
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    evecs = zeros(ComplexF64, pulsation_length, N_eigvals)
    resids = zero(ω_init)
    retcodes = String[];

    rr = [];

    function bivariate_func(dx, x, params)
        ω0 = x[1]
        ξ_prime_BC = fast_integrate(m, p, ε, ν, cs, cs_prime, r, ω0; shift=1e-18);
        dx[1] = ξ_prime_BC
        return nothing
    end

    for i in 1:N_eigvals
        print_progress ? println("Root finding eigenvalue $i / $N_eigvals") : nothing
        u0 = [real(ω_init[i])]
        prob = NonlinearProblem(bivariate_func, u0)
        sol = solve(prob, NeutronStarOscillations.FrequencyDomain.NL_solver, termination_condition = NeutronStarOscillations.FrequencyDomain.NL_termination_condition,
        reltol = NL_reltol, maxiters = NL_maxiter)
        ω[i] = sol.u[1]

        # now compute eigenvector and boundary condition residual
        rr, ξ, X, ξ_prime_BC = integrate(m, p, ε, ν, cs, cs_prime, r, ω[i]; shift = 1e-18)
        evecs[:, i] .= ξ
        resids[i] = ξ_prime_BC
        push!(retcodes, string(sol.retcode))
    end
    print_progress ? println("Computed all eigenvalues") : nothing
    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Must provide filename to save eigenvalues"))
        else
            save_eigensystem(fname, ω_init ./ (2π * kHz_to_km)^2, ω ./ (2π * kHz_to_km)^2, rr, evecs, resids, retcodes)
        end
    else
        return ω_init ./ (2π * kHz_to_km)^2, ω ./ (2π * kHz_to_km)^2, rr, evecs, resids, retcodes
    end
end

function save_eigensystem(fname::String, ω_init::Vector{<:Number}, shoot_eigvals::Vector{<:Number}, radius::Vector{Float64}, evecs::LinearAlgebra.Matrix{<:Number}, resids::Vector{<:Number}, retcodes::Vector{String})
    h5open(fname, "w") do file
        file["ω_init"] = ω_init;
        file["shoot_eigvals"] = shoot_eigvals;
        file["radius"] = radius;
        file["eigvecs"] = evecs;
        file["resids"] = resids;
        file["retcodes"] = retcodes;
    end;
end

function load_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64)
    fname = FrequencyDomain.shoot_fname(star, h_shoot);
    h5f = h5open(fname, "r")
        ω_init = h5f["ω_init"][:];
        shoot_eigvals = h5f["shoot_eigvals"][:];
        radius = h5f["radius"][:];
        evecs = h5f["eigvecs"][:, :];
        resids = h5f["resids"][:];
        retcodes = h5f["retcodes"][:];
    close(h5f)
    return ω_init, shoot_eigvals, radius, evecs, resids, retcodes
end

# write ODE as A2(r) ξ''(r) + A1(r) ξ'(r) + A0(r) ξ(r) = ω^2 ξ(r)
A0(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = (2*exp(ν)*(m^2 - 16*π^2*r^6*p^2 - m*(r + 8*π*r^3*p) + cs^2*(r^2 - 5*r*m + 5*m^2 - 2*π*r^4*p + 2*π*r^3*(-r + 4*m + 8*π*r^3*p)*ε) + r*cs*(r - 2*m)*(-2*r + 5*m + 4*π*r^3*p)*cs_prime))/(r^3*(r - 2*m))
A1(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = (exp(ν)*(m + 4*π*r^3*p + cs^2*(5*m - 2*(r + 2*π*r^3*ε)) - 2*r*cs*(r - 2*m)*cs_prime))/r^2
A2(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64)::Float64 = -((exp(ν)*cs^2*(r - 2*m))/r)


# RK4 equations. Integrate second-order ODE in space. To reduce to first-order, define X ≡ ξ'[r], ω2 = ω^2
function X_prime(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, ξ, X, ω2)
    A1 = Shoot.A0(m, p, ε, ν, cs, cs_prime, r);
    A2 = Shoot.A1(m, p, ε, ν, cs, cs_prime, r);
    A3 = Shoot.A2(m, p, ε, ν, cs, cs_prime, r);
    return (-(A2*X) + (-A1 + ω2)*ξ)/A3
end
ξ_prime(X) = X

# homogeneous boundary condition to enforce Lagrangian variation of pressure Δp = 0 at surface. Nonlinear rootfinder finds ω^2 such that this quantity is zero.
ξ_prime_BC(m::Float64, p::Float64, r::Float64, ξ) = ((-2*r + 5*m + 4*π*r^3*p)*ξ)/(r*(r - 2*m))

# TOV equations
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function fast_integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω2; shift::Float64=1e-18)
    h = diff(r)[1]; # step size for the TOV solution
    H = 2.0h;

    # initial conditions
    ξi = 0.0
    Xi = 1e-10

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(Number, 2);
    k2 = zeros(Number, 2);
    k3 = zeros(Number, 2);
    k4 = zeros(Number, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξi;
        X0 = Xi;
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, ξ0, X0, ω2) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0
    
        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, ξ1, X1, ω2);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, ξ2, X2, ω2);


        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, ξ3, X3, ω2);

        # update values
        ξi = ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        Xi = X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    p0 = p[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξi;
    X0 = Xi;

    BC_residual = X0 - ξ_prime_BC(m0, p0, r0, ξ0)

    return BC_residual
end

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, ω2; shift::Float64=1e-18)
    h = diff(r)[1];
    H = 2.0h;

    # initial conditions
    ξ = [0.0];
    X = [1e-10];
    rr = [r[1]];

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(ComplexF64, 2);
    k2 = zeros(ComplexF64, 2);
    k3 = zeros(ComplexF64, 2);
    k4 = zeros(ComplexF64, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξ[end];
        X0 = X[end];
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, ξ0, X0, ω2) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0

        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, ξ1, X1, ω2);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, ξ2, X2, ω2);

        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, ξ3, X3, ω2);

        # update values
        push!(rr, r0 + H);
        push!(ξ, ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]));
        push!(X, X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]));

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    p0 = p[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξ[end];
    X0 = X[end];

    BC_residual = X0 - ξ_prime_BC(m0, p0, r0, ξ0)

    return rr, ξ, X, BC_residual
end

end

end

module Eckart
module Matrix
using LinearAlgebra
using Dierckx
using HDF5
using NeutronStarOscillations
using ...FrequencyDomain

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;


function get_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)

    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    matrix_fname = FrequencyDomain.matrix_fname(star, nPoints); 
    ηTimesL = star.η * star.L;
    ζTimesL = star.ζ * star.L;
    FrequencyDomain.Eckart.Matrix.get_eigensystem(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, ηTimesL, ζTimesL, N_eigvals; fname=matrix_fname, save_to_file=true);
end

function get_eigenvalues(star::NeutronStarOscillations.Star, nPoints::Int64, h_TOV::Float64, N_eigvals::Int64)

    # first run TOV solution at higher resolution than h
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r_matrix = range(start = 0.0, stop = Rs, length = nPoints) |> collect;

    # interpolate TOV solution
    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    mInterp(r::Float64)::Float64 = m_spline(r);
    pInterp(r::Float64)::Float64 = p_spline(r);
    εInterp(r::Float64)::Float64 = ε_spline(r);
    νInterp(r::Float64)::Float64 = ν_spline(r);
    csInterp(r::Float64)::Float64 = sqrt(star.dp_dε(εInterp(r)));
    csPrimeInterp(r::Float64)::Float64 = FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(mInterp(r), pInterp(r), εInterp(r), csInterp(r), r, star.d2p_dε2);

    matrix_fname = FrequencyDomain.matrix_fname(star, nPoints); 
    ηTimesL = star.η * star.L;
    ζTimesL = star.ζ * star.L;
    return FrequencyDomain.Eckart.Matrix.get_eigenvalues(mInterp, pInterp, εInterp, νInterp, csInterp, csPrimeInterp, r_matrix, ηTimesL, ζTimesL, N_eigvals);
end

function get_eigensystem(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, η::Float64, ζ::Float64, N_eigvals::Int64; fname::String = "", save_to_file::Bool=false)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    N = length(m_arr) - 1;
    A = zeros(2 * N, 2 * N);

    # upper right block of A = -i
    A[1:N,(N+1):end] .= LinearAlgebra.Matrix(-1.0I, N, N)

    # fill matrix
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h, η, ζ);
    
    # compute eigensystem
    eigsyst = eigen(A);
    eig_vals = -im * eigsyst.values; # eigenvalues are λ = iω, so ω = -iλ
    eig_vecs = eigsyst.vectors;
    omega_filtered_idx = filter_and_sort(eig_vals; atol=atol_pf); # prune the eigenvalues

    omega_filtered = eig_vals[omega_filtered_idx];
    eigvecs_filtered = eig_vecs[:, omega_filtered_idx];
    perm = sortperm(omega_filtered, by=imag, rev=true); # sort by imaginary part in descending order
    omega = omega_filtered[perm][1:N_eigvals] ./ (2π * kHz_to_km);
    vecs = eigvecs_filtered[:, perm][1:N, 1:N_eigvals]; # [1:N] because only the first N entries correspond to ξ
    # don't solve for ξ at center so prepend this to arrays
    evecs = vcat(zeros(size(vecs, 2))', vecs)
    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Filename must be provided to save eigenvalues."))
        end
        save_eigensystem(fname, omega, evecs, r)
    else
        return omega, evecs, r
    end
end

function save_eigensystem(fname::String, filtered_eigvals::Vector{ComplexF64}, evecs::LinearAlgebra.Matrix{ComplexF64}, r::Vector{Float64})
    h5open(fname, "w") do file
        file["eigvals"] = filtered_eigvals;
        file["eigvecs"] = evecs;
        file["radius"] = r;
    end
end

function load_eigensystem(star::NeutronStarOscillations.Star, nPoints::Int64)
    fname = FrequencyDomain.matrix_fname(star, nPoints)
    h5f = h5open(fname, "r")
        eigvals = h5f["eigvals"][:];
        eigvecs = h5f["eigvecs"][:,:];
        radius = h5f["radius"][:];
    close(h5f)
    return eigvals, eigvecs, radius
end

function get_eigenvalues(m::Function, p::Function, ε::Function, ν::Function, cs::Function, cs_prime::Function, r::Vector{Float64}, η::Float64, ζ::Float64, N_eigvals::Int64; fname::String = "", save_to_file::Bool=false)
    m_arr = @. m(r);
    p_arr = @. p(r);
    ε_arr = @. ε(r);
    ν_arr = @. ν(r);
    cs_arr = @. cs(r);
    cs_prime_arr = @. cs_prime(r);
    
    h = diff(r)[1]; # step size for the downsampled TOV solution
    isapprox(h * ones(length(r)-1), diff(r)) || throw(ArgumentError("The radial step size in the downsampled TOV solution is not constant."))

    N = length(m_arr) - 1;
    A = zeros(2 * N, 2 * N);

    # upper right block of A = -i
    A[1:N,(N+1):end] .= LinearAlgebra.Matrix(-1.0I, N, N)

    # fill matrix
    fill_matrix!(A, m_arr, p_arr, ε_arr, ν_arr, cs_arr, cs_prime_arr, r, h, η, ζ);

    # compute and pruneeigenvalues
    eig_vals = -im * eigvals(A);
    omega_filtered_idx = filter_and_sort(eig_vals; atol=atol_pf); # prune the eigenvalues
    omega_filtered = eig_vals[omega_filtered_idx];
    perm = sortperm(omega_filtered, by=imag, rev=true); # sort by imaginary part in descending order
    omega = omega_filtered[perm][1:N_eigvals] ./ (2π * kHz_to_km);
    
    return omega
end


#=

    Here we detail our conventions. We have two equations we are discretizing:

    -Ξ(r) = iωξ(r) [trivial time reduction equaiton Ξ = ∂_{t}ξ in the frequency domain]
    C2 ξ''(r) + C1 ξ'(r) + C0 ξ(r) + B2 Ξ''(r) + B1 Ξ'(r) + B0 Ξ(r) = iωΞ(r)

    Let x := [ξ(r=0), ξ(r=h), ξ(r=2h), ..., ξ(r=R), Ξ(r=0), Ξ(r=h), Ξ(r=2h), ..., Ξ(r=R)]. We thus have the following eigenvalue problem

            Ax = iωx,

    where A is a 2n x 2n matrix, where n is the number of grid points. Let A be composed of blocks such that A = [A11 A12; A21 A22], where A11, A12, A21 and A22 are n x n matrices.

    By inspection, A11=0 and A12 = -I. A21 and A22 are filled with the coefficients of the ξ and Ξ terms, respectively. In the functions below, we therefore fill only the final n rows of the matrix A.
    Note that the solution at r=0 is Ξ(0) = ξ(0) = 0 due to regularity, so n here is the number of grid points minus one.

=#

D2ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = -((exp(ν)*(r - 2*m)*cs^2)/r)
D1ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = (exp(ν)*(4*π*r^3*p - 2*r*cs*((1 + 2*π*r^2*ε)*cs + r*cs_prime) + m*(1 + 5*cs^2 + 4*r*cs*cs_prime)))/r^2
D0ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = (-2*exp(ν)*(m^2*(-1 - 5*cs*(cs - 2*r*cs_prime)) + r^2*(-cs^2 + 2*r*(π*r*(8*π*r^2*p^2 + (p + ε - 8*π*r^2*p*ε)*cs^2) + (1 - 2*π*r^2*p)*cs*cs_prime)) + r*m*(1 + 5*cs^2 + r*(8*π*r*(p - ε*cs^2) + (-9 + 8*π*r^2*p)*cs*cs_prime))))/(r^3*(r - 2*m))

D2Ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = -0.3333333333333333*(exp(ν/2.)*(r - 2*m)*(3*ζ + 4*η)*cs^2)/r
D1Ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = (exp(ν/2.)*(3*ζ + 4*η)*(4*π*r^3*p*(1 + cs^2) - 2*r*cs*((1 + 2*π*r^2*ε)*cs + r*cs_prime) + m*(1 + 6*cs^2 + 4*r*cs*cs_prime)))/(3.0*r^2)
D0Ξ_coeff(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64)::Float64 = (exp(ν/2.)*(m^2*(-15*ζ + 4*η + 5*(3*ζ + 4*η)*cs^2 + 4*r*(-15*ζ + 4*η)*cs*cs_prime) + 2*r*m*(2*η*(-1 - 8*(1 + π*r^2*(6*p + ε))*cs^2 - 2*r*(3 + 8*π*r^2*p)*cs*cs_prime) + 3*ζ*(1 - 4*cs^2 + r*(4*π*r*(-3*p + (-3*p + 2*ε)*cs^2) + (9 - 8*π*r^2*p)*cs*cs_prime))) + 2*r^2*(3*ζ*(cs^2 + 2*r*(π*r*(2*p*(1 - 2*π*r^2*p) + (p - 4*π*r^2*p*(p - 2*ε) - ε)*cs^2) + (-1 + 2*π*r^2*p)*cs*cs_prime)) + 4*η*(cs^2 + r*(-2*π*r*p*(1 + 4*π*r^2*p) + 4*π*r*(2*p*(1 - π*r^2*(p - 2*ε)) + ε)*cs^2 + (1 + 4*π*r^2*p)*cs*cs_prime)))))/(3.0*r^3*(r - 2*m))

p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)


# (i,i-1), (i, i), (i, i+1) entries in the matrix A for a discretized second-order ODE D2 f'' + D1 f' + D0 f = 0. Applies to both ξ and Ξ in their respective blocks.
A_block_i_m1(D0::Float64, D1::Float64, D2::Float64, h::Float64) = -0.5*(-2*D2 + h*D1)/(h^2)
A_block_ii(D0::Float64, D1::Float64, D2::Float64, h::Float64) = (-2*D2)/h^2 + D0
A_block_i_p1(D0::Float64, D1::Float64, D2::Float64, h::Float64) = (2*D2 + h*D1)/(2.0*h^2)

# at the final grid point, enforce Δp = 0 at surface. Get an expression like aξ + bξ' = 0. Discretize this with centered differences, solve for ξ[R+h] and eliminate from equation discretized
# at final grid point, which is also discretized with centered differences.
A_block_n_minus1(D0::Float64, D1::Float64, D2::Float64, h::Float64) = (2*D2)/h^2
A_block_nn(D0::Float64, D1::Float64, D2::Float64, m::Float64, p::Float64, r::Float64, h::Float64) = (h^2*r*D0*(r - 2*m) + h^2*D1*(-2*r + 5*m + 4*π*r^3*p) + 2*D2*(-2*h*r - r^2 + 5*h*m + 2*r*m + 4*h*π*r^3*p))/(h^2*r*(r - 2*m))

function fill_matrix!(A::AbstractArray{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, h::Float64, η::Float64, ζ::Float64)
    n = length(m)-1 # A is an 2n x 2n matrix.

    for i in 2:(n - 1)
        mi = m[i+1]; pressure_i = p[i+1]; εi = ε[i+1]; νi = ν[i+1]; csi = cs[i+1]; cs_prime_i = cs_prime[i+1]; ri = r[i+1];

        D0ξ = D0ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
        D1ξ = D1ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
        D2ξ = D2ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

        D0Ξ = D0Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
        D1Ξ = D1Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
        D2Ξ = D2Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

        A[n + i, i-1] = A_block_i_m1(D0ξ, D1ξ, D2ξ, h)
        A[n + i, i] = A_block_ii(D0ξ, D1ξ, D2ξ, h)
        A[n + i, i+1] = A_block_i_p1(D0ξ, D1ξ, D2ξ, h)

        A[n + i, n + i-1] = A_block_i_m1(D0Ξ, D1Ξ, D2Ξ, h)
        A[n + i, n + i] = A_block_ii(D0Ξ, D1Ξ, D2Ξ, h)
        A[n + i, n + i+1] = A_block_i_p1(D0Ξ, D1Ξ, D2Ξ, h)
    end

    # fill first row manually
    # A[1, 1] = Aii(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2], h)
    # A[1, 2] = Ai_plus1(m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2], h)

    i=1;
    mi, pressure_i, εi, νi, csi, cs_prime_i, ri = m[2], p[2], ε[2], ν[2], cs[2], cs_prime[2], r[2]

    D0ξ = D0ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D1ξ = D1ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D2ξ = D2ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

    D0Ξ = D0Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D1Ξ = D1Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D2Ξ = D2Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

    A[n + i, i] = A_block_ii(D0ξ, D1ξ, D2ξ, h)
    A[n + i, i+1] = A_block_i_p1(D0ξ, D1ξ, D2ξ, h)

    A[n + i, n + i] = A_block_ii(D0Ξ, D1Ξ, D2Ξ, h)
    A[n + i, n + i+1] = A_block_i_p1(D0Ξ, D1Ξ, D2Ξ, h)

    # fill last row manually
    # A[end, end-1] = Ai_minus1(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)
    # A[end, end] = Aii(m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end], h)

    i=n;
    mi, pressure_i, εi, νi, csi, cs_prime_i, ri = m[end], p[end], ε[end], ν[end], cs[end], cs_prime[end], r[end]

    D0ξ = D0ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D1ξ = D1ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D2ξ = D2ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

    D0Ξ = D0Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D1Ξ = D1Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);
    D2Ξ = D2Ξ_coeff(mi, pressure_i, εi, νi, csi, cs_prime_i, ri, η, ζ);

    # # ghost BCs
    # A[n + i, i-1] = A_block_i_m1(D0ξ, D1ξ, D2ξ, h)
    # A[n + i, i] = A_block_ii(D0ξ, D1ξ, D2ξ, h)

    # A[n + i, n + i-1] = A_block_i_m1(D0Ξ, D1Ξ, D2Ξ, h)
    # A[n + i, n + i] = A_block_ii(D0Ξ, D1Ξ, D2Ξ, h)

    # enforcing lagrangian pressure variation Δp = 0 at surface
    A[n + i, i-1] = A_block_n_minus1(D0ξ, D1ξ, D2ξ, h)
    A[n + i, i] = A_block_nn(D0ξ, D1ξ, D2ξ, mi, pressure_i, ri, h)

    A[n + i, n + i-1] = A_block_n_minus1(D0Ξ, D1Ξ, D2Ξ, h)
    A[n + i, n + i] = A_block_nn(D0Ξ, D1Ξ, D2Ξ, mi, pressure_i, ri, h)
end


const atol_pf = 1e-10; # absolute tolerance for deciding whether two eigenvalues are approximately equal in real and imaginary parts
# remove duplicate eigenvalues which have same real part up to a sign. then order by imaginary part descending
function filter_and_sort(z::Vector{ComplexF64}; atol::Float64=1e-8)
    keep_idx = Int64[]
    
    for (i, num) in enumerate(z)
        if real(num) < 0
            # Look for approximate positive counterpart
            has_pair = any(isapprox(real(num), -real(other); atol=atol) &&
                           isapprox(imag(num), imag(other); atol=atol)
                           for other in z)
            if has_pair
                continue  # skip this negative-real member
            end
        end
        push!(keep_idx, i)
    end
    
    return keep_idx
end

end

module Shoot
using LinearAlgebra
using Dierckx
using HDF5
using NonlinearSolve
using StaticArrays
using NeutronStarOscillations
using ...FrequencyDomain
using ...Eckart
using ...PerfectFluid
using Printf

const c = 299792458;
const sec_to_km = c * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;


function get_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64, h_TOV::Float64, N_eigvals::Int64; nPointsMatrix::Int64, print_progress::Bool=false)
    # extract frequencies from matrix method
    matrix_freqs = NeutronStarOscillations.FrequencyDomain.Eckart.Matrix.get_eigenvalues(star, nPointsMatrix, h_TOV, N_eigvals)
    ω_init = matrix_freqs * (2π * kHz_to_km)

    # compute high resolution TOV solution and interpolate to desired computational grid
    r_TOV, m_TOV, p_TOV, ε_TOV, ν_TOV = NeutronStarOscillations.TOV.Explicit.solve(star, h_TOV; save_to_file=false);

    Rs = r_TOV[end];
    r = range(start = 0.0, stop = Rs, step = h_shoot) |> collect;

    m_spline = Spline1D(r_TOV, m_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    p_spline = Spline1D(r_TOV, p_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ε_spline = Spline1D(r_TOV, ε_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    ν_spline = Spline1D(r_TOV, ν_TOV; k=FrequencyDomain.spline_order, s=FrequencyDomain.s, bc=FrequencyDomain.bc);
    m = @. m_spline(r);
    p = @. p_spline(r);
    ε = @. ε_spline(r);
    ν = @. ν_spline(r);
    cs = @. sqrt(star.dp_dε(ε_spline(r)));
    cs_prime = @. FrequencyDomain.PerfectFluid.Matrix.cs_prime_func(m, p, ε, cs, r, star.d2p_dε2);
    cs_prime[1] = 0.0;

    # compute eigenvalues and eigenvectors via matrix + shooting method
    fname = FrequencyDomain.shoot_fname(star, h_shoot);
    Eckart.Shoot.iterate_frequencies(m, p, ε, ν, cs, cs_prime, r, star.η, star.ζ, star.L, ω_init, star.NL_reltol, star.NL_abstol, star.NL_maxiter; fname=fname, save_to_file=true, print_progress = print_progress);
end

function iterate_frequencies(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, L::Float64, ω_init::Vector{<:Number}, NL_reltol::Float64, NL_abstol::Float64, NL_maxiter::Int64; fname::String = "", save_to_file::Bool=false, print_progress::Bool=false)
    # forgot to include length scale in functions for Eckart fluid, so the inputs η and ζ are actually η * L and ζ * L
    ηTimesL = η * L
    ζTimesL = ζ * L

    # now use root finding to solve for ω with initial condition given by (low resolution) matrix method
    N_eigvals = length(ω_init);
    ω = zero(ω_init)

    TOV_length = length(r); # minus four comes from removing the additional four steps in the background solution that allows us to evolve the pulsation equations until the surface of the star (since RK4 requires background solution at r + h/2 and r + h)
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    evecs = zeros(ComplexF64, pulsation_length, N_eigvals)
    resids = zero(ω_init)
    retcodes = String[];

    rr = [];

    function bivariate_func(dx, x, params)
        ω0 = x[1] + im * x[2]
        ξ_prime_BC = fast_integrate(m, p, ε, ν, cs, cs_prime, r, ηTimesL, ζTimesL, ω0; shift=1e-18);
        dx[1] = real(ξ_prime_BC)
        dx[2] = imag(ξ_prime_BC)
        return nothing
    end

    for i in 1:N_eigvals
        print_progress ? println("Root finding eigenvalue $i / $N_eigvals") : nothing
        u0 = [real(ω_init[i]), imag(ω_init[i])]
        prob = NonlinearProblem(bivariate_func, u0)
        sol = solve(prob, NeutronStarOscillations.FrequencyDomain.NL_solver, termination_condition = NeutronStarOscillations.FrequencyDomain.NL_termination_condition,
        reltol = NL_reltol, maxiters = NL_maxiter)
        ω[i] = sol.u[1] + im * sol.u[2]

        # now compute eigenvector and boundary condition residual
        rr, ξ, X, ξ_prime_BC = integrate(m, p, ε, ν, cs, cs_prime, r, ηTimesL, ζTimesL, ω[i]; shift = 1e-18)
        evecs[:, i] .= ξ
        resids[i] = ξ_prime_BC
        push!(retcodes, string(sol.retcode))
    end
    print_progress ? println("Computed all eigenvalues") : nothing
    FrequencyDomain.normalize_eigvecs!(evecs)
    if save_to_file
        if fname == ""
            throw(ArgumentError("Must provide filename to save eigenvalues"))
        else
            save_eigensystem(fname, ω_init ./ (2π * kHz_to_km), ω ./ (2π * kHz_to_km), rr, evecs, resids, retcodes)
        end
    else
        return ω_init ./ (2π * kHz_to_km), ω ./ (2π * kHz_to_km), rr, evecs, resids, retcodes
    end
end

function save_eigensystem(fname::String, ω_init::Vector{<:Number}, shoot_eigvals::Vector{<:Number}, radius::Vector{Float64}, evecs::LinearAlgebra.Matrix{<:Number}, resids::Vector{<:Number}, retcodes::Vector{String})
    h5open(fname, "w") do file
        file["ω_init"] = ω_init;
        file["shoot_eigvals"] = shoot_eigvals;
        file["radius"] = radius;
        file["eigvecs"] = evecs;
        file["resids"] = resids;
        file["retcodes"] = retcodes;
    end;
end

function load_eigensystem(star::NeutronStarOscillations.Star, h_shoot::Float64)
    fname = FrequencyDomain.shoot_fname(star, h_shoot);
    h5f = h5open(fname, "r")
        ω_init = h5f["ω_init"][:];
        shoot_eigvals = h5f["shoot_eigvals"][:];
        radius = h5f["radius"][:];
        evecs = h5f["eigvecs"][:, :];
        resids = h5f["resids"][:];
        retcodes = h5f["retcodes"][:];
    close(h5f)
    return ω_init, shoot_eigvals, radius, evecs, resids, retcodes
end


C1(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = -((exp(ν)*(r - 2*m)*cs^2)/r)
C2(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = (exp(ν)*(4*π*r^3*p - 2*r*cs*((1 + 2*π*r^2*ε)*cs + r*cs_prime) + m*(1 + 5*cs^2 + 4*r*cs*cs_prime)))/r^2
C3(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = (-2*exp(ν)*(m^2*(-1 - 5*cs*(cs - 2*r*cs_prime)) + r^2*(-cs^2 + 2*r*(π*r*(8*π*r^2*p^2 + (p + ε - 8*π*r^2*p*ε)*cs^2) + (1 - 2*π*r^2*p)*cs*cs_prime)) + r*m*(1 + 5*cs^2 + r*(8*π*r*(p - ε*cs^2) + (-9 + 8*π*r^2*p)*cs*cs_prime))))/(r^3*(r - 2*m))
B1(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = -0.3333333333333333*(exp(ν/2.)*(r - 2*m)*(3*ζ + 4*η)*cs^2)/r
B2(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = (exp(ν/2.)*(3*ζ + 4*η)*(4*π*r^3*p*(1 + cs^2) - 2*r*cs*((1 + 2*π*r^2*ε)*cs + r*cs_prime) + m*(1 + 6*cs^2 + 4*r*cs*cs_prime)))/(3.0*r^2)
B3(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64) = (exp(ν/2.)*(m^2*(-15*ζ + 4*η + 5*(3*ζ + 4*η)*cs^2 + 4*r*(-15*ζ + 4*η)*cs*cs_prime) + 2*r*m*(2*η*(-1 - 8*(1 + π*r^2*(6*p + ε))*cs^2 - 2*r*(3 + 8*π*r^2*p)*cs*cs_prime) + 3*ζ*(1 - 4*cs^2 + r*(4*π*r*(-3*p + (-3*p + 2*ε)*cs^2) + (9 - 8*π*r^2*p)*cs*cs_prime))) + 2*r^2*(3*ζ*(cs^2 + 2*r*(π*r*(2*p*(1 - 2*π*r^2*p) + (p - 4*π*r^2*p*(p - 2*ε) - ε)*cs^2) + (-1 + 2*π*r^2*p)*cs*cs_prime)) + 4*η*(cs^2 + r*(-2*π*r*p*(1 + 4*π*r^2*p) + 4*π*r*(2*p*(1 - π*r^2*(p - 2*ε)) + ε)*cs^2 + (1 + 4*π*r^2*p)*cs*cs_prime)))))/(3.0*r^3*(r - 2*m))

# RK4 equations. Integrate second order ODE in space. To reduce to first order, define X ≡ ξ'[r]
X_prime(m::Float64, p::Float64, ε::Float64, ν::Float64, cs::Float64, cs_prime::Float64, r::Float64, η::Float64, ζ::Float64, ξ, X, ω) = (-((C2(m, p, ε, ν, cs, cs_prime, r, η, ζ) - im*B2(m, p, ε, ν, cs, cs_prime, r, η, ζ)*ω)*X) + (-C3(m, p, ε, ν, cs, cs_prime, r, η, ζ) + ω*(im*B3(m, p, ε, ν, cs, cs_prime, r, η, ζ) + ω))*ξ)/(C1(m, p, ε, ν, cs, cs_prime, r, η, ζ) - im*B1(m, p, ε, ν, cs, cs_prime, r, η, ζ)*ω)
ξ_prime(X) = X

# homogeneous boundary condition to enforce Lagrangian variation of pressure Δp = 0 at surface. Nonlinear rootfinder finds ω^2 such that this quantity is zero.
ξ_prime_BC(m::Float64, p::Float64, r::Float64, ξ) = ((-2*r + 5*m + 4*π*r^3*p)*ξ)/(r*(r - 2*m))

# TOV equations
p_prime(m::Float64, p::Float64, ε::Float64, r::Float64)::Float64 = -0.5*((-1 + r/(r - 2*m) + (8*π*r^3*p)/(r - 2*m))*(p + ε))/r
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Function)::Float64 = d2p_dε2(ε) * p_prime(m, p, ε, r) / (2 * cs^3)
cs_prime_func(m::Float64, p::Float64, ε::Float64, cs::Float64, r::Float64, d2p_dε2::Float64)::Float64 = d2p_dε2 * p_prime(m, p, ε, r) / (2 * cs^3)

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function fast_integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, ω; shift::Float64=1e-18)
    h = diff(r)[1]; # step size for the TOV solution
    H = 2.0h;

    # initial conditions
    ξi = 0.0+0.0im;
    Xi = 1e-10+0.0im;

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(Number, 2);
    k2 = zeros(Number, 2);
    k3 = zeros(Number, 2);
    k4 = zeros(Number, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξi;
        X0 = Xi;
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, η, ζ, ξ0, X0, ω) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0
    
        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, η, ζ, ξ1, X1, ω);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, η, ζ, ξ2, X2, ω);


        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, η, ζ, ξ3, X3, ω);

        # update values
        ξi = ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        Xi = X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    p0 = p[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξi;
    X0 = Xi;

    BC_residual = X0 - ξ_prime_BC(m0, p0, r0, ξ0)

    return BC_residual
end

# in this function we take as argument the TOV solution with radial step h. RK4 requires background solution at half radial steps, so to ensure this can be done we will solve the pulsation equation at radial steps of 2h
function integrate(m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, cs_prime::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, ω; shift::Float64=1e-18)
    h = diff(r)[1];
    H = 2.0h;

    # initial conditions
    ξ = [0.0+0.0im];
    X = [1e-10+0.0im];
    rr = [r[1]];

    TOV_length = length(r);
    pulsation_length = floor(Int, (TOV_length - 1) / 2) + 1;

    # [ξ, X]
    k1 = zeros(ComplexF64, 2);
    k2 = zeros(ComplexF64, 2);
    k3 = zeros(ComplexF64, 2);
    k4 = zeros(ComplexF64, 2);

    regularity = 0.0; # used to enforce X'[0] = 0
    for i in 1:pulsation_length-1
        TOV_sol_idx = 2 * i - 1

        # compute k1 (background is at location r[TOV_sol_idx])
        r0 = r[TOV_sol_idx];
        p0 = p[TOV_sol_idx];
        m0 = m[TOV_sol_idx];
        ν0 = ν[TOV_sol_idx];
        ε0 = ε[TOV_sol_idx];
        cs0 = cs[TOV_sol_idx];
        cs0_prime = cs_prime[TOV_sol_idx];

        ξ0 = ξ[end];
        X0 = X[end];
 
        k1[1] = ξ_prime(X0);
        k1[2] = X_prime(m0, p0, ε0, ν0, cs0, cs0_prime, r0+shift, η, ζ, ξ0, X0, ω) * regularity; # effect of shift & regularity is to enforce that X'[0] = 0

        # compute k2 (background is at step r[TOV_sol_idx] + h/2) where h = time step in pulsation solution = 2 * h, where h is the time step in TOV solution
        r1 = r[TOV_sol_idx+1];
        p1 = p[TOV_sol_idx+1];
        m1 = m[TOV_sol_idx+1];
        ν1 = ν[TOV_sol_idx+1];
        ε1 = ε[TOV_sol_idx+1];
        cs1 = cs[TOV_sol_idx+1];
        cs1_prime = cs_prime[TOV_sol_idx+1];

        ξ1 = ξ0 + H/2*k1[1];
        X1 = X0 + H/2*k1[2];

        k2[1] = ξ_prime(X1);
        k2[2] = X_prime(m1, p1, ε1, ν1, cs1, cs1_prime, r1, η, ζ, ξ1, X1, ω);

        # compute k3 (background is at step r[TOV_sol_idx] + h/2)
        r2 = r1;
        p2 = p1;
        m2 = m1;
        ν2 = ν1;
        ε2 = ε1;
        cs2 = cs1;
        cs2_prime = cs1_prime

        ξ2 = ξ0 + H/2*k2[1];
        X2 = X0 + H/2*k2[2];

        k3[1] = ξ_prime(X2);
        k3[2] = X_prime(m2, p2, ε2, ν2, cs2, cs2_prime, r2, η, ζ, ξ2, X2, ω);

        # compute k4 (background is at step r[TOV_sol_idx] + h)
        r3 = r[TOV_sol_idx+2];
        p3 = p[TOV_sol_idx+2];
        m3 = m[TOV_sol_idx+2];
        ν3 = ν[TOV_sol_idx+2];
        ε3 = ε[TOV_sol_idx+2];
        cs3 = cs[TOV_sol_idx+2];
        cs3_prime = cs_prime[TOV_sol_idx+2];

        ξ3 = ξ0 + H*k3[1];
        X3 = X0 + H*k3[2];

        k4[1] = ξ_prime(X3);
        k4[2] = X_prime(m3, p3, ε3, ν3, cs3, cs3_prime, r3, η, ζ, ξ3, X3, ω);

        # update values
        push!(rr, r0 + H);
        push!(ξ, ξ0 + H/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]));
        push!(X, X0 + H/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]));

        # remove regularity control after first time step
        shift = 0.0;
        regularity = 1.0;
    end

    # compute BC residual at the last step
    TOV_sol_idx = 2 * pulsation_length - 1
    r0 = r[TOV_sol_idx];
    p0 = p[TOV_sol_idx];
    m0 = m[TOV_sol_idx];
    ξ0 = ξ[end];
    X0 = X[end];

    BC_residual = X0 - ξ_prime_BC(m0, p0, r0, ξ0)

    return rr, ξ, X, BC_residual
end
end
end
end
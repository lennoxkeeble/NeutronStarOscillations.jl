include("params.jl")
res = 1e-2;
mode = 0;
KO = 0.0;
CFL = 0.05; # Courant-Friedrichs-Lewy factor — (Float64)
total_time_ms = 0.05; # total integration time [ms] — (Float64)
dt_save_ms = total_time_ms / 1000.0; # time interval between saved data points [ms] — (Float64)
save_every = 200; # save to file after 'save_every' time steps have been stored in memory (i.e., after every Δt = save_every * dt_save_ms) — (Int64)

cowling = false;
# star = NeutronStarOscillations.Star(eps_central, kappa, n, 0.0, 0.0, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, mode, cowling, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);

star = NeutronStarOscillations.Star(eps_central, kappa, n, η, ζ, 0.0, 0.0, 0.0, L, ptol, ptol_TD, data_path, fig_path, TOV_iter_tol, TOV_max_iter, TOV_max_steps, TOV_initial_r, NL_reltol, NL_abstol, NL_maxiter, N_eigvals, xi_ID, xi_dt_ID, du_ID, du_dr_ID, du_dt_ID, de_ID, de_dr_ID, de_dt_ID, KO, CFL, dt_save_ms, h_save, save_every, total_time_ms);
@time time_integrate(star, res, cowling; print_progress = true);

var = "xi";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi\,[\mathrm{km}]";

time = 1000 * total_time_ms / 1000.0; plot_td_var(
    star, res, var, time, cowling;
    desample_factor=1,
    annotate_time=true,
    xlabel = xlabel,
    ylabel = ylabel
);

var = "xi";
xlabel = L"r\,[\mathrm{km}]";
ylabel = L"\xi";
animate_td_var(
    star, res, var, cowling;
    desample_factor=1,
    stop_time=star.T, # stop time of animation in ms
    animation_length=10.0, # length of animation in seconds
    framerate=20, # frames per second
    xlabel = xlabel,
    ylabel = ylabel,
    fix_ylims = true, # whether to fix y-limits of plot
    lim_y_min = -1.0,
    lim_y_max = 1.0
)

# frequency domain comparison
include("FFTFunctions.jl")
function extract_td_var_time_series(r::Float64, var::String, star::NeutronStarOscillations.Star, h::Float64)
    sol = load_td_solution(star, h, cowling)
    r_idx = argmin(@. abs(sol["solution/r"][:] - r))
    signal = sol["solution/"*var][:, r_idx]
    close(sol)
    return signal
end

function extract_single_FFT(t::Vector{Float64}, signal::Vector{Float64})
    nPointsFFT = length(signal)
    fs = nPointsFFT / t[end]; # sampling frequency
    F, freqs = FourierFunctions.compute_real_FFT(signal, fs, nPointsFFT);
    return freqs, F
end

function extract_averaged_FFT(num_xvals::Int64, var::String, star::NeutronStarOscillations.Star, h::Float64)
    sol = load_td_solution(star, h, cowling)
    r = sol["solution/r"][:]
    t = sol["solution/t"][:]
    close(sol)

    if length(r) < num_xvals
        xvals = r[:];
    else
        ds = Int(floor((length(r)-2) / num_xvals))
        xvals = r[2:ds:end-1];
    end

    xvals = [4.0]

    signal = extract_td_var_time_series(xvals[1], var, star, h)
    freqs, F = extract_single_FFT(t, signal)

    for i in 2:length(xvals)
        ff, FF = extract_single_FFT(t, signal)
        F .+= FF
    end

    return freqs, abs.(F) / maximum(abs.(F))
end

function get_FFT(star::NeutronStarOscillations.Star, num_xvals::Int64)
    star_type = NeutronStarOscillations.get_star_type(star)
    
    if star_type == "BDNK"
        var = "du"
    else
        var = "xi"
    end
    
    return extract_averaged_FFT(num_xvals, var, star, res)
end

xlabel = L"ω \; [\mathrm{kHz}]";
ylabel = L"|\mathcal{F}|";
labels = ["KO=0", "KO=0.5"];

lim_x_min = 0.0;
lim_x_max = 20.0;
lim_y_min = 1e-4;
lim_y_max = 1e1;

function plot_star_FD(star, num_xvals)
    freqs_1, F_1 = get_FFT(star, num_xvals);
    
    # fd_freqs = load_fd_freqs(param_idx);
    # vlines = real.(fd_freqs);
    vlines = [];

    NeutronStarOscillations.QuickPlots.plot11([freqs_1], [F_1];
        xlabel = xlabel,
        ylabel = ylabel,
        lim_x_min = lim_x_min,
        lim_x_max = lim_x_max,
        lim_y_min = lim_y_min,
        lim_y_max = lim_y_max,
        labels = labels,
        yscale = log10,
        legend = true,
        position = :cb,
        vlines = vlines)
end

num_xvals = 5;
plot_star_FD(star, num_xvals)
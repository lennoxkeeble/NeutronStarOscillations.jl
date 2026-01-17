include("../../params.jl");

# solve TOV equations for a very low pressure tolerance
BDNK_star.ptol = 1e-14;
h = 1e-5;
@time NeutronStarOscillations.TOV.Explicit.solve(BDNK_star, h); # explicit RK4 solver (fourth-order accurate)

# load TOV solution
r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.load(BDNK_star, h);

# extract mass and radius as a function of pressure near the surface
ptol_facts = [10^(-1.0i) for i=2:10];
ptol_vals = @. BDNK_star.pc_SI * ptol_facts;

stellar_radii = zero(ptol_vals);
stellar_masses = zero(ptol_vals);
fitted_stellar_radii = zero(ptol_vals);
fitted_stellar_masses = zero(ptol_vals);

for i in eachindex(ptol_vals)
    idx = argmin(@. abs(p - ptol_vals[i]));
    stellar_radii[i] = r[idx];
    stellar_masses[i] = m[idx];
end

# extract mass and radius from fitting solution near the surface terminated at each pressure tolerance.
# fit to y = m * (R - x)^a for m and R
function fit(x::Vector{Float64}, y::Vector{Float64}, a::Float64)
    # Transform y based on the linearization: y^(1/a) = (m^(1/a)R) - (m^(1/a)) * x
    y_trans = y .^ (1/a)

    # Create the design matrix for linear regression: [1, x]
    A = [ones(length(x)) x]

    # Solve the linear system (Least Squares)
    coeffs = A \ y_trans
    intercept, slope = coeffs[1], coeffs[2]

    # Recover m and R from the slope and intercept
    # Slope is expected to be negative: slope = -m^(1/a)
    prefactor = (-slope)^a
    radius = intercept / (-slope)
    return prefactor, radius
end

function get_fitted_MR(star::NeutronStarOscillations.Star, ptol::Float64, h::Float64; fit_Δr::Float64=20.0h)
    # load solution with very low pressure tolerance
    r, m, p, ε, ν = NeutronStarOscillations.TOV.Explicit.load(star, h);

    # truncate solution at given pressure tolerance
    max_idx = argmin(@. abs(p - ptol));
    r = r[1:max_idx];
    m = m[1:max_idx];
    p = p[1:max_idx];
    ε = ε[1:max_idx];
    cs2 = star.dp_dε.(ε);

    # truncate solution near stellar surface for fitting
    mask = @. (r >= r[end] - fit_Δr);

    r_surface = r[mask];
    ε_surface = ε[mask];
    p_surface = p[mask];
    cs2_surface = cs2[mask];

    # extract mass from surface fits of ε and p
    ε_mass(coeff, radius) = (coeff^(1/star.n)*star.kappa*(1 + star.n)*radius^2)/(1 + 2*coeff^(1/star.n)*star.kappa*radius + 2*coeff^(1/star.n)*star.kappa*star.n*radius)
    p_mass(coeff, radius) = (star.kappa*(1 + star.n)*radius^2*(coeff/star.kappa)^(1/(1 + star.n)))/(1 + 2*star.kappa*radius*(coeff/star.kappa)^(1/(1 + star.n)) + 2*star.kappa*star.n*radius*(coeff/star.kappa)^(1/(1 + star.n)))

    ε_coeff, ε_radius = fit(r_surface, ε_surface, star.n); ε_mass_km = ε_mass(ε_coeff, ε_radius)
    p_coeff, p_radius = fit(r_surface, p_surface, star.n + 1); p_mass_km = p_mass(p_coeff, p_radius)
    cs2_coeff, cs2_radius = fit(r_surface, cs2_surface, 1.0);
    cs_coeff = sqrt(cs2_coeff);

    # check that mass and radius from fits to energy density and pressure agree
    if !isapprox(ε_mass_km, p_mass_km, atol=h)
        error("Masses from ε and p fits do not agree: ε mass = $(ε_mass_km) Msun, p mass = $(p_mass_km) Msun")
    end

    if !isapprox(ε_radius, p_radius, atol=h)
        error("Radii from ε and p fits do not agree: ε radius = $(ε_radius) km, p radius = $(p_radius) km")
    end

    return ε_mass_km, ε_radius
end

for i in eachindex(ptol_vals)
    fitted_mass_km, fitted_radius = get_fitted_MR(BDNK_star, ptol_vals[i], h; fit_Δr=50.0h)
    fitted_stellar_masses[i] = fitted_mass_km
    fitted_stellar_radii[i] = fitted_radius
end

# compute error in mass and radii as a function of pressure tolerance, taking the TOV solution at the lowest pressure tolerance as the benchmark
mass_errors = @. abs(1 - stellar_masses / m[end]);
radius_errors = @. abs(1 - stellar_radii / r[end]);
fitted_mass_errors = @. abs(1 - fitted_stellar_masses / m[end]);
fitted_radius_errors = @. abs(1 - fitted_stellar_radii / r[end]);

# plot errors
xlabel = L"p_{\mathrm{tol}} / p_{\mathrm{c}}"
ylabel = L"\mathrm{Relative\ Error}"
x_plot = [ptol_facts, ptol_facts, ptol_facts, ptol_facts];
y_plot = [mass_errors, radius_errors, fitted_mass_errors, fitted_radius_errors];
linestyles = [:solid, :dash, :solid, :dash];
colors = [:tomato, :tomato, :royalblue, :royalblue];
labels = ["Mass (RK4)", "Radius (RK4)", "Mass (Fitted)", "Radius (Fitted)"];
NeutronStarOscillations.QuickPlots.plot11(
    x_plot, y_plot;
    colors = colors,
    labels = labels,
    linestyles = linestyles,
    linewidths = NeutronStarOscillations.QuickPlots.default_linewidths,
    alphas = NeutronStarOscillations.QuickPlots.default_alphas,
    xlabel = xlabel,
    ylabel = ylabel,
    xscale = log10,
    yscale = log10,
    legend = false,
    position = :rt, # position of legend
    framevisible = true, # whether to draw box around legend
    scatter_lines = false # whether to plot scatter lines instead of normal lines
);
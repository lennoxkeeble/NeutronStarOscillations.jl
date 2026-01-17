include("../../params.jl");
h = 1e-4;
NeutronStarOscillations.TOV.Explicit.solve(BDNK_star, h); # explicit RK4 solver (fourth-order accurate)
NeutronStarOscillations.TOV.Implicit.solve(BDNK_star, h); # implicit Crank Nicholson solver (second-order accurate)

r_implicit, m_implicit, p_implicit, ε_implicit, ν_implicit = NeutronStarOscillations.TOV.Implicit.load(BDNK_star, h);
r_explicit, m_explicit, p_explicit, ε_explicit, ν_explicit = NeutronStarOscillations.TOV.Explicit.load(BDNK_star, h);

atol = 1e-5;
isapprox(r_implicit, r_explicit, atol=atol) && isapprox(m_implicit, m_explicit, atol=atol) && isapprox(p_implicit, p_explicit, atol=atol) && isapprox(ε_implicit, ε_explicit, atol=atol) && isapprox(ν_implicit, ν_explicit, atol=atol)

println("Mass = $(m_explicit[end] / NeutronStarOscillations.Msun_to_km) solar masses, R = $(r_explicit[end]) km")

# plots
vars = ["mass", "pressure", "energy_density", "nu"]; 
xlabel = L"r\,[\mathrm{km}]";
ylabels = [L"m\,[\mathrm{km}]", L"p\,[\mathrm{km}^{-2}]", L"\epsilon\,[\mathrm{km}^{-2}]", L"\nu"];
legend = false;
plot_desample_factor = 100; # desampling factor for plots to reduce number of points plotted

for i in eachindex(vars)
    var = vars[i];
    ylabel = ylabels[i];
    NeutronStarOscillations.TOV.plot_var(var, BDNK_star, h, "Explicit"; 
        desample_factor = plot_desample_factor,
        xlabel = xlabel,
        ylabel = ylabel,
        legend = legend,
    );
end
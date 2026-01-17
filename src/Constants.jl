#=

    Useful physical constants and conversion factors.

=#

const c = 299792458;
const sec_to_cm = c * 1e2;
const sec_to_km = c * 1e-3;
const G = 6.6743015e-11;
const M_sun = 1.988416e30;
const kg_to_m = G / c^2;
const cm_to_km = 1e-5;
const kg_to_g = 1e3;
const gram_per_cm3_to_km_minus2 = kg_to_m * 1e-3 / kg_to_g / cm_to_km^3
const e_charge = 1.602176634e-19;
const kB = 1.380649e-23;

const dyne_per_cm2_to_km_minus2 = gram_per_cm3_to_km_minus2 / sec_to_cm^2
const Msun_to_km = M_sun * kg_to_m * 1e-3;
const kHz_to_km = 1e3 / sec_to_km;
const MeV_to_km = 1e6 * G * c^(-4) * e_charge;
const fm_to_km = 1e-15 * 1e-3;
const MeV_per_fm3_to_dyne_per_cm2 = 1.602176634e33
const MeV_per_fm3_to_km_minus2 = MeV_per_fm3_to_dyne_per_cm2 * dyne_per_cm2_to_km_minus2;
const MeV_to_K = 10^6 * e_charge / kB; # MeV to Kelvin conversion factor
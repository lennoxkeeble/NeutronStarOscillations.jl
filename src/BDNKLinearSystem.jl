#=

    Module comprising of the discretized time evolution equations for the BDNK system.

=#

module BDNKLinearSystem

# equation for v3 ≡ ∂_{t}u3 used to eliminate time derivatives from EFEs and conservation of energy-momentum equations
Forward_v3(i::Int64, u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64)::Float64 = (3*u1[i] - 4*u1[i+1] + u1[2+i])/h - (6*exp(ν[i]/2.)*u2[i])/(L*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (3*exp(ν[i]/2.)*(r[i] - 2*m[i])*(3*u3[i] - 4*u3[i+1] + u3[2 + i]))/(8.0*h*L*π*r[i]^2*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*v2[i])/(p[i] + ε[i]) + (3*exp(ν[i]/2.)*u3[i]*(r[i]^(-2) - 8*π*ε[i]))/(4.0*L*π*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*u1[i]*(2*r[i] - 5*m[i] + 4*π*r[i]^3*ε[i]))/(r[i]*(r[i] - 2*m[i])) + (2*(m[i] + 4*π*r[i]^3*p[i])*u1[i])/(r[i]*(r[i] - 2*m[i])*cs[i]^2)
Backward_v3(i::Int64, u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64)::Float64 = -((u1[i-2] - 4*u1[i-1] + 3*u1[i])/h) - (6*exp(ν[i]/2.)*u2[i])/(L*(3*ζ + 4*η)*τε*(p[i] + ε[i])) + (3*exp(ν[i]/2.)*(r[i] - 2*m[i])*(u3[i-2] - 4*u3[i-1] + 3*u3[i]))/(8.0*h*L*π*r[i]^2*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*v2[i])/(p[i] + ε[i]) + (3*exp(ν[i]/2.)*u3[i]*(r[i]^(-2) - 8*π*ε[i]))/(4.0*L*π*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*u1[i]*(2*r[i] - 5*m[i] + 4*π*r[i]^3*ε[i]))/(r[i]*(r[i] - 2*m[i])) + (2*(m[i] + 4*π*r[i]^3*p[i])*u1[i])/(r[i]*(r[i] - 2*m[i])*cs[i]^2)
Interior_v3(i::Int64, u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64)::Float64 = (u1[i-1] - u1[i+1])/h - (6*exp(ν[i]/2.)*u2[i])/(L*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (3*exp(ν[i]/2.)*(r[i] - 2*m[i])*(u3[i-1] - u3[i+1]))/(8.0*h*L*π*r[i]^2*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*v2[i])/(p[i] + ε[i]) + (3*exp(ν[i]/2.)*u3[i]*(r[i]^(-2) - 8*π*ε[i]))/(4.0*L*π*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*u1[i]*(2*r[i] - 5*m[i] + 4*π*r[i]^3*ε[i]))/(r[i]*(r[i] - 2*m[i])) + (2*(m[i] + 4*π*r[i]^3*p[i])*u1[i])/(r[i]*(r[i] - 2*m[i])*cs[i]^2)

function compute_v3!(u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64, h::Float64)
    N = length(v3)
    for i = 2:N-1
        v3[i] = BDNKLinearSystem.Interior_v3(i, u1, u2, u3, v1, v2, m, p, ε, ν, cs, r, η, ζ, τε, τP, τQ, L, h)
    end

    v3[1] = BDNKLinearSystem.Forward_v3(1, u1, u2, u3, v1, v2, m, p, ε, ν, cs, r, η, ζ, τε, τP, τQ, L, h)
    v3[N] = BDNKLinearSystem.Backward_v3(N, u1, u2, u3, v1, v2, m, p, ε, ν, cs, r, η, ζ, τε, τP, τQ, L, h)
end

# equation for v3 ≡ ∂_{t}u3 used to eliminate time derivatives from EFEs and conservation of energy-momentum equations
v3(i::Int64, u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, w1::Vector{Float64}, w2::Vector{Float64}, w3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = -2*w1[i] - (6*exp(ν[i]/2.)*u2[i])/(L*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*v2[i])/(p[i] + ε[i]) + (3*exp(ν[i]/2.)*(r[i] - 2*m[i])*w3[i])/(4.0*L*π*r[i]^2*(3*ζ + 4*η)*τε*(p[i] + ε[i])) + (3*exp(ν[i]/2.)*u3[i]*(r[i]^(-2) - 8*π*ε[i]))/(4.0*L*π*(3*ζ + 4*η)*τε*(p[i] + ε[i])) - (2*u1[i]*(2*r[i] - 5*m[i] + 4*π*r[i]^3*ε[i]))/(r[i]*(r[i] - 2*m[i])) + (2*(m[i] + 4*π*r[i]^3*p[i])*u1[i])/(r[i]*(r[i] - 2*m[i])*cs[i]^2)

function compute_v3!(u1::Vector{Float64}, u2::Vector{Float64}, u3::Vector{Float64}, w1::Vector{Float64}, w2::Vector{Float64}, w3::Vector{Float64}, v1::Vector{Float64}, v2::Vector{Float64}, v3::Vector{Float64}, m::Vector{Float64}, p::Vector{Float64}, ε::Vector{Float64}, ν::Vector{Float64}, cs::Vector{Float64}, r::Vector{Float64}, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)
    for i in eachindex(r)
        v3[i] = BDNKLinearSystem.v3(i, u1, u2, u3, w1, w2, w3, v1, v2, m, p, ε, ν, cs, r, η, ζ, τε, τP, τQ, L)
    end
end

module Center

Eq(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, h::Float64, k::Float64)::Float64 = (CNSystem.E1*((4*u1n[2] - u1n[3])/(2.0*h) + (4*u1np1[2] - u1np1[3])/(2.0*h)))/2. + (CNSystem.E2*(u2n[1] + u2np1[1]))/2. + (CNSystem.E3*((2*u2n[1] - 5*u2n[2] + 4*u2n[3] - u2n[4])/h^2 + (2*u2np1[1] - 5*u2np1[2] + 4*u2np1[3] - u2np1[4])/h^2))/4. + (CNSystem.E4*((4*v1n[2] - v1n[3])/(2.0*h) + (4*v1np1[2] - v1np1[3])/(2.0*h)))/2. + (CNSystem.E6*(-v2n[1] + v2np1[1]))/k + (CNSystem.E5*(v2n[1] + v2np1[1]))/2.

pCenterEq2pU1_2(CNSystem, h::Float64, k::Float64)::Float64 = CNSystem.E1/h
pCenterEq2pU1_3(CNSystem, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.E1/h

pCenterEq2pU2_1(CNSystem, h::Float64, k::Float64)::Float64 = CNSystem.E2/2. + CNSystem.E3/(2.0*h^2)
pCenterEq2pU2_2(CNSystem, h::Float64, k::Float64)::Float64 = (-5*CNSystem.E3)/(4.0*h^2)
pCenterEq2pU2_3(CNSystem, h::Float64, k::Float64)::Float64 = CNSystem.E3/h^2
pCenterEq2pU2_4(CNSystem, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.E3/h^2

pCenterEq2pV1_2(CNSystem, h::Float64, k::Float64)::Float64 = CNSystem.E4/h
pCenterEq2pV1_3(CNSystem, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.E4/h

pCenterEq2pV2_1(CNSystem, h::Float64, k::Float64)::Float64 = CNSystem.E5/2. + CNSystem.E6/k

end

module RHS
using ..Center

# set up right hand side of linear system. Our convention is that the first matrix_solve_nPoints entries correspond to the trivial time reduction equation for u1 at each spatial point, the next matrix_solve_nPoints entries correspond to the trivial time reduction equation for u2 at each spatial point, the next matrix_solve_nPoints entries correspond to Eq1 at each spatial point, and the last matrix_solve_nPoints entries correspond to Eq2 at each spatial point, where Eq1 and Eq2 are the nontrivial BDNK equations.
function compute_RHS!(linear_system_RHS::Vector{Float64}, u1_np1::Vector{Float64}, u2_np1::Vector{Float64}, u3_np1::Vector{Float64}, v1_np1::Vector{Float64}, v2_np1::Vector{Float64}, u1_n::Vector{Float64}, u2_n::Vector{Float64}, u3_n::Vector{Float64}, v1_n::Vector{Float64}, v2_n::Vector{Float64}, CNSystem, h::Float64, k::Float64, nPointsSpace::Int64)
    # awkward indexing since we solve for (u1, v1, w2, u3, w3) everywhere apart from center but solve for (w1, u2, v2) everywhere including center
    total_trivial_u1_eqs = nPointsSpace-1
    total_trivial_u2_eqs = nPointsSpace
    total_BDNK_eq1s = nPointsSpace # one equation at the center which constrains u2, v2 at center
    total_BDNK_eq2s = nPointsSpace-1
    total_ODE_eqs = nPointsSpace-1

    total_eqs_block_1 = 0
    total_eqs_block_2 = total_eqs_block_1 + total_trivial_u1_eqs
    total_eqs_block_3 = total_eqs_block_2 + total_trivial_u2_eqs
    total_eqs_block_4 = total_eqs_block_3 + total_BDNK_eq1s
    total_eqs_block_5 = total_eqs_block_4 + total_BDNK_eq2s

    # center
    j = 1;
    linear_system_RHS[total_eqs_block_2 + j] = -(u2_np1[j] - u2_n[j] - k * (v2_np1[j] + v2_n[j]) / 2.0); # trivial time reduction crank nicholson equation.
    linear_system_RHS[total_eqs_block_3 + j] = -Center.Eq(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, h, k)
    
    for j in 2:nPointsSpace-1
        linear_system_RHS[total_eqs_block_1 + (j-1)] = -(u1_np1[j] - u1_n[j] - k * (v1_np1[j] + v1_n[j]) / 2.0); # trivial time reduction crank nicholson equation.
        linear_system_RHS[total_eqs_block_2 + j] = -(u2_np1[j] - u2_n[j] - k * (v2_np1[j] + v2_n[j]) / 2.0); # trivial time reduction crank nicholson equation.
        linear_system_RHS[total_eqs_block_3 + j] = -InteriorEqs.Eq1(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h, k)
        linear_system_RHS[total_eqs_block_4 + (j-1)] = -InteriorEqs.Eq2(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h, k)
        linear_system_RHS[total_eqs_block_5 + (j-1)] = -InteriorEqs.Eq3(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h)
    end

    # surface (backwards differences)
    j = nPointsSpace
    linear_system_RHS[total_eqs_block_1 + (j-1)] = -(u1_np1[j] - u1_n[j] - k * (v1_np1[j] + v1_n[j]) / 2.0); # trivial time reduction crank nicholson equation.
    linear_system_RHS[total_eqs_block_2 + j] = -(u2_np1[j] - u2_n[j] - k * (v2_np1[j] + v2_n[j]) / 2.0); # trivial time reduction crank nicholson equation.
    linear_system_RHS[total_eqs_block_3 + j] = -SurfaceEqs.Eq1(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h, k)
    linear_system_RHS[total_eqs_block_4 + (j-1)] = -SurfaceEqs.Eq2(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h, k)
    linear_system_RHS[total_eqs_block_5 + (j-1)] = -SurfaceEqs.Eq3(u1_np1, u2_np1, u3_np1, v1_np1, v2_np1, u1_n, u2_n, u3_n, v1_n, v2_n, CNSystem, j, h)
end

module InteriorEqs

Eq1(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (CNSystem.A1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.A2[j]*((-u1n[j-1] + u1n[j+1])/(2.0*h) + (-u1np1[j-1] + u1np1[j+1])/(2.0*h)))/2. + (CNSystem.A3[j]*((u1n[j-1] - 2*u1n[j] + u1n[j+1])/h^2 + (u1np1[j-1] - 2*u1np1[j] + u1np1[j+1])/h^2))/2. + (CNSystem.A7[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.A8[j]*((-u2n[j-1] + u2n[j+1])/(2.0*h) + (-u2np1[j-1] + u2np1[j+1])/(2.0*h)))/2. + (CNSystem.A9[j]*((u2n[j-1] - 2*u2n[j] + u2n[j+1])/h^2 + (u2np1[j-1] - 2*u2np1[j] + u2np1[j+1])/h^2))/2. + (CNSystem.A13[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.A14[j]*((-u3n[j-1] + u3n[j+1])/(2.0*h) + (-u3np1[j-1] + u3np1[j+1])/(2.0*h)))/2. + (CNSystem.A6[j]*(-v1n[j] + v1np1[j]))/k + (CNSystem.A4[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.A5[j]*((-v1n[j-1] + v1n[j+1])/(2.0*h) + (-v1np1[j-1] + v1np1[j+1])/(2.0*h)))/2. + (CNSystem.A12[j]*(-v2n[j] + v2np1[j]))/k + (CNSystem.A10[j]*(v2n[j] + v2np1[j]))/2. + (CNSystem.A11[j]*((-v2n[j-1] + v2n[j+1])/(2.0*h) + (-v2np1[j-1] + v2np1[j+1])/(2.0*h)))/2.


Eq2(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (CNSystem.B1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.B2[j]*((-u1n[j-1] + u1n[j+1])/(2.0*h) + (-u1np1[j-1] + u1np1[j+1])/(2.0*h)))/2. + (CNSystem.B3[j]*((u1n[j-1] - 2*u1n[j] + u1n[j+1])/h^2 + (u1np1[j-1] - 2*u1np1[j] + u1np1[j+1])/h^2))/2. + (CNSystem.B7[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.B8[j]*((-u2n[j-1] + u2n[j+1])/(2.0*h) + (-u2np1[j-1] + u2np1[j+1])/(2.0*h)))/2. + (CNSystem.B9[j]*((u2n[j-1] - 2*u2n[j] + u2n[j+1])/h^2 + (u2np1[j-1] - 2*u2np1[j] + u2np1[j+1])/h^2))/2. + (CNSystem.B13[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.B14[j]*((-u3n[j-1] + u3n[j+1])/(2.0*h) + (-u3np1[j-1] + u3np1[j+1])/(2.0*h)))/2. + (CNSystem.B6[j]*(-v1n[j] + v1np1[j]))/k + (CNSystem.B4[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.B5[j]*((-v1n[j-1] + v1n[j+1])/(2.0*h) + (-v1np1[j-1] + v1np1[j+1])/(2.0*h)))/2. + (CNSystem.B12[j]*(-v2n[j] + v2np1[j]))/k + (CNSystem.B10[j]*(v2n[j] + v2np1[j]))/2. + (CNSystem.B11[j]*((-v2n[j-1] + v2n[j+1])/(2.0*h) + (-v2np1[j-1] + v2np1[j+1])/(2.0*h)))/2.


Eq3(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64)::Float64 = (CNSystem.C1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.C2[j]*((-u1n[j-1] + u1n[j+1])/(2.0*h) + (-u1np1[j-1] + u1np1[j+1])/(2.0*h)))/2. + (CNSystem.C4[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.C5[j]*((-u2n[j-1] + u2n[j+1])/(2.0*h) + (-u2np1[j-1] + u2np1[j+1])/(2.0*h)))/2. + (CNSystem.C7[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.C8[j]*((-u3n[j-1] + u3n[j+1])/(2.0*h) + (-u3np1[j-1] + u3np1[j+1])/(2.0*h)))/2. + (CNSystem.C3[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.C6[j]*(v2n[j] + v2np1[j]))/2.

end

module SurfaceEqs

Eq1(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (CNSystem.A1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.A3[j]*((-u1n[j-3] + 4*u1n[j-2] - 5*u1n[j-1] + 2*u1n[j])/h^2 + (-u1np1[j-3] + 4*u1np1[j-2] - 5*u1np1[j-1] + 2*u1np1[j])/h^2))/2. + (CNSystem.A2[j]*((u1n[j-2] - 4*u1n[j-1] + 3*u1n[j])/(2.0*h) + (u1np1[j-2] - 4*u1np1[j-1] + 3*u1np1[j])/(2.0*h)))/2. + (CNSystem.A7[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.A9[j]*((-u2n[j-3] + 4*u2n[j-2] - 5*u2n[j-1] + 2*u2n[j])/h^2 + (-u2np1[j-3] + 4*u2np1[j-2] - 5*u2np1[j-1] + 2*u2np1[j])/h^2))/2. + (CNSystem.A8[j]*((u2n[j-2] - 4*u2n[j-1] + 3*u2n[j])/(2.0*h) + (u2np1[j-2] - 4*u2np1[j-1] + 3*u2np1[j])/(2.0*h)))/2. + (CNSystem.A13[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.A14[j]*((u3n[j-2] - 4*u3n[j-1] + 3*u3n[j])/(2.0*h) + (u3np1[j-2] - 4*u3np1[j-1] + 3*u3np1[j])/(2.0*h)))/2. + (CNSystem.A6[j]*(-v1n[j] + v1np1[j]))/k + (CNSystem.A4[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.A5[j]*((v1n[j-2] - 4*v1n[j-1] + 3*v1n[j])/(2.0*h) + (v1np1[j-2] - 4*v1np1[j-1] + 3*v1np1[j])/(2.0*h)))/2. + (CNSystem.A12[j]*(-v2n[j] + v2np1[j]))/k + (CNSystem.A10[j]*(v2n[j] + v2np1[j]))/2. + (CNSystem.A11[j]*((v2n[j-2] - 4*v2n[j-1] + 3*v2n[j])/(2.0*h) + (v2np1[j-2] - 4*v2np1[j-1] + 3*v2np1[j])/(2.0*h)))/2.


Eq2(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (CNSystem.B1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.B3[j]*((-u1n[j-3] + 4*u1n[j-2] - 5*u1n[j-1] + 2*u1n[j])/h^2 + (-u1np1[j-3] + 4*u1np1[j-2] - 5*u1np1[j-1] + 2*u1np1[j])/h^2))/2. + (CNSystem.B2[j]*((u1n[j-2] - 4*u1n[j-1] + 3*u1n[j])/(2.0*h) + (u1np1[j-2] - 4*u1np1[j-1] + 3*u1np1[j])/(2.0*h)))/2. + (CNSystem.B7[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.B9[j]*((-u2n[j-3] + 4*u2n[j-2] - 5*u2n[j-1] + 2*u2n[j])/h^2 + (-u2np1[j-3] + 4*u2np1[j-2] - 5*u2np1[j-1] + 2*u2np1[j])/h^2))/2. + (CNSystem.B8[j]*((u2n[j-2] - 4*u2n[j-1] + 3*u2n[j])/(2.0*h) + (u2np1[j-2] - 4*u2np1[j-1] + 3*u2np1[j])/(2.0*h)))/2. + (CNSystem.B13[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.B14[j]*((u3n[j-2] - 4*u3n[j-1] + 3*u3n[j])/(2.0*h) + (u3np1[j-2] - 4*u3np1[j-1] + 3*u3np1[j])/(2.0*h)))/2. + (CNSystem.B6[j]*(-v1n[j] + v1np1[j]))/k + (CNSystem.B4[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.B5[j]*((v1n[j-2] - 4*v1n[j-1] + 3*v1n[j])/(2.0*h) + (v1np1[j-2] - 4*v1np1[j-1] + 3*v1np1[j])/(2.0*h)))/2. + (CNSystem.B12[j]*(-v2n[j] + v2np1[j]))/k + (CNSystem.B10[j]*(v2n[j] + v2np1[j]))/2. + (CNSystem.B11[j]*((v2n[j-2] - 4*v2n[j-1] + 3*v2n[j])/(2.0*h) + (v2np1[j-2] - 4*v2np1[j-1] + 3*v2np1[j])/(2.0*h)))/2.

Eq3(u1np1::Vector{Float64}, u2np1::Vector{Float64}, u3np1::Vector{Float64}, v1np1::Vector{Float64}, v2np1::Vector{Float64}, u1n::Vector{Float64}, u2n::Vector{Float64}, u3n::Vector{Float64}, v1n::Vector{Float64}, v2n::Vector{Float64}, CNSystem, j::Int64, h::Float64)::Float64 = (CNSystem.C1[j]*(u1n[j] + u1np1[j]))/2. + (CNSystem.C2[j]*((u1n[j-2] - 4*u1n[j-1] + 3*u1n[j])/(2.0*h) + (u1np1[j-2] - 4*u1np1[j-1] + 3*u1np1[j])/(2.0*h)))/2. + (CNSystem.C4[j]*(u2n[j] + u2np1[j]))/2. + (CNSystem.C5[j]*((u2n[j-2] - 4*u2n[j-1] + 3*u2n[j])/(2.0*h) + (u2np1[j-2] - 4*u2np1[j-1] + 3*u2np1[j])/(2.0*h)))/2. + (CNSystem.C7[j]*(u3n[j] + u3np1[j]))/2. + (CNSystem.C8[j]*((u3n[j-2] - 4*u3n[j-1] + 3*u3n[j])/(2.0*h) + (u3np1[j-2] - 4*u3np1[j-1] + 3*u3np1[j])/(2.0*h)))/2. + (CNSystem.C3[j]*(v1n[j] + v1np1[j]))/2. + (CNSystem.C6[j]*(v2n[j] + v2np1[j]))/2.

end
end

module Jacobian
using ..Center

# set up jacobian. The rows correspond to the equations in the same order as the RHS above, while the columns correspond to the variables in the order u1_np1[2:N], u2_np1[2:N], v1_np1[2:N], v2_np1[2:N].
function compute_jacobian!(jacobian::AbstractArray{Float64}, CNSystem, h::Float64, k::Float64, nPointsSpace::Int64)
    # awkward indexing since we solve for (u1, v1, u3, w3) everywhere apart from center but solve for (w1, u2, w2, v2) everywhere including center
    total_trivial_u1_eqs = nPointsSpace-1
    total_trivial_u2_eqs = nPointsSpace
    total_BDNK_eq1s = nPointsSpace # one equation at the center which constrains u2, v2 at center
    total_BDNK_eq2s = nPointsSpace-1
    total_ODE_eqs = nPointsSpace-1

    total_eqs_block_1 = 0
    total_eqs_block_2 = total_eqs_block_1 + total_trivial_u1_eqs
    total_eqs_block_3 = total_eqs_block_2 + total_trivial_u2_eqs
    total_eqs_block_4 = total_eqs_block_3 + total_BDNK_eq1s
    total_eqs_block_5 = total_eqs_block_4 + total_BDNK_eq2s

    total_u1_vars = nPointsSpace-1
    total_u2_vars = nPointsSpace
    total_u3_vars = nPointsSpace-1
    total_u_vars = total_u1_vars + total_u2_vars + total_u3_vars

    total_v1_vars = nPointsSpace-1
    total_v2_vars = nPointsSpace
    total_v_vars = total_v1_vars + total_v2_vars

    ## center
    j = 1;
    # first equation: trivial time reduction for u2
    jacobian[total_eqs_block_2 + j, total_u1_vars + j] = 1.0; # ∂(Eq) / ∂(u2_np1[1+n, 2])
    jacobian[total_eqs_block_2 + j, total_u_vars + total_v1_vars + j] = -k / 2.0; # ∂(Eq) / ∂(v2_np1[1+n, 2])

    # fourth equation: equation for v2 at center
    jacobian[total_eqs_block_3 + j, 1] = Center.pCenterEq2pU1_2(CNSystem, h, k); # ∂(Eq) / ∂(u1_np1[1+n, 2])
    jacobian[total_eqs_block_3 + j, 2] = Center.pCenterEq2pU1_3(CNSystem, h, k); # ∂(Eq) / ∂(u1_np1[1+n, 3])

    jacobian[total_eqs_block_3 + j, total_u1_vars + 1] = Center.pCenterEq2pU2_1(CNSystem, h, k); # ∂(Eq) / ∂(u2_np1[1+n, 1])
    jacobian[total_eqs_block_3 + j, total_u1_vars + 2] = Center.pCenterEq2pU2_2(CNSystem, h, k); # ∂(Eq) / ∂(u2_np1[1+n, 2])
    jacobian[total_eqs_block_3 + j, total_u1_vars + 3] = Center.pCenterEq2pU2_3(CNSystem, h, k); # ∂(Eq) / ∂(u2_np1[1+n, 3])
    jacobian[total_eqs_block_3 + j, total_u1_vars + 4] = Center.pCenterEq2pU2_4(CNSystem, h, k); # ∂(Eq) / ∂(u2_np1[1+n, 4])

    jacobian[total_eqs_block_3 + j, total_u_vars + 1] = Center.pCenterEq2pV1_2(CNSystem, h, k); # ∂(Eq) / ∂(v1_np1[1+n, 2])
    jacobian[total_eqs_block_3 + j, total_u_vars + 2] = Center.pCenterEq2pV1_3(CNSystem, h, k); # ∂(Eq) / ∂(v1_np1[1+n, 3])

    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + 1] = Center.pCenterEq2pV2_1(CNSystem, h, k); # ∂(Eq) / ∂(v2_np1[1+n, 1])

    j = 2
    # first equation: trivial time reduction for u1
    jacobian[total_eqs_block_1 + (j-1), (j-1)] = 1.0; # ∂(Eq) / ∂(u1_np1[1+n, j])
    jacobian[total_eqs_block_1 + (j-1), total_u_vars + (j-1)] = -k / 2.0; # ∂(Eq) / ∂(v1_np1[1+n, j])

    # second equation: trivial time reduction for u2
    jacobian[total_eqs_block_2 + j, total_u1_vars + j] = 1.0; # ∂(Eq) / ∂(u2_np1[1+n, j])
    jacobian[total_eqs_block_2 + j, total_u_vars + total_v1_vars + j] = -k / 2.0; # ∂(Eq) / ∂(v2_np1[1+n, j])

    # sixth equation: Interior Eq1
    jacobian[total_eqs_block_3 + j, (j-1)] = Interior.pEq1pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
    jacobian[total_eqs_block_3 + j, (j-1) + 1] = Interior.pEq1pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u1_vars + j - 1] = Interior.pEq1pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u1_vars + j] = Interior.pEq1pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
    jacobian[total_eqs_block_3 + j, total_u1_vars + j + 1] = Interior.pEq1pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq1pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq1pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1)] = Interior.pEq1pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) + 1] = Interior.pEq1pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j - 1] = Interior.pEq1pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j] = Interior.pEq1pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j + 1] = Interior.pEq1pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])

    # seventh equation: Interior Eq2
    jacobian[total_eqs_block_4 + (j-1), (j-1)] = Interior.pEq2pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
    jacobian[total_eqs_block_4 + (j-1), (j-1) + 1] = Interior.pEq2pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j - 1] = Interior.pEq2pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j] = Interior.pEq2pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j + 1] = Interior.pEq2pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq2pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq2pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1)] = Interior.pEq2pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) + 1] = Interior.pEq2pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j - 1] = Interior.pEq2pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j] = Interior.pEq2pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j + 1] = Interior.pEq2pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])

    # eigth equation: Interior Eq2
    jacobian[total_eqs_block_5 + (j-1), (j-1)] = Interior.pEq3pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
    jacobian[total_eqs_block_5 + (j-1), (j-1) + 1] = Interior.pEq3pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j - 1] = Interior.pEq3pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j] = Interior.pEq3pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j + 1] = Interior.pEq3pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq3pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq3pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1)] = Interior.pEq3pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) + 1] = Interior.pEq3pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j - 1] = Interior.pEq3pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j] = Interior.pEq3pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j + 1] = Interior.pEq3pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])

    for j in 3:nPointsSpace-1
        # first equation: trivial time reduction for u1
        jacobian[total_eqs_block_1 + (j-1), (j-1)] = 1.0; # ∂(Eq) / ∂(u1_np1[1+n, j])
        jacobian[total_eqs_block_1 + (j-1), total_u_vars + (j-1)] = -k / 2.0; # ∂(Eq) / ∂(v1_np1[1+n, j])

        # second equation: trivial time reduction for u2
        jacobian[total_eqs_block_2 + j, total_u1_vars + j] = 1.0; # ∂(Eq) / ∂(u2_np1[1+n, j])
        jacobian[total_eqs_block_2 + j, total_u_vars + total_v1_vars + j] = -k / 2.0; # ∂(Eq) / ∂(v2_np1[1+n, j])

        # sixth equation: Interior Eq1
        jacobian[total_eqs_block_3 + j, (j-1) - 1] = Interior.pEq1pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
        jacobian[total_eqs_block_3 + j, (j-1)] = Interior.pEq1pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
        jacobian[total_eqs_block_3 + j, (j-1) + 1] = Interior.pEq1pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

        jacobian[total_eqs_block_3 + j, total_u1_vars + j - 1] = Interior.pEq1pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
        jacobian[total_eqs_block_3 + j, total_u1_vars + j] = Interior.pEq1pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
        jacobian[total_eqs_block_3 + j, total_u1_vars + j + 1] = Interior.pEq1pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

        jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) - 1] = Interior.pEq1pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
        jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq1pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
        jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq1pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

        jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) - 1] = Interior.pEq1pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
        jacobian[total_eqs_block_3 + j, total_u_vars + (j-1)] = Interior.pEq1pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
        jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) + 1] = Interior.pEq1pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

        jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j - 1] = Interior.pEq1pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
        jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j] = Interior.pEq1pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
        jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j + 1] = Interior.pEq1pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])

        # seventh equation: Interior Eq2
        jacobian[total_eqs_block_4 + (j-1), (j-1) - 1] = Interior.pEq2pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
        jacobian[total_eqs_block_4 + (j-1), (j-1)] = Interior.pEq2pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
        jacobian[total_eqs_block_4 + (j-1), (j-1) + 1] = Interior.pEq2pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j - 1] = Interior.pEq2pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j] = Interior.pEq2pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j + 1] = Interior.pEq2pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 1] = Interior.pEq2pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq2pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
        jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq2pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

        jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) - 1] = Interior.pEq2pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
        jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1)] = Interior.pEq2pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
        jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) + 1] = Interior.pEq2pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

        jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j - 1] = Interior.pEq2pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
        jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j] = Interior.pEq2pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
        jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j + 1] = Interior.pEq2pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])

        # eigth equation: Interior Eq2
        jacobian[total_eqs_block_5 + (j-1), (j-1) - 1] = Interior.pEq3pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
        jacobian[total_eqs_block_5 + (j-1), (j-1)] = Interior.pEq3pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])
        jacobian[total_eqs_block_5 + (j-1), (j-1) + 1] = Interior.pEq3pU1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j - 1] = Interior.pEq3pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j] = Interior.pEq3pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])
        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j + 1] = Interior.pEq3pU2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 1] = Interior.pEq3pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Interior.pEq3pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])
        jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) + 1] = Interior.pEq3pU3_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j+1])

        jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) - 1] = Interior.pEq3pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
        jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1)] = Interior.pEq3pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])
        jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) + 1] = Interior.pEq3pV1_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

        jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j - 1] = Interior.pEq3pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
        jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j] = Interior.pEq3pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
        jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j + 1] = Interior.pEq3pV2_jp1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j+1])
    end

    # surface point j = N
    j = nPointsSpace
    # first equation: trivial time reduction for u1
    jacobian[total_eqs_block_1 + (j-1), (j-1)] = 1.0; # ∂(Eq) / ∂(u1_np1[1+n, j])
    jacobian[total_eqs_block_1 + (j-1), total_u_vars + (j-1)] = -k / 2.0; # ∂(Eq) / ∂(v1_np1[1+n, j])

    # second equation: trivial time reduction for u2
    jacobian[total_eqs_block_2 + j, total_u1_vars + j] = 1.0; # ∂(Eq) / ∂(u2_np1[1+n, j])
    jacobian[total_eqs_block_2 + j, total_u_vars + total_v1_vars + j] = -k / 2.0; # ∂(Eq) / ∂(v2_np1[1+n, j])

    # sixth equation: Interior Eq1
    jacobian[total_eqs_block_3 + j, (j-1) - 3] = Surface.pEq1pU1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-3])
    jacobian[total_eqs_block_3 + j, (j-1) - 2] = Surface.pEq1pU1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-2])
    jacobian[total_eqs_block_3 + j, (j-1) - 1] = Surface.pEq1pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, (j-1)] = Surface.pEq1pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u1_vars + j - 3] = Surface.pEq1pU2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-3])
    jacobian[total_eqs_block_3 + j, total_u1_vars + j - 2] = Surface.pEq1pU2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-2])
    jacobian[total_eqs_block_3 + j, total_u1_vars + j - 1] = Surface.pEq1pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u1_vars + j] = Surface.pEq1pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) - 3] = Surface.pEq1pU3_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-3])
    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) - 2] = Surface.pEq1pU3_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-2])
    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1) - 1] = Surface.pEq1pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u1_vars + total_u2_vars + (j-1)] = Surface.pEq1pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])

    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) - 3] = Surface.pEq1pV1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-3])
    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) - 2] = Surface.pEq1pV1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-2])
    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1) - 1] = Surface.pEq1pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u_vars + (j-1)] = Surface.pEq1pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j+1])

    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j - 3] = Surface.pEq1pV2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-3])
    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j - 2] = Surface.pEq1pV2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-2])
    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j - 1] = Surface.pEq1pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_3 + j, total_u_vars + total_v1_vars + j] = Surface.pEq1pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])

    # seventh equation: Surface Eq2
    jacobian[total_eqs_block_4 + (j-1), (j-1) - 3] = Surface.pEq2pU1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-3])
    jacobian[total_eqs_block_4 + (j-1), (j-1) - 2] = Surface.pEq2pU1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-2])
    jacobian[total_eqs_block_4 + (j-1), (j-1) - 1] = Surface.pEq2pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), (j-1)] = Surface.pEq2pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])

    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j - 3] = Surface.pEq2pU2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-3])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j - 2] = Surface.pEq2pU2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-2])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j - 1] = Surface.pEq2pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + j] = Surface.pEq2pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])

    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 3] = Surface.pEq2pU3_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-3])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 2] = Surface.pEq2pU3_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-2])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 1] = Surface.pEq2pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Surface.pEq2pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])

    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) - 3] = Surface.pEq2pV1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-3])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) - 2] = Surface.pEq2pV1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-2])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1) - 1] = Surface.pEq2pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + (j-1)] = Surface.pEq2pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])

    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j - 3] = Surface.pEq2pV2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-3])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j - 2] = Surface.pEq2pV2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-2])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j - 1] = Surface.pEq2pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_4 + (j-1), total_u_vars + total_v1_vars + j] = Surface.pEq2pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])

    # eigth equation: Surface Eq2
    jacobian[total_eqs_block_5 + (j-1), (j-1) - 3] = Surface.pEq3pU1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-3])
    jacobian[total_eqs_block_5 + (j-1), (j-1) - 2] = Surface.pEq3pU1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-2])
    jacobian[total_eqs_block_5 + (j-1), (j-1) - 1] = Surface.pEq3pU1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), (j-1)] = Surface.pEq3pU1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u1_np1[1+n, j])

    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j - 3] = Surface.pEq3pU2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-3])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j - 2] = Surface.pEq3pU2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-2])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j - 1] = Surface.pEq3pU2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + j] = Surface.pEq3pU2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u2_np1[1+n, j])

    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 3] = Surface.pEq3pU3_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-3])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 2] = Surface.pEq3pU3_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-2])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1) - 1] = Surface.pEq3pU3_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u1_vars + total_u2_vars + (j-1)] = Surface.pEq3pU3_j(CNSystem, j, h, k); # ∂(Eq) / ∂(u3_np1[1+n, j])

    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) - 3] = Surface.pEq3pV1_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-3])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) - 2] = Surface.pEq3pV1_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-2])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1) - 1] = Surface.pEq3pV1_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + (j-1)] = Surface.pEq3pV1_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v1_np1[1+n, j])

    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j - 3] = Surface.pEq3pV2_jm3(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-3])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j - 2] = Surface.pEq3pV2_jm2(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-2])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j - 1] = Surface.pEq3pV2_jm1(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j-1])
    jacobian[total_eqs_block_5 + (j-1), total_u_vars + total_v1_vars + j] = Surface.pEq3pV2_j(CNSystem, j, h, k); # ∂(Eq) / ∂(v2_np1[1+n, j])
end
module Interior
# interior Eq 1
pEq1pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A3[j]/(2.0*h^2) - CNSystem.A2[j]/(4.0*h)
pEq1pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A1[j]/2. - CNSystem.A3[j]/h^2
pEq1pU1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A3[j]/(2.0*h^2) + CNSystem.A2[j]/(4.0*h)

pEq1pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A9[j]/(2.0*h^2) - CNSystem.A8[j]/(4.0*h)
pEq1pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A7[j]/2. - CNSystem.A9[j]/h^2
pEq1pU2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A9[j]/(2.0*h^2) + CNSystem.A8[j]/(4.0*h)

pEq1pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.A14[j]/h
pEq1pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A13[j]/2.
pEq1pU3_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A14[j]/(4.0*h)

pEq1pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.A5[j]/h
pEq1pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A4[j]/2. + CNSystem.A6[j]/k
pEq1pV1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A5[j]/(4.0*h)

pEq1pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.A11[j]/h
pEq1pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A10[j]/2. + CNSystem.A12[j]/k
pEq1pV2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A11[j]/(4.0*h)

# interior Eq 2
pEq2pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B3[j]/(2.0*h^2) - CNSystem.B2[j]/(4.0*h)
pEq2pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B1[j]/2. - CNSystem.B3[j]/h^2
pEq2pU1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B3[j]/(2.0*h^2) + CNSystem.B2[j]/(4.0*h)

pEq2pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B9[j]/(2.0*h^2) - CNSystem.B8[j]/(4.0*h)
pEq2pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B7[j]/2. - CNSystem.B9[j]/h^2
pEq2pU2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B9[j]/(2.0*h^2) + CNSystem.B8[j]/(4.0*h)

pEq2pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.B14[j]/h
pEq2pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B13[j]/2.
pEq2pU3_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B14[j]/(4.0*h)

pEq2pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.B5[j]/h
pEq2pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B4[j]/2. + CNSystem.B6[j]/k
pEq2pV1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B5[j]/(4.0*h)

pEq2pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.B11[j]/h
pEq2pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B10[j]/2. + CNSystem.B12[j]/k
pEq2pV2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B11[j]/(4.0*h)

# interior Eq 3
pEq3pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.C2[j]/h
pEq3pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C1[j]/2.
pEq3pU1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C2[j]/(4.0*h)

pEq3pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.C5[j]/h
pEq3pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C4[j]/2.
pEq3pU2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C5[j]/(4.0*h)

pEq3pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.25*CNSystem.C8[j]/h
pEq3pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C7[j]/2.
pEq3pU3_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C8[j]/(4.0*h)

pEq3pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C3[j]/2.
pEq3pV1_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0

pEq3pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C6[j]/2.
pEq3pV2_jp1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
end

module Surface
# surface Eq 1
pEq1pU1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.5*CNSystem.A3[j]/h^2
pEq1pU1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (2*CNSystem.A3[j])/h^2 + CNSystem.A2[j]/(4.0*h)
pEq1pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (-5*CNSystem.A3[j])/(2.0*h^2) - CNSystem.A2[j]/h
pEq1pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A1[j]/2. + CNSystem.A3[j]/h^2 + (3*CNSystem.A2[j])/(4.0*h)

pEq1pU2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.5*CNSystem.A9[j]/h^2
pEq1pU2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (2*CNSystem.A9[j])/h^2 + CNSystem.A8[j]/(4.0*h)
pEq1pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (-5*CNSystem.A9[j])/(2.0*h^2) - CNSystem.A8[j]/h
pEq1pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A7[j]/2. + CNSystem.A9[j]/h^2 + (3*CNSystem.A8[j])/(4.0*h)

pEq1pU3_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq1pU3_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A14[j]/(4.0*h)
pEq1pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.A14[j]/h)
pEq1pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A13[j]/2. + (3*CNSystem.A14[j])/(4.0*h)

pEq1pV1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq1pV1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A5[j]/(4.0*h)
pEq1pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.A5[j]/h)
pEq1pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A4[j]/2. + (3*CNSystem.A5[j])/(4.0*h) + CNSystem.A6[j]/k

pEq1pV2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq1pV2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A11[j]/(4.0*h)
pEq1pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.A11[j]/h)
pEq1pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.A10[j]/2. + (3*CNSystem.A11[j])/(4.0*h) + CNSystem.A12[j]/k

# surface Eq 2
pEq2pU1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.5*CNSystem.B3[j]/h^2
pEq2pU1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (2*CNSystem.B3[j])/h^2 + CNSystem.B2[j]/(4.0*h)
pEq2pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (-5*CNSystem.B3[j])/(2.0*h^2) - CNSystem.B2[j]/h
pEq2pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B1[j]/2. + CNSystem.B3[j]/h^2 + (3*CNSystem.B2[j])/(4.0*h)

pEq2pU2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -0.5*CNSystem.B9[j]/h^2
pEq2pU2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (2*CNSystem.B9[j])/h^2 + CNSystem.B8[j]/(4.0*h)
pEq2pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = (-5*CNSystem.B9[j])/(2.0*h^2) - CNSystem.B8[j]/h
pEq2pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B7[j]/2. + CNSystem.B9[j]/h^2 + (3*CNSystem.B8[j])/(4.0*h)

pEq2pU3_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq2pU3_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B14[j]/(4.0*h)
pEq2pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.B14[j]/h)
pEq2pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B13[j]/2. + (3*CNSystem.B14[j])/(4.0*h)

pEq2pV1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq2pV1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B5[j]/(4.0*h)
pEq2pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.B5[j]/h)
pEq2pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B4[j]/2. + (3*CNSystem.B5[j])/(4.0*h) + CNSystem.B6[j]/k

pEq2pV2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq2pV2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B11[j]/(4.0*h)
pEq2pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.B11[j]/h)
pEq2pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.B10[j]/2. + (3*CNSystem.B11[j])/(4.0*h) + CNSystem.B12[j]/k

# surface Eq 3
pEq3pU1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pU1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C2[j]/(4.0*h)
pEq3pU1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.C2[j]/h)
pEq3pU1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C1[j]/2. + (3*CNSystem.C2[j])/(4.0*h)

pEq3pU2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pU2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C5[j]/(4.0*h)
pEq3pU2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.C5[j]/h)
pEq3pU2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C4[j]/2. + (3*CNSystem.C5[j])/(4.0*h)

pEq3pU3_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pU3_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C8[j]/(4.0*h)
pEq3pU3_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = -(CNSystem.C8[j]/h)
pEq3pU3_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C7[j]/2. + (3*CNSystem.C8[j])/(4.0*h)

pEq3pV1_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV1_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV1_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV1_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C3[j]/2.

pEq3pV2_jm3(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV2_jm2(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV2_jm1(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = 0
pEq3pV2_j(CNSystem, j::Int64, h::Float64, k::Float64)::Float64 = CNSystem.C6[j]/2.
end

end

end
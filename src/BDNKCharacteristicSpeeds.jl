#=

    Module comprising of functions for computing the characteristic speeds of the BDNK system

=#

module BDNKCharacteristicSpeeds

c_plus_sq(Λ0::Vector{Float64}, Λ1::Vector{Float64}, Λ2::Vector{Float64})::Vector{Float64} = @. (Λ1 + sqrt(Λ1^2 - Λ0)) / (Λ2);
c_minus_sq(Λ0::Vector{Float64}, Λ1::Vector{Float64}, Λ2::Vector{Float64})::Vector{Float64} = @. (Λ1 - sqrt(Λ1^2 - Λ0)) / (Λ2);

function compute(p::Float64, ε::Float64, cs::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)
    λ0 = Λ0(p, ε, cs, η, ζ, τε, τP, τQ, L)
    λ1 = Λ1(p, ε, cs, η, ζ, τε, τP, τQ, L)
    λ2 = Λ2(p, ε, cs, η, ζ, τε, τP, τQ, L)

    return λ0, λ1, λ2
end

Λ0(p::Float64, ε::Float64, cs::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (4*L^4*(3*ζ + 4*η)^4*(-1 + τP)*τQ^2*τε*(p + ε)^2*cs^4)/81.

Λ1(p::Float64, ε::Float64, cs::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (L^2*(3*ζ + 4*η)^2*(τε + τQ*(τP + τε))*(p + ε)*cs^2)/9.

Λ2(p::Float64, ε::Float64, cs::Float64, η::Float64, ζ::Float64, τε::Float64, τP::Float64, τQ::Float64, L::Float64)::Float64 = (2*L^2*(3*ζ + 4*η)^2*τQ*τε*(p + ε))/9.

end
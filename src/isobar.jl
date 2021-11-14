
function J_I(σ,pars::NamedTuple{(:m,:Γ)})
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1/(m^2 - σ + 1im*m*Γ*FF)
end

function J_II(σ,pars::NamedTuple{(:m,:Γ)})
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1/(m^2 - σ - 1im*m*Γ*FF)
end

pole_position(R::NamedTuple{(:m,:Γ),T} where T) = R.m^2-1im*R.m*R.Γ


abstract type AbstractLinesShape end

@with_kw struct BW <: AbstractLinesShape
    m::Float64
    Γ::Float64
end
J_I( σ, pars::BW) = 1.0 / (pars.m^2 - σ - 1im*pars.m*pars.Γ)
J_II(σ, pars::BW) = 1.0 / (pars.m^2 - σ + 1im*pars.m*pars.Γ)

@with_kw struct BW_Swave <: AbstractLinesShape
    m::Float64
    g::Float64
    ma::Float64
    mb::Float64
end
J_I( σ, pars::BW_Swave) = 1 / (pars.m^2 - σ - 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)
J_II(σ, pars::BW_Swave) = 1 / (pars.m^2 - σ + 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)

@with_kw struct BW_norm <: AbstractLinesShape
    m::Float64
    Γ::Float64
	c::Float64
end
J_I( σ, pars::BW_norm) = pars.c / (pars.m^2 - σ - 1im*pars.m*pars.Γ)
J_II(σ, pars::BW_norm) = pars.c / (pars.m^2 - σ + 1im*pars.m*pars.Γ)


function obj2nt(amp::AbstractLinesShape)
	type = typeof(amp)
	names = fieldnames(type)
	nt = NamedTuple{names}([getproperty(amp, n) for n in names])
	(; type=string(type), nt...)
end

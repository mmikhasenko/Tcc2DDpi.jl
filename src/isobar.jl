
function Jᴵ(σ,pars::NamedTuple{(:m,:Γ)})
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1/(m^2 - σ + 1im*m*Γ*FF)
end

function Jᴵᴵ(σ,pars::NamedTuple{(:m,:Γ)})
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
Jᴵ( σ, pars::BW) = σ == pars.m^2-1im*pars.m*pars.Γ ? 0.0 : 1.0 / (pars.m^2 - σ + 1im*pars.m*pars.Γ)
Jᴵᴵ(σ, pars::BW) = σ == pars.m^2+1im*pars.m*pars.Γ ? 0.0 : 1.0 / (pars.m^2 - σ - 1im*pars.m*pars.Γ)

pole_position(R::BW) = R.m^2-1im*R.m*R.Γ

@with_kw struct BW_Swave <: AbstractLinesShape
    m::Float64
    g::Float64
    ma::Float64
    mb::Float64
end
Jᴵ( σ, pars::BW_Swave) = 1 / (pars.m^2 - σ + 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)
Jᴵᴵ(σ, pars::BW_Swave) = 1 / (pars.m^2 - σ - 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)

@with_kw struct BW_norm <: AbstractLinesShape
    m::Float64
    Γ::Float64
	n::Float64
end
Jᴵ( σ, pars::BW_norm) = pars.n / (pars.m^2 - σ + 1im*pars.m*pars.Γ)
Jᴵᴵ(σ, pars::BW_norm) = pars.n / (pars.m^2 - σ - 1im*pars.m*pars.Γ)

@with_kw struct ZeroBW <: AbstractLinesShape
    m::Float64
    Γ::Float64
end
Jᴵ(σ::Number,pars::ZeroBW) = 0.0
Jᴵᴵ(σ::Number,pars::ZeroBW) = 0.0


function obj2nt(amp::AbstractLinesShape)
	type = typeof(amp)
	names = fieldnames(type)
	nt = NamedTuple{names}([getproperty(amp, n) for n in names])
	(; type=string(type), nt...)
end

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


struct BW_Swave
    m::Float64
    g::Float64
    ma::Float64
    mb::Float64
end
J_I( σ, pars::BW_Swave) = 1 / (pars.m^2 - σ - 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)
J_II(σ, pars::BW_Swave) = 1 / (pars.m^2 - σ + 1im*pars.g^2*sqrt(λ(σ,pars.ma^2,pars.mb^2))/σ)

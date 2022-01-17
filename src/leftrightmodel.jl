# πDD
struct ChannelWithIntegrationMethod{
        T <: AbstractxDD,
        H <:AbstractDalitzMapping} <: AbstractxDD
	channel::T
	mapdalitzmethod::H
end


# default: take everything from the channel
decay_matrix_element_squared(d::ChannelWithIntegrationMethod,s,σ3,σ2) =
	decay_matrix_element_squared(d.channel,s,σ3,σ2)
masses(set::ChannelWithIntegrationMethod) = masses(set.channel)
mapdalitzmethod(set::ChannelWithIntegrationMethod) = set.mapdalitzmethod


decay_matrix_element_squared(
    set::ChannelWithIntegrationMethod{T,HookSqrtDalitzMapping{3}} where T<:πDD,s,σ3,σ2) = 
        covertapply(πDD_𝔐²_nonana3,set.channel,s,σ3,σ2)
decay_matrix_element_squared(
    set::ChannelWithIntegrationMethod{T,HookSqrtDalitzMapping{2}} where T<:πDD,s,σ3,σ2) = 
        covertapply(πDD_𝔐²_nonana2,set.channel,s,σ3,σ2)


decay_matrix_element_squared(
    set::ChannelWithIntegrationMethod{T,HookSqrtDalitzMapping{3}} where T<:γDD,s,σ3,σ2) = 
    covertapply(γDD_𝔐²_nonana3,set.channel,s,σ3,σ2)
decay_matrix_element_squared(
    set::ChannelWithIntegrationMethod{T,HookSqrtDalitzMapping{2}} where T<:γDD,s,σ3,σ2) = 
    covertapply(γDD_𝔐²_nonana2,set.channel,s,σ3,σ2)
        


function projecttocos23(d, s, cos23)
    !(-1 < cos23 < 1) && return 0.0
    
    σ1_0, σ1_e = (d.ms[2]+d.ms[3])^2, (√s-d.ms[1])^2
        
    M²of1(σ1) = real(
        X2DDpi.decay_matrix_element_squared(d,s,
            X2DDpi.σ3of1(σ1, cos23, d.ms^2, s),
            X2DDpi.σ2of1(σ1, cos23, d.ms^2, s))
        ) *
        sqrt(λ(σ1,d.ms[2]^2,d.ms[3]^2)*λ(s,σ1,d.ms[1]^2)) / σ1
    return quadgk(M²of1, σ1_0, σ1_e)[1] / (2π*s)
end

function projectto1(d, s, σ1)

    σ1_0, σ1_e = (d.ms[2]+d.ms[3])^2, (√s-d.ms[1])^2
    !(σ1_0 < σ1 < σ1_e) && return 0.0
    # 
	σ3_0, σ3_e = X2DDpi.σ3of1_pm(σ1, d.ms^2, s)
    
    M²of3(σ3) = real(
        X2DDpi.decay_matrix_element_squared(d,s,
            σ3,sum(d.ms^2)+s-σ3-σ1))
    return quadgk(M²of3, σ3_0, σ3_e)[1] / (2π*s)
end
#

function projectto1(d, σ1, lineshape, smin, smax)
    σ1_0 = (d.ms[2]+d.ms[3])^2
    σ1 < σ1_0 && return 0.0
    # 
    M²of3(s, σ3) = real(
        X2DDpi.decay_matrix_element_squared(d,s,
            σ3,sum(d.ms^2)+s-σ3-σ1)) * abs2(lineshape(s))
    function M²of3(x) # x = [xs, xσ3]
        s = smin + x[1]*(smax-smin)
        # 
        σ1_e = (√s-d.ms[1])^2
        σ1 > σ1_e && return 0.0
        # 
        σ3_0, σ3_e = X2DDpi.σ3of1_pm(σ1, d.ms^2, s)
        σ3 = σ3_0 + x[2]*(σ3_e-σ3_0)
        M²of3(s, σ3)*(σ3_e-σ3_0)*(smax-smin) / (2π*s)
    end
    return cuhre((x,f)->(f[1]=M²of3(x)), 2, 1)[1][1]
end


function projectto3(d, s, σ3)

    σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    !(σ3_0 < σ3 < σ3_e) && return 0.0
    # 
	σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
    
    M²of2(σ2) = real(X2DDpi.decay_matrix_element_squared(d,s,σ3,σ2))
    return quadgk(M²of2, σ2_0, σ2_e)[1] / (2π*s)
end
# 
function projectto3(d, σ3, lineshape, smin, smax)
    σ3_0 = (d.ms[1]+d.ms[2])^2
    σ3 < σ3_0 && return 0.0
    # 
    M²of2(s, σ2) = real(X2DDpi.decay_matrix_element_squared(d,s,σ3,σ2)) * abs2(lineshape(s))
    function M²of2(x)
        s = smin + x[1]*(smax-smin)
        # 
        σ3_e = (√s-d.ms[3])^2
        σ3 > σ3_e && return 0.0
        # 
        σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
        σ2 = σ2_0 + x[2]*(σ2_e-σ2_0)
        M²of2(s, σ2)*(σ2_e-σ2_0)*(smax-smin) / (2π*s)
    end 
    return cuhre((x,f)->(f[1]=M²of2(x)), 2, 1)[1][1]
end

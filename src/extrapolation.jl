const RegressionResult = NamedTuple{(:α, :β, :SE_α, :SE_β)}
import Base: *
*(nt::RegressionResult, f::Number) = RegressionResult(collect(nt) .* f)

function Base.show(io::IO, nt::RegressionResult)
    print(io,
        "Interpolation to zero:",
        "\n\tintercept: ", round(nt.β, digits=2),
        "\n\tinterpolation error: ", round(nt.SE_β, digits=2),
        "\n\tslope: ", round(nt.α, digits=2),"\n")
end


"""
    linear_regression(xv, yv)

Perform linear regression on the given data points `xv` and `yv`.

# Arguments
- `xv::Vector`: A vector of x-values.
- `yv::Complex`: A vector of corresponding y-values.

# Returns
- `α::Complex`: The slope of the linear regression.
- `β::Complex`: The intercept of the linear regression.
- `SE_α_over_α::Complex`: The standard error of the slope as δα_re/|α| + j*δα_im/|α|.
- `SE_β_over_β::Complex`: The standard error of the intercept as δβ_re/|α| + j*δβ_im/|α|.

# Example
```julia
x = [1.0, 2.0, 3.0, 4.0, 5.0]
y = [2.1, 4.0, 6.1, 8.0, 10.2]
result = linear_regression(x, y)

# Note
The function computes the intercept and slope using the formulae derived from the method of least squares.
    It also computes the standard errors for the intercept and slope to provide a measure of the uncertainty associated with these estimates. 
"""
function linear_regression(xv, yv)
    Ex = mean(xv)
    Ey = mean(yv)
    Exy = mean(xv .* yv)
    Ex² = mean(xv .^ 2)
    Dx = Ex² - Ex^2

    # α is a slope, β is the intercept 
    α, β = [Exy - Ex * Ey, Ey * Ex² - Ex * Exy] ./ Dx

    # Compute residuals
    residuals = yv .- (β .+ α .* xv)

    # Compute residual sum of squares
    RSS_real = sum(real.(residuals) .^ 2)
    RSS_imag = sum(imag.(residuals) .^ 2)

    # Compute standard error of the regression
    n = length(xv)
    dof = n - 2
    SER_real = sqrt(RSS_real / dof)
    SER_imag = sqrt(RSS_imag / dof)

    # Compute standard error of the slope
    SE_β_real = SER_real / sqrt(sum((xv .- Ex) .^ 2))
    SE_β_imag = SER_imag / sqrt(sum((xv .- Ex) .^ 2))

    # Compute standard error of the intercept
    SE_α_real = SE_β_real * sqrt(1/n + Ex^2 / sum((xv .- Ex) .^ 2))
    SE_α_imag = SE_β_imag * sqrt(1/n + Ex^2 / sum((xv .- Ex) .^ 2))

    # SE complex
    SE_α = SE_α_real + 1im*SE_α_imag
    SE_β = SE_β_real + 1im*SE_β_imag

    return RegressionResult((; α, β, SE_α, SE_β))
end


function effective_range_scan(model, settings)
    @unpack ϵ0 = settings
    @unpack ϵf_grid = settings
    @unpack circular_sum_n = settings
    @unpack model_δm0 = settings
    # 
    method(r) = ComplexBranchPointExpansion(CircularSum(r, circular_sum_n))

    scan = map(ϵf_grid) do ϵf
        circle_r = ϵf * ϵ0
        effrangepars =
            effectiverangeexpansion(
                Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, model_δm0),
                Δe -> k3b(Eᵦˣ⁺ + Δe),
                method(circle_r))
        efe = (; tophysicsunits(effrangepars)..., effrangepars..., ϵf, circle_r)
        efe
    end

    return scan
end

function linear_regression_in_sqrt_ϵf(scan, quantity)
    ϵfv = getproperty.(scan, :ϵf)
    sqrt_ϵfv = sqrt.(ϵfv)
    # 
    calv = getproperty.(scan, quantity)
    regression = linear_regression(sqrt_ϵfv, calv)
    # 
    return regression
end

function effective_range_extrapolation(scan)
    r_regression = linear_regression_in_sqrt_ϵf(scan, :r)
    inva_regression = linear_regression_in_sqrt_ϵf(scan, :a⁻¹)
    # 
    (; r_regression, inva_regression)
end

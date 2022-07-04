# Analysis of Tcc⁺ → D⁰D⁰π⁺ amplitude

The code is used for several projects:

- [LHCb-ANA-2020-041 (internal twiki page)](https://twiki.cern.ch/twiki/bin/viewauth/LHCbPhysics/X2DDstar)
- [LHCb-PAPER-2021-031 (public page)](https://lhcbproject.web.cern.ch/Publications/LHCbProjectPublic/LHCb-PAPER-2021-032.html) published by [Nature Physics](https://www.nature.com/articles/s41567-022-01614-y)
- [LHCb-PAPER-2021-032 (public page)](https://lhcbproject.web.cern.ch/Publications/LHCbProjectPublic/LHCb-PAPER-2021-031.html) published by [Nature Communication](https://www.nature.com/articles/s41467-022-30206-w)
- [Arxiv:2203.04622](https://arxiv.org/abs/2203.04622) Accurate computation of the effective-range parameters

The code implements an analytic continuation of the amplitude to the complex plane and calculation of the effective range parameters.

The main types exported are:

```julia
struct πDD <: AbstractxDD  # both isospin channels
struct γDD <: AbstractxDD
```

with a specific `decay_matrix_element_squared(d::AbstractxDD,s,σ3,σ2)` for every type.

The decay amplitude integrated over the phase space is called via

```julia
ρ_thr(d::AbstractxDD, e)
```

The latter can be precalculated using the `interpolated` method,

```julia
interpolated(channel::AbstractxDD, cutoff::Real; estep=0.01)
```

at the value of `cutoff` the expression is merged to the two-body phase-space function. The default value of the discrete grid is 0.01 MeV.

The precalculated structure can be used to calculate the dispersion integral:

```julia
dispersive(d::interpolated, e)
```

returns a complex value, for which `Im(dispersive) = ρ_thr`.

## Analysis scripts

1. **Default analysis**

   - `model.jl` build the default model, cache the tables
   - `computeasymptotics.jl` save ρInf
   - `effectiverange.jl` value of the scattering parameters using the matching method
   - `plotinverseamplitude.jl` matching plot for the inverse amplitude
   - `plotwidthdependance.jl` profile for different values Γ
   - `DpDz_spectrum.jl` projection to D+D0
   - `Dpi_spectrum.jl` projection to Dpi
   - `polesearch.jl` default pole position
   - `polescang.jl` the width saturation curve
   - `plotcomplexplane.jl` plot the complex plane for the paper
   - `interference_matrix.jl` contribution of the decay-chain interference
   - `pimDpDp.jl` why π⁻D⁺D⁺ is negligible
   - `loopdiagrams.jl` plot loop diagrams for the supplemental material
   - `Dxcomplexmass.jl` Tcc parameters in the model with complex D* mass

2. **Angular dependence**

   - `dalitz_plots.jl` dalitz plot for different quantum numbers
   - `angular_distribution.jl` anglular distribution (sympy) along the D* band
   - `profile_peakandtail.jl` fraction of the tail in different models
   - `DpionL.jl` Sensiticity of the Dpi spectrum on quantum numbers

3. **Isospin dependence**

   - `pole_vs_g1g2.jl`

4. Plotting complex plane in with `plotlyjs` for [EP news](https://ep-news.web.cern.ch/content/lhcb-discovers-double-charm-tetraquark)

   - `complexplaneongrid.jl` Compute the values on the complex grid using the default model
   - `DxD_model_without_interference.jl` Compute the values on the complex grid using a simplified model
   - `plotcomplexplane_3d_plotlyjs.jl` Visualize the values

5. Effective range computation using the cauchy integrals

   - `effectiverangecauchy.jl`
   - `plotcompositeness.jl`

6. Additional plots for talks
   - `howwidthismade.jl` contributions to the width
   - `xxx.jl` DDpi momentum spectrum

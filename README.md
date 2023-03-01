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

   - [`S-000-model.jl`](scripts/S-000-model.jl) build the default model, cache the tables
   - [`S-001-computeasymptotics.jl`](scripts/S-001-computeasymptotics.jl) save ρInf
   - [`S-002-effectiverange.jl`](scripts/S-002-effectiverange.jl) value of the scattering parameters using the matching method
   - [`S-003-plotinverseamplitude.jl`](scripts/S-003-plotinverseamplitude.jl) matching plot for the inverse amplitude
   - [`S-004-plotwidthdependance.jl`](scripts/S-004-plotwidthdependance.jl) profile for different values Γ
   - [`S-005-DpDz_spectrum.jl`](scripts/S-005-DpDz_spectrum.jl) projection to D+D0
   - [`S-006-Dpi_spectrum.jl`](scripts/S-006-Dpi_spectrum.jl) projection to Dpi
   - [`S-007-polesearch.jl`](scripts/S-007-polesearch.jl) default pole position
   - [`S-008-polescang.jl`](scripts/S-008-polescang.jl) the width saturation curve
   - [`S-009-plotcomplexplane.jl`](scripts/S-009-plotcomplexplane.jl) plot the complex plane for the paper
   - [`S-010-interference_matrix.jl`](scripts/S-010-interference_matrix.jl) contribution of the decay-chain interference
   - [`S-011-pimDpDp.jl`](scripts/S-011-pimDpDp.jl) why π⁻D⁺D⁺ is negligible
   - [`S-012-loopdiagrams.jl`](scripts/S-012-loopdiagrams.jl) plot loop diagrams for the supplemental material
   - [`S-013-Dxcomplexmass.jl`](scripts/S-013-Dxcomplexmass.jl) Tcc parameters in the model with complex D* mass

2. **Angular dependence**

   - [`S-100-dalitz_plots.jl`](scripts/S-100-dalitz_plots.jl) dalitz plot for different quantum numbers
   - [`S-101-angular_distribution.jl`](scripts/S-101-angular_distribution.jl) anglular distribution (sympy) along the D* band
   - [`S-102-profile_peakandtail.jl`](scripts/S-102-profile_peakandtail.jl) fraction of the tail in different models
   - [`S-103-DpionL.jl`](scripts/S-103-DpionL.jl) Sensiticity of the Dpi spectrum on quantum numbers

3. **Isospin dependence**

   - [`S-200-pole_vs_g1g2.jl`](scripts/S-200-pole_vs_g1g2.jl)

4. **Complex plane** for [EP news](https://ep-news.web.cern.ch/content/lhcb-discovers-double-charm-tetraquark)

   - [`S-300-complexplaneongrid.jl`](scripts/S-300-complexplaneongrid.jl) Compute the values on the complex grid using the default model
   - [`S-301-DxD_model_without_interference.jl`](scripts/S-301-DxD_model_without_interference.jl) Compute the values on the complex grid using a simplified model
   - [`S-302-plotcomplexplane_3d_plotlyjs.jl`](scripts/S-302-plotcomplexplane_3d_plotlyjs.jl) Visualize the values
   - [`S-303-integration_domain.jl`](scripts/S-303-integration_domain.jl) Validate that no singularities happen in the integration domain.

5. **Scattering parameters using the cauchy integrals**

   - [`S-400-effectiverangecauchy.jl`](scripts/S-400-effectiverangecauchy.jl)
   - [`S-401-plotcompositeness.jl`](scripts/S-401-plotcompositeness.jl)
   - [`S-402-effective_range_with_fit.jl`](scripts/S-402-effective_range_with_fit.jl)

6. **Additional plots for talks**

   - [`S-900-howwidthismade.jl`](scripts/S-900-howwidthismade.jl) contributions to the width
   - [`S-901-xxx.jl`](scripts/S-901-xxx.jl) DDpi momentum spectrum. Femtoscopy
   - [`N-902-triangle_singularity.jl`](notebooks/N-902-triangle_singularity.jl) Position of the traingle singularity in the complex energy plane

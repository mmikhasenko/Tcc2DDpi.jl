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

## Mass and energy conventions

Phase space and dispersion use **absolute masses in GeV**. These routines
use the masses stored in each channel, so they also apply to channels with
masses different from the nominal Tcc system. This matches the mass-parametric
API in the sibling `../xDDPhaseSpace.jl` package.

| API | Units |
| --- | --- |
| `ρ_thr(channel, m)`, `ρ_tb(channel, m)` | Absolute mass `m` in GeV |
| `interpolated(channel, m_cutoff; mstep=1e-5)` | Cutoff and grid spacing in GeV |
| `dispersive(ichannel, m)` | Absolute real or complex mass in GeV |
| `denominator_I/II(amp, e, δm)` | Offsets from D⁰D*⁺ threshold in MeV |
| `pole_position(amp, δm)` | Input and returned pole offsets in MeV |
| `k3b(e)`, `kNR(e)` | Energy offset in MeV; momentum in GeV |
| Dalitz integration and projection arguments `s`, `σᵢ` | Squared masses in GeV² |

For ordinary numbers, `e2m(e) = mD⁰ + mDˣ⁺ + e/1000` and `m2e` is its
inverse. The amplitude converts at the boundary; its correction function
continues to receive energy in MeV. The phase-space/interpolation/dispersion
kernels do not use that reference threshold. For another mass system, construct
channels with its masses and call those kernels with absolute masses directly.
The Tcc amplitude and momentum helpers retain their Tcc conventions.
`WithThrE(e, reference_mass)` from xDDPhaseSpace can express a different
reference for conversion with `e2m`; convert it before mass-API calls. It is
not a general arithmetic or unit-tracking type.

```julia
using X2DDpi

channel = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰),
    BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
e_cutoff, e_step = 4.0, 0.003  # MeV, as in settings.json
ichannel = interpolated(channel, e2m(e_cutoff); mstep=e_step/1000)
e = -0.3592  # MeV
ρ = ρ_thr(ichannel, e2m(e))
Σ = dispersive(ichannel, e2m(e))
D = denominator_I(Amplitude(ichannel), e, e)
```

The interpolator returns zero below the three-body threshold, uses a linear
mass grid below `m_cutoff`, and matches to a rescaled two-body phase space above
it. Complex masses bypass interpolation and evaluate the underlying channel.
The once-subtracted dispersion integral integrates over squared mass, with
subtraction at `s=0`. Its upper-rim imaginary part approaches `ρ_thr`; real
arguments use a `1e-9im` GeV regulator, equivalent to the former `1e-6im` MeV.

### Migrating callers and caches

Version 0.2 changes the meaning of the phase-space and dispersion arguments.
Replace `ρ_thr(channel, e)` with `ρ_thr(channel, e2m(e))`, and do the same for
`ρ_tb` and `dispersive`. For interpolation, replace
`interpolated(channel, cutoff; estep=step)` with
`interpolated(channel, e2m(cutoff); mstep=step/1000)` when the old values are MeV.
The draft PR's intermediate `e_step` keyword is also replaced by `mstep` in GeV.
Settings use `e_cutoff` and `e_step`, both in MeV. Convert them explicitly when
constructing an interpolator. `mstep` matches xDDPhaseSpace's keyword and units.

New serialized interpolators store `m_cutoff` and `knots` in GeV; `coefs` and
`cutoffratio` retain their normalization. `interpolated(d2nt(entry))` also reads
historical Tcc caches with `cutoff` and knots in MeV, converting coordinates on
load. That legacy format always uses the historical D⁰D*⁺ reference threshold;
it must not be used to encode an arbitrary mass reference. Field order is
irrelevant, including nested masses and lineshapes. Historical photon channels
without explicit couplings receive the analysis defaults `(1.0, -3.77)`.

X2DDpi retains its own interpolator and dispersion implementation and the
historical-cache adapter. It uses xDDPhaseSpace's mass-parametric channel
kernels; the two packages' interpolator types and cache formats are not aliases.

### Validation

With the sibling xDDPhaseSpace checkout available at `../xDDPhaseSpace.jl`:

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project test/runtests.jl
```

The suite covers the mass/MeV boundary, channels with different masses,
serialization order and legacy caches, interpolation versus direct phase space,
and real/complex dispersion and amplitude continuation.

To refresh the published notebooks and documentation, see the [local site build instructions](website/README.md).

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

   - [`S-400-effectiverangecauchy.jl`](scripts/S-400-effectiverangecauchy.jl) compute the eff. range parameters for 5 models
   - [`S-401-plotcompositeness.jl`](scripts/S-401-plotcompositeness.jl) compute compositeness from r/a, plot X and Z.
   - [`S-402-effective_range_with_fit.jl`](scripts/S-402-effective_range_with_fit.jl) fit a⁻¹, r, and N parameters
   - [`S-403-effect-of-pion-exchange.jl`](scripts/S-403-effect-of-pion-exchange.jl) comparison of the (Dˣ⁺ + Dˣ⁺) model with the (one-Dˣ⁺) model.
   - [`S-404-effective_range_cusps_in_R.jl`](scripts/S-404-effective_range_cusps_in_R.jl) effect of grid discrimination on the value of the eff. range
   - [`S-405-validate-effrange-expansion.jl`](scripts/S-405-validate-effrange-expansion.jl) plot EFE over the amplitude (reuse S-400 output)
   - [`N-406-symbolic-r.jl`](notebooks/N-406-symbolic-r.jl) evaluate effective range symbolically for two-body coupled channels.
   - [`S-407-fit-intervals.jl`](scripts/S-407-fit-intervals.jl) Fit model with OPE in different intervals
   - [`S-408-second-trhreshold.jl`](scripts/S-408-second-trhreshold.jl) Effective range for the second threshold
   - [`S-409-higher-term-expansion.jl`](scripts/S-409-higher-term-expansion.jl) Effective range for the second threshold

6. **Additional plots for talks**

   - [`S-900-howwidthismade.jl`](scripts/S-900-howwidthismade.jl) contributions to the width
   - [`S-901-xxx.jl`](scripts/S-901-xxx.jl) DDpi momentum spectrum. Femtoscopy
   - [`N-902-triangle_singularity.jl`](notebooks/N-902-triangle_singularity.jl) Position of the traingle singularity in the complex energy plane

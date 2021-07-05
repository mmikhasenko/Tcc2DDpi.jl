
# masses
const mπ⁰ = 0.1349768
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mD⁺ = 1.86965
# 
const mDˣ⁺ = 2.01026; const ΓDˣ⁺ = 83.4e-6
const mDˣ⁰ = 2.00685; const ΓDˣ⁰ = 55.2e-6
const mγ = 0.0
# 
const μDˣ⁺D⁰ = mDˣ⁺*mD⁰ / (mDˣ⁺+mD⁰)

const fm_times_mev = 197.3269804

# coupling constants
const h² = 20.13e-3
const f² = 282.42
const μ₀ = -3.77
const μ₊ = 1.0

# conversion functinos
e2m(e) = (mD⁰+mDˣ⁺)+e*1e-3
m2e(m) = (m-mD⁰-mDˣ⁺)*1e3


const ΔNLL_90CL = 1.352
const ΔNLL_95CL = 1.921
const ΔNLL_68CL = 0.5

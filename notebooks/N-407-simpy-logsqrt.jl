### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ 15be4d9e-4d5b-11ee-3721-c7a9a885bcb2
using SymPy

# ╔═╡ 9157d5bf-7f1d-437f-80a4-38930e70fdd4
md"""
# Teylor for the LogSqrt example

The notebook demonstrate the teylor expantion in powers of √x of a function with logarithmic branch point.

$f(x) = \log(\sqrt x+1)$

The function has a sqrt branch point at zero, and two sheet associated to it.
The logarithmic branch point is located on the secod sheet at x=1.

The function can be expanded at x=0

$f(x) = \sum_0^\infty r_i x^i + \sqrt x \sum_{i=0}^\infty n_i x^i,\,\,\text{ for } x<1$

The convergence is limited by the position of the logarithmic branch point.

The expansion coefficient are calculated by circular integral around the sqrt branch point.

$I_t(\epsilon) = \oint_{|r|=\epsilon} \frac{f(x)}{x^t} \mathrm{d}x\,,$

The integrals can be related to the coeffecient of the expansion, up to $\sqrt{\epsilon}$.

```math
\begin{align}
n_0 &\leftarrow I_0  + O(\sqrt{\epsilon})\,,\\
r_0 &\leftarrow \int \frac{f(x)-n_0 \sqrt{x}}{x} \mathrm{d} x = I_1 + O(\sqrt{\epsilon})\\
r_1 &\leftarrow \int \frac{f(x)-n_0 \sqrt{x}}{x^2} \mathrm{d} x = I_2  + O(\sqrt{\epsilon})\,,\\
\end{align}
```

The coefficient of the effective range expansion,

$f = N(1/a+r x - i \sqrt{x}) + o(x)$

are related to the taylor series coefficient as,

```math
\begin{align}
N &= i n_0\,,\\
1/a &= r_0 / n_0\,,\\
r &= r_1 / n_0\,. 
\end{align}
```

Below, we find how these scale with $\epsilon$
"""

# ╔═╡ fc75fb06-afa0-430e-9ce4-5a20cc2f0a5c
begin
	@syms x
	f(x) = log(sqrt(x)+1)
end

# ╔═╡ 82a9f81b-562c-4526-9956-6c8a9e79a141
md"""
## Integrals
"""

# ╔═╡ 7d9e676c-3f9d-4d8f-af33-0755f2cd8957
f(x)

# ╔═╡ acbdefd0-d93a-4eec-a904-66bfe153eb40
sympy.Integral(f(x),x).doit()

# ╔═╡ fc5f1c6e-2981-4e63-b90a-1edcf1413de7
sympy.Integral(f(x)/x,x).doit()

# ╔═╡ 8ddcf343-61ae-4f0d-af4a-6d56a524e183
sympy.Integral(f(x)/x^2,x).doit() |> factor# |> simplify |> numerator

# ╔═╡ 7a499a42-c9f5-4d7b-92c2-af451b013b24
md"""
Mathematica gives [nicer expression](https://www.wolframalpha.com/input?i=integral+of+log%28sqrt%28x%29%2B1%29%2Fx%5E2). The discontinuity reads:
```
F_2 => 2/sqrt(x) + (x-1)/x*(log(1+sqrt(x))-log(1-sqrt(x)))
```
"""

# ╔═╡ 0711e3e4-4b31-4e45-bec2-21ede7b18201
md"""
## Discontinuity

is computed on two sides of the branch cut. For the upper(right) side, one takes the sqrt function with the + sign, for the opposite part, the sign flips, $-\sqrt x$

$\text{Disc} F = F(\sqrt x) - F(-\sqrt x)$

**Important:** the integral is equal to

$I = \frac{- \text{Disc}\,F}{2 \pi i}$

due to the direction of the integration (anti-clockwise).
"""

# ╔═╡ d3624fc6-5749-4dfb-b7ba-0945a5a849b8
md"""
#### Normalization integral, $I_0$
"""

# ╔═╡ d6e0a3cd-61eb-4aa4-bafe-adff8f2326ee
DiscN = 2sqrt(x) + (x-1)*(log(1+sqrt(x))-log(1-sqrt(x)))

# ╔═╡ 921de76d-d94f-438d-89b1-fe5fb1749fff
DiscN0 = 4*x*sqrt(x)/3

# ╔═╡ 3a047bca-6299-4df0-a018-1c5f0eedf447
N_over_i = DiscN / DiscN0

# ╔═╡ 2e131b50-4a4d-4afa-91c5-edcb700ce414
series(N_over_i, sqrt(x), 0, 5) 

# ╔═╡ 3afea485-b0a4-4bd0-9883-0080bb0fb766
md"""
#### Effective range integral, $I_2$

We split the expression for $r$ into three terms,

$r = 2i (\text{Disc}_R + \text{Corr}_N + V) / N$

- the $\text{Disc}_R$: comes from the integral of $(f-\sqrt x)/x^2 \to I_2 + 2/\sqrt{x}$
- the correction $\text{Corr}_N$: comes from the mismatch of $N$ and $1$. It is equal to $2(n_0 - 1)  /\sqrt x$.
"""

# ╔═╡ ab495608-e89b-47ae-b50b-2c6a4fed4cf8
begin
	corrN = 4/sqrt(x)*(DiscN / DiscN0 - 1)
	series(corrN, x, 0, 3)
end

# ╔═╡ d08b93c5-6e90-4672-ad00-32d6c5ecc635
begin
	DiscR = 2/sqrt(x) + (x-1)/x*(log(1+sqrt(x))-log(1-sqrt(x)))
	series(DiscR, sqrt(x), 0, 5) 
end

# ╔═╡ 1ef01a61-cb5b-49e6-b124-c0d47a1eddb9
begin
	const πi = sympy.pi*Sym(1im)
	V = 2πi/2
end;

# ╔═╡ 6c97f079-0a55-47f7-9546-0668198fb9e2
Disc = V + DiscR + corrN ;

# ╔═╡ c4dcd9dc-99c4-4e11-b670-277c08946e26
r  = 2*Disc/(2sympy.pi)/N_over_i; 

# ╔═╡ 9b6e1598-7144-4a3f-9b07-3ac4fbd8ab31
series(r, sqrt(x), 0, 5) 

# ╔═╡ 71a3c42f-c60a-4a7d-ac41-0b27e2061217
md"""
## Result

$r = i$ is the correct result. The following expansion term is $-32\sqrt{\epsilon}/15\pi$!
"""

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─9157d5bf-7f1d-437f-80a4-38930e70fdd4
# ╠═15be4d9e-4d5b-11ee-3721-c7a9a885bcb2
# ╠═fc75fb06-afa0-430e-9ce4-5a20cc2f0a5c
# ╟─82a9f81b-562c-4526-9956-6c8a9e79a141
# ╠═7d9e676c-3f9d-4d8f-af33-0755f2cd8957
# ╠═acbdefd0-d93a-4eec-a904-66bfe153eb40
# ╠═fc5f1c6e-2981-4e63-b90a-1edcf1413de7
# ╠═8ddcf343-61ae-4f0d-af4a-6d56a524e183
# ╟─7a499a42-c9f5-4d7b-92c2-af451b013b24
# ╟─0711e3e4-4b31-4e45-bec2-21ede7b18201
# ╟─d3624fc6-5749-4dfb-b7ba-0945a5a849b8
# ╠═d6e0a3cd-61eb-4aa4-bafe-adff8f2326ee
# ╠═921de76d-d94f-438d-89b1-fe5fb1749fff
# ╠═3a047bca-6299-4df0-a018-1c5f0eedf447
# ╠═2e131b50-4a4d-4afa-91c5-edcb700ce414
# ╟─3afea485-b0a4-4bd0-9883-0080bb0fb766
# ╠═ab495608-e89b-47ae-b50b-2c6a4fed4cf8
# ╠═d08b93c5-6e90-4672-ad00-32d6c5ecc635
# ╠═1ef01a61-cb5b-49e6-b124-c0d47a1eddb9
# ╠═6c97f079-0a55-47f7-9546-0668198fb9e2
# ╠═c4dcd9dc-99c4-4e11-b670-277c08946e26
# ╠═9b6e1598-7144-4a3f-9b07-3ac4fbd8ab31
# ╟─71a3c42f-c60a-4a7d-ac41-0b27e2061217

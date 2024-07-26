### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 282214f4-2e36-4fd0-bc20-689b668e0c7f
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ ecf16049-a279-4f8b-8e85-7b59b0d09037
using BioSequences: LongRNA

# ╔═╡ 0c8d61cd-1b2d-42c4-bb54-574e0936e453
using Distributions: Gamma

# ╔═╡ a60b60ea-3f73-4044-b0be-794bbfc67c7b
using Makie: @L_str

# ╔═╡ 2f0e747f-1388-4a7c-86a6-46c3e4b49caa
using NaNStatistics: nanmean

# ╔═╡ c641f453-7622-4b5f-99c3-d19763f46dfb
using NaNStatistics: nanstd

# ╔═╡ 2f78dd2a-3112-4043-902b-9e8452f26e25
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 7f3e3858-7a2f-4e96-ab20-5ecf3e1f6950
using Statistics: cor

# ╔═╡ 401e15f2-c991-40ce-a066-849534c7b6d9
using Statistics: mean

# ╔═╡ 4c66655e-c41a-4c8e-a303-295c2418cd44
md"""
# Imports
"""

# ╔═╡ f9e4aacf-d826-4c91-ae1e-b6a3f88a14b8
import PlutoUI

# ╔═╡ a552d4fb-d7b8-4c5d-8bbe-62e1c97487d2
import CairoMakie

# ╔═╡ 0a39d5a9-e1b9-4916-a97c-c1db8e485616
import Makie

# ╔═╡ 42979e48-f4dc-4605-9e63-192f82bea15a
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ b58239bb-12ea-4263-a24c-77ad1dd54732
import SamApp2024

# ╔═╡ f9044b51-72f9-458f-878f-88cbd590a271
import StatsBase

# ╔═╡ dd74de67-3877-40be-941b-e01861966149
PlutoUI.TableOfContents()

# ╔═╡ 175edee8-046a-4f00-8a8c-ab38b737a130
md"""
# Figures
"""

# ╔═╡ 3279530f-d231-438e-a033-83595cdac31d


# ╔═╡ Cell order:
# ╠═4c66655e-c41a-4c8e-a303-295c2418cd44
# ╠═282214f4-2e36-4fd0-bc20-689b668e0c7f
# ╠═f9e4aacf-d826-4c91-ae1e-b6a3f88a14b8
# ╠═a552d4fb-d7b8-4c5d-8bbe-62e1c97487d2
# ╠═0a39d5a9-e1b9-4916-a97c-c1db8e485616
# ╠═42979e48-f4dc-4605-9e63-192f82bea15a
# ╠═b58239bb-12ea-4263-a24c-77ad1dd54732
# ╠═f9044b51-72f9-458f-878f-88cbd590a271
# ╠═ecf16049-a279-4f8b-8e85-7b59b0d09037
# ╠═0c8d61cd-1b2d-42c4-bb54-574e0936e453
# ╠═a60b60ea-3f73-4044-b0be-794bbfc67c7b
# ╠═2f0e747f-1388-4a7c-86a6-46c3e4b49caa
# ╠═c641f453-7622-4b5f-99c3-d19763f46dfb
# ╠═2f78dd2a-3112-4043-902b-9e8452f26e25
# ╠═7f3e3858-7a2f-4e96-ab20-5ecf3e1f6950
# ╠═401e15f2-c991-40ce-a066-849534c7b6d9
# ╠═dd74de67-3877-40be-941b-e01861966149
# ╠═175edee8-046a-4f00-8a8c-ab38b737a130
# ╠═3279530f-d231-438e-a033-83595cdac31d

### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f08d8cfc-711e-47c7-a5a9-f70e3c0cb1ac
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 35b53fe4-c79d-48b0-bbe6-b8c6d469c3c0
using BioSequences: LongRNA

# ╔═╡ b738e180-90d8-419b-95bd-c383e4b19012
using DataFrames: DataFrame

# ╔═╡ 4c6704f4-8802-4f5c-80c5-25d8520fcc25
using Makie: @L_str

# ╔═╡ 7f4e47e0-980e-440c-9e81-ea25a0149f34
using NaNStatistics: nansum

# ╔═╡ 8ff5bd07-2672-4565-b14f-deadc0a3b2d3
using Random: bitrand

# ╔═╡ c4f91ffb-cc6e-4347-804a-887ca45be249
using Statistics: cor

# ╔═╡ a718fd95-2963-4b35-bddd-d3b53a2faf52
using Statistics: mean

# ╔═╡ 4244cffe-81f2-41ff-950d-97e973947a32
using StatsBase: countmap

# ╔═╡ d42dae48-fdc2-11ef-191c-b7b143d538a5
md"# Imports"

# ╔═╡ dd8fcacd-6feb-439f-8495-861731ef26db
import Makie

# ╔═╡ a80399d9-96e5-45a1-8d0c-0ab364781e3e
import CairoMakie

# ╔═╡ 3553a742-fbfb-4553-abe0-17a4044f6eae
import FASTX

# ╔═╡ 059b139e-3d06-41d5-a89a-6990db4ab9fa
import Infernal

# ╔═╡ e71e2863-3a8e-4854-812e-ec14b5458350
import SamApp2024

# ╔═╡ 57681d3b-7d23-4b9c-ab1c-d24b72845654
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 4aa377c9-05fc-45df-9ab2-0500952b41f7
import Rfam

# ╔═╡ 59e38478-1323-4047-b591-c0a788e898c8
import PlutoUI

# ╔═╡ 41fb8639-2472-4341-98a1-b8f11475e31e
import Unitful

# ╔═╡ 030731fb-dd82-4d0f-a9d7-026af8378d20
import ViennaRNA

# ╔═╡ 1c2478bc-bdfc-4596-a053-b5736cbe6504
import StatsBase

# ╔═╡ 6bc48dfb-9009-4aee-8f47-3623b8834ee7
import KernelDensity

# ╔═╡ 6155ff9f-af5f-4c1a-8e6e-965f11efd285
PlutoUI.TableOfContents()

# ╔═╡ e0367f59-2f02-4a09-8018-e71317695f3b
md"# Load data"

# ╔═╡ aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
dms_data = SamApp2024.load_dms_data_20250303()

# ╔═╡ Cell order:
# ╠═d42dae48-fdc2-11ef-191c-b7b143d538a5
# ╠═f08d8cfc-711e-47c7-a5a9-f70e3c0cb1ac
# ╠═dd8fcacd-6feb-439f-8495-861731ef26db
# ╠═a80399d9-96e5-45a1-8d0c-0ab364781e3e
# ╠═3553a742-fbfb-4553-abe0-17a4044f6eae
# ╠═059b139e-3d06-41d5-a89a-6990db4ab9fa
# ╠═e71e2863-3a8e-4854-812e-ec14b5458350
# ╠═57681d3b-7d23-4b9c-ab1c-d24b72845654
# ╠═4aa377c9-05fc-45df-9ab2-0500952b41f7
# ╠═59e38478-1323-4047-b591-c0a788e898c8
# ╠═41fb8639-2472-4341-98a1-b8f11475e31e
# ╠═030731fb-dd82-4d0f-a9d7-026af8378d20
# ╠═1c2478bc-bdfc-4596-a053-b5736cbe6504
# ╠═6bc48dfb-9009-4aee-8f47-3623b8834ee7
# ╠═35b53fe4-c79d-48b0-bbe6-b8c6d469c3c0
# ╠═b738e180-90d8-419b-95bd-c383e4b19012
# ╠═4c6704f4-8802-4f5c-80c5-25d8520fcc25
# ╠═7f4e47e0-980e-440c-9e81-ea25a0149f34
# ╠═8ff5bd07-2672-4565-b14f-deadc0a3b2d3
# ╠═c4f91ffb-cc6e-4347-804a-887ca45be249
# ╠═a718fd95-2963-4b35-bddd-d3b53a2faf52
# ╠═4244cffe-81f2-41ff-950d-97e973947a32
# ╠═6155ff9f-af5f-4c1a-8e6e-965f11efd285
# ╠═e0367f59-2f02-4a09-8018-e71317695f3b
# ╠═aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d

### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 91f4edd8-290d-4270-83c1-f7c6281e9f68
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ b33c48ad-c1a2-4393-a711-534378768f36
using BioSequences: LongRNA

# ╔═╡ d1863573-2e64-4cc4-a7c7-153d891a1420
using DataFrames: DataFrame

# ╔═╡ 770b8573-ee3f-47ed-a2f4-1699765b7081
using Distributions: Gamma

# ╔═╡ e51de4c5-db65-46ee-ae7e-419a8af4d178
using Distributions: logpdf

# ╔═╡ 960f924a-046d-4834-a00f-4b45f4449232
using Distributions: pdf

# ╔═╡ d80f3a16-6dea-4816-8d1e-2346e3efc34c
using Distributions: Poisson

# ╔═╡ b3707ee7-4f35-4f78-862d-bdef6e08ad14
using LinearAlgebra: Diagonal

# ╔═╡ aea74daa-6d3d-4c2c-8285-c541e1375772
using LinearAlgebra: eigen

# ╔═╡ 98c030fe-e0f2-4949-90e6-da7dd9b9d245
using Makie: @L_str

# ╔═╡ bd112693-84c9-4f18-927f-b307319006d0
using NaNStatistics: nansum

# ╔═╡ 0c7ec2d1-f6fe-4a95-bc81-e320666d880c
using NaNStatistics: nanmean

# ╔═╡ 2665e9c5-2e4f-41ec-a9f7-58b3b96a9ff3
using NaNStatistics: nanstd

# ╔═╡ 9299b1d4-a573-4808-9e26-396c2e04b8cc
using Random: bitrand

# ╔═╡ 530dbe76-e3a2-4658-ab86-daf6588cf372
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 0c2a4cd9-a42a-42fc-b85b-46e5987a77cc
using Statistics: cor

# ╔═╡ 8da6f416-1a00-4cc8-a9c0-cdd76ec364b0
using Statistics: mean

# ╔═╡ d705aa70-ad0e-47aa-a31b-714984d3b65c
using StatsBase: countmap

# ╔═╡ dba2cd0a-5edb-4581-970b-c8d7d84cd331
import CairoMakie

# ╔═╡ 90bb26b9-08ca-4de7-a383-f9d2e2c200f1
import CSV

# ╔═╡ 870f93fb-1969-4654-abd3-6aedf1215cec
import FASTX

# ╔═╡ 8bba31d0-4abe-4c3a-a776-0355f256fa29
import HDF5

# ╔═╡ a3f85e35-0381-4e8d-8442-bfb177da795f
import Infernal

# ╔═╡ b0d27526-33de-40f5-bb91-27f36154816c
import KernelDensity

# ╔═╡ 89b54d54-f2a9-4cca-93aa-756042ccb8b2
import Makie

# ╔═╡ 28c40841-a3cf-42be-aac4-9076c881e6de
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 1813f520-ca4b-4675-8520-89b2be2513a7
import Rfam

# ╔═╡ 3597a655-3b5e-47ef-8731-0447cc7cfac3
import SamApp2024

# ╔═╡ 67f6774a-0417-49ba-a7fc-27b6dd543d38
import StatsBase

# ╔═╡ dcf0701e-3fcc-4199-8285-803305ee38b2
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ d1387ef6-6497-4f0b-873c-4de5f8aa0715
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 0297ebf2-0822-4ae0-b4e1-73b8ca1c4452
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ Cell order:
# ╠═91f4edd8-290d-4270-83c1-f7c6281e9f68
# ╠═dba2cd0a-5edb-4581-970b-c8d7d84cd331
# ╠═90bb26b9-08ca-4de7-a383-f9d2e2c200f1
# ╠═870f93fb-1969-4654-abd3-6aedf1215cec
# ╠═8bba31d0-4abe-4c3a-a776-0355f256fa29
# ╠═a3f85e35-0381-4e8d-8442-bfb177da795f
# ╠═b0d27526-33de-40f5-bb91-27f36154816c
# ╠═89b54d54-f2a9-4cca-93aa-756042ccb8b2
# ╠═28c40841-a3cf-42be-aac4-9076c881e6de
# ╠═1813f520-ca4b-4675-8520-89b2be2513a7
# ╠═3597a655-3b5e-47ef-8731-0447cc7cfac3
# ╠═67f6774a-0417-49ba-a7fc-27b6dd543d38
# ╠═b33c48ad-c1a2-4393-a711-534378768f36
# ╠═d1863573-2e64-4cc4-a7c7-153d891a1420
# ╠═770b8573-ee3f-47ed-a2f4-1699765b7081
# ╠═e51de4c5-db65-46ee-ae7e-419a8af4d178
# ╠═960f924a-046d-4834-a00f-4b45f4449232
# ╠═d80f3a16-6dea-4816-8d1e-2346e3efc34c
# ╠═b3707ee7-4f35-4f78-862d-bdef6e08ad14
# ╠═aea74daa-6d3d-4c2c-8285-c541e1375772
# ╠═98c030fe-e0f2-4949-90e6-da7dd9b9d245
# ╠═bd112693-84c9-4f18-927f-b307319006d0
# ╠═0c7ec2d1-f6fe-4a95-bc81-e320666d880c
# ╠═2665e9c5-2e4f-41ec-a9f7-58b3b96a9ff3
# ╠═9299b1d4-a573-4808-9e26-396c2e04b8cc
# ╠═530dbe76-e3a2-4658-ab86-daf6588cf372
# ╠═0c2a4cd9-a42a-42fc-b85b-46e5987a77cc
# ╠═8da6f416-1a00-4cc8-a9c0-cdd76ec364b0
# ╠═d705aa70-ad0e-47aa-a31b-714984d3b65c
# ╠═dcf0701e-3fcc-4199-8285-803305ee38b2
# ╠═d1387ef6-6497-4f0b-873c-4de5f8aa0715
# ╠═0297ebf2-0822-4ae0-b4e1-73b8ca1c4452

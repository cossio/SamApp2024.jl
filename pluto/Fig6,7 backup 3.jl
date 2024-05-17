### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 77855f93-2b64-45bc-a307-df7f6e6187b3
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
using BioSequences: LongRNA

# ╔═╡ 55735555-c2f7-41f0-becd-cbad5717d7be
using DataFrames: DataFrame

# ╔═╡ eff20728-e3ee-470a-afc0-3e26a30c11be
using Distributions: Gamma

# ╔═╡ 746ffc2b-1be1-4f10-8d79-12364a8a8251
using Distributions: logpdf

# ╔═╡ b0997272-f623-43d8-a5ff-c1b1d6be99a4
using Distributions: pdf

# ╔═╡ a1807f68-6f16-4114-98ab-b075b45b7201
using Distributions: Poisson

# ╔═╡ 483d3154-5222-40b8-a380-aea037e4618c
using LinearAlgebra: Diagonal

# ╔═╡ 2132164b-51c2-4415-9d07-022a24011f3d
using LinearAlgebra: eigen

# ╔═╡ fb8e8cbd-050d-4d0e-81ac-32528f628a0e
using Makie: @L_str

# ╔═╡ 8973e48f-df81-4749-b71d-aea5ac4614b3
using NaNStatistics: nanmean

# ╔═╡ 9141d7f1-194b-4199-b954-9cd0e13b2f15
using NaNStatistics: nansum

# ╔═╡ 96e03cbd-dc6f-4923-bb56-106ca66c1867
using Random: bitrand

# ╔═╡ e73b3886-162a-4a77-a28b-cdf269109b98
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2abcc581-8722-4cf7-bc09-8bf98a9b8648
using Statistics: cor

# ╔═╡ 6d6e3dc8-28eb-496d-9622-8128422ed2ed
using Statistics: mean

# ╔═╡ 7edc30fd-a2d6-4818-855c-c7f69c9f589b
using Statistics: middle

# ╔═╡ 4974c2e2-058d-41ca-924d-16709e4a58e6
using StatsBase: countmap

# ╔═╡ 94d99837-7415-4acb-b5d3-3b1dec5af05e
using Unitful: ustrip

# ╔═╡ 45907f4d-29ec-4be7-97e8-bfcb4695416b
import CairoMakie

# ╔═╡ f743167e-8552-4002-b9ca-c995cf3e9829
import CSV

# ╔═╡ 6be9532a-762f-4830-8d7f-175545fe75c9
import FASTX

# ╔═╡ 02b62a0f-159c-470d-91d0-673898ca277b
import HDF5

# ╔═╡ 651a9155-c7fc-4858-ae65-b433d1d116cc
import Infernal

# ╔═╡ c6371fcc-bb84-4ae4-a2d8-972b362c5350
import KernelDensity

# ╔═╡ 297c712c-e7e6-4d2f-a9c2-59b29665e6e3
import Makie

# ╔═╡ c53d715e-40d3-44cb-b85a-9c4c61d99819
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ f513852e-ce6a-41f0-98aa-155e06244aaf
import Rfam

# ╔═╡ ca499e53-296f-4cf3-9df1-5070e22fd6f0
import SamApp2024

# ╔═╡ 15974106-a5c6-4867-acfb-cdc69f5a57be
import StatsBase

# ╔═╡ b471828b-0822-4a08-a63c-5fa01e8d90b2
import ViennaRNA

# ╔═╡ 13f58dca-f8c6-40d9-bf2c-d20995181775
import ViennaRNA_jll

# ╔═╡ Cell order:
# ╠═77855f93-2b64-45bc-a307-df7f6e6187b3
# ╠═45907f4d-29ec-4be7-97e8-bfcb4695416b
# ╠═f743167e-8552-4002-b9ca-c995cf3e9829
# ╠═6be9532a-762f-4830-8d7f-175545fe75c9
# ╠═02b62a0f-159c-470d-91d0-673898ca277b
# ╠═651a9155-c7fc-4858-ae65-b433d1d116cc
# ╠═c6371fcc-bb84-4ae4-a2d8-972b362c5350
# ╠═297c712c-e7e6-4d2f-a9c2-59b29665e6e3
# ╠═c53d715e-40d3-44cb-b85a-9c4c61d99819
# ╠═f513852e-ce6a-41f0-98aa-155e06244aaf
# ╠═ca499e53-296f-4cf3-9df1-5070e22fd6f0
# ╠═15974106-a5c6-4867-acfb-cdc69f5a57be
# ╠═b471828b-0822-4a08-a63c-5fa01e8d90b2
# ╠═13f58dca-f8c6-40d9-bf2c-d20995181775
# ╠═c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
# ╠═55735555-c2f7-41f0-becd-cbad5717d7be
# ╠═eff20728-e3ee-470a-afc0-3e26a30c11be
# ╠═746ffc2b-1be1-4f10-8d79-12364a8a8251
# ╠═b0997272-f623-43d8-a5ff-c1b1d6be99a4
# ╠═a1807f68-6f16-4114-98ab-b075b45b7201
# ╠═483d3154-5222-40b8-a380-aea037e4618c
# ╠═2132164b-51c2-4415-9d07-022a24011f3d
# ╠═fb8e8cbd-050d-4d0e-81ac-32528f628a0e
# ╠═8973e48f-df81-4749-b71d-aea5ac4614b3
# ╠═9141d7f1-194b-4199-b954-9cd0e13b2f15
# ╠═96e03cbd-dc6f-4923-bb56-106ca66c1867
# ╠═e73b3886-162a-4a77-a28b-cdf269109b98
# ╠═2abcc581-8722-4cf7-bc09-8bf98a9b8648
# ╠═6d6e3dc8-28eb-496d-9622-8128422ed2ed
# ╠═7edc30fd-a2d6-4818-855c-c7f69c9f589b
# ╠═4974c2e2-058d-41ca-924d-16709e4a58e6
# ╠═94d99837-7415-4acb-b5d3-3b1dec5af05e

### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 388ad1ca-19c9-4407-87cf-e3cb2f5553c2
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 5ed177a0-2ac4-41ac-a3c6-3610382b448f
using BioSequences: LongRNA

# ╔═╡ acaf92ee-3f1c-42a7-a7ef-f350383efd22
using DataFrames: DataFrame

# ╔═╡ e5dec066-bf0f-4e1f-95da-8cc088100a72
using Distributions: Gamma

# ╔═╡ 875ac2a1-9f75-4af7-88c4-c67f01816ede
using Distributions: logpdf

# ╔═╡ 1799ff95-f532-454a-9a8e-fc7724034f5a
using Distributions: pdf

# ╔═╡ d3cc6cb0-46e3-4dc2-9d19-084fb63357e7
using Distributions: Poisson

# ╔═╡ 626193e6-c344-48f6-a18c-d00aeed27795
using LinearAlgebra: Diagonal

# ╔═╡ db1673d0-5c89-4438-b719-e38486ee2a38
using LinearAlgebra: eigen

# ╔═╡ e4d9e6a6-5c42-4c37-a34d-68e493788e59
using Makie: @L_str

# ╔═╡ 134a1e15-1125-45bc-9ff6-41f49548fe0e
using NaNStatistics: nansum

# ╔═╡ ae6295d4-c810-4ea4-af5e-4d25ad98f7d1
using Random: bitrand

# ╔═╡ 10e33667-2e6a-4d8a-884c-d43cf50232c7
using Statistics: cor

# ╔═╡ aec5a29a-ce77-4977-a15b-4cd2bf973ff4
using Statistics: mean

# ╔═╡ 7035958e-0b50-4736-9a1c-ca42145cf73d
using StatsBase: countmap

# ╔═╡ 9fc798b2-aa30-4c84-b19a-f425ce5f25d5
import CairoMakie

# ╔═╡ 327dc207-cc7f-4ed9-9f1c-7da5b32cef98
import CSV

# ╔═╡ 1469e5a8-0303-4ebc-b7ef-ec86423572c8
import HDF5

# ╔═╡ 7fb99cac-a2f3-4d9b-b1a5-bf95577878ad
import FASTX

# ╔═╡ cd5f8b04-582f-4376-acfa-8a901e385dbb
import Infernal

# ╔═╡ 685e0d41-b639-40b5-96e8-2d8c56b7b9a0
import KernelDensity

# ╔═╡ 3430b1d4-9605-48d6-b1e9-0f0702214817
import Makie

# ╔═╡ bcf1ddb2-74be-44e8-809f-66b7f27855a3
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ eb8830a1-3ea2-4757-9a2b-5a19249baae0
import Rfam

# ╔═╡ d7c904eb-c1ab-4f72-a0b9-3c53e0b86a59
import SamApp2024

# ╔═╡ f1e092ee-12ac-4680-9cb0-80a7fec907a6
import StatsBase

# ╔═╡ 1ae2d1fc-db46-4c34-a445-82c86953385c
import Logomaker

# ╔═╡ 0a8dbd37-3697-47c9-8a1e-b8017978082e
# CM model from Rfam (this has the noisy floor!)
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162");

# ╔═╡ c856268c-a2e2-43ef-ab74-f76494ef235e
# trimmed (no inserts) aligned fasta
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");
# these are already aligned and without inserts

# ╔═╡ 9cf5ab26-52db-46e5-b582-8d3d8280fe7b
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ aa150c41-9322-4fbd-9221-29f747f94b2c
begin
	# APC-corrected contact map
	Cmap = SamApp2024.effective_contacts(SamApp2024.rbm2022(), SamApp2024.onehot(RF00162_hits_sequences));
	
	# zero-out diagonal
	for i in axes(Cmap, 1)
	    Cmap[i,i] = 0
	end
end

# ╔═╡ 805a7767-6ee3-4414-8a49-a272fdbcf60d
begin
	plotted_matrix = copy(Cmap)
	for i = axes(plotted_matrix, 1), j = axes(plotted_matrix, 2)
	    if i < j
	        plotted_matrix[i,j] = Cmap[i,j]
	    else
	        plotted_matrix[i,j] = 0
	    end
	end
end

# ╔═╡ fb8646db-e078-4ad0-92b9-54ad63dcffed
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1], yreversed=true, width=400, height=400, xticks=0:10:size(Cmap,1), yticks=0:10:size(Cmap,2), xlabel="sequence position", ylabel="sequence position")
	hm = Makie.heatmap!(ax, plotted_matrix, colormap=["white", "white", "gray", "dimgray", "black", "black"])
	Makie.hidespines!(ax, :t, :r)
	Makie.Colorbar(fig[1,2], hm, height=200, label="Epistasis score")
	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Contacts bw.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═388ad1ca-19c9-4407-87cf-e3cb2f5553c2
# ╠═9fc798b2-aa30-4c84-b19a-f425ce5f25d5
# ╠═327dc207-cc7f-4ed9-9f1c-7da5b32cef98
# ╠═1469e5a8-0303-4ebc-b7ef-ec86423572c8
# ╠═7fb99cac-a2f3-4d9b-b1a5-bf95577878ad
# ╠═cd5f8b04-582f-4376-acfa-8a901e385dbb
# ╠═685e0d41-b639-40b5-96e8-2d8c56b7b9a0
# ╠═3430b1d4-9605-48d6-b1e9-0f0702214817
# ╠═bcf1ddb2-74be-44e8-809f-66b7f27855a3
# ╠═eb8830a1-3ea2-4757-9a2b-5a19249baae0
# ╠═d7c904eb-c1ab-4f72-a0b9-3c53e0b86a59
# ╠═f1e092ee-12ac-4680-9cb0-80a7fec907a6
# ╠═1ae2d1fc-db46-4c34-a445-82c86953385c
# ╠═5ed177a0-2ac4-41ac-a3c6-3610382b448f
# ╠═acaf92ee-3f1c-42a7-a7ef-f350383efd22
# ╠═e5dec066-bf0f-4e1f-95da-8cc088100a72
# ╠═875ac2a1-9f75-4af7-88c4-c67f01816ede
# ╠═1799ff95-f532-454a-9a8e-fc7724034f5a
# ╠═d3cc6cb0-46e3-4dc2-9d19-084fb63357e7
# ╠═626193e6-c344-48f6-a18c-d00aeed27795
# ╠═db1673d0-5c89-4438-b719-e38486ee2a38
# ╠═e4d9e6a6-5c42-4c37-a34d-68e493788e59
# ╠═134a1e15-1125-45bc-9ff6-41f49548fe0e
# ╠═ae6295d4-c810-4ea4-af5e-4d25ad98f7d1
# ╠═10e33667-2e6a-4d8a-884c-d43cf50232c7
# ╠═aec5a29a-ce77-4977-a15b-4cd2bf973ff4
# ╠═7035958e-0b50-4736-9a1c-ca42145cf73d
# ╠═0a8dbd37-3697-47c9-8a1e-b8017978082e
# ╠═c856268c-a2e2-43ef-ab74-f76494ef235e
# ╠═9cf5ab26-52db-46e5-b582-8d3d8280fe7b
# ╠═aa150c41-9322-4fbd-9221-29f747f94b2c
# ╠═805a7767-6ee3-4414-8a49-a272fdbcf60d
# ╠═fb8646db-e078-4ad0-92b9-54ad63dcffed

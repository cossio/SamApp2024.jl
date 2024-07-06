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

# ╔═╡ e743ce81-29e8-4bac-b22e-79e123db3326
using NaNStatistics: nancor

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

# ╔═╡ d443b9aa-7449-4447-983c-18263e61be12
md"""
# Imports
"""

# ╔═╡ 59e3b0b4-c955-4985-bc79-bbc9220539a6
import PlutoUI

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

# ╔═╡ 18452b0f-3c62-492c-a41b-6a9043ce86b1
PlutoUI.TableOfContents()

# ╔═╡ 4206fb3f-99c4-4e4b-bbc9-a97891a2f321
md"""
# Load data
"""

# ╔═╡ dcf0701e-3fcc-4199-8285-803305ee38b2
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ d1387ef6-6497-4f0b-873c-4de5f8aa0715
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 0297ebf2-0822-4ae0-b4e1-73b8ca1c4452
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ c011571f-74e7-4bae-be5e-918ad678397c
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ e39876f7-5399-4e4a-a1d5-101355741f28
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ d527106f-7e5b-49a5-930a-fd847a782542
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 3ac302e6-b7c6-474d-af2b-c12420ebd0be
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 64cb3148-1d0a-4f7d-8859-ded1e2fc5283
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ f9b70f45-5ce5-460a-b7bd-46a9cca2c8f3
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 4e123003-aede-4e18-956f-7f87a85d5083
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 6bdf0eb7-511a-4900-bc51-66e960cd391c
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ 6994d41d-0fdf-4dd6-bc8e-122d745ad34c
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 44042f22-0ae9-44a4-ad22-7d316e976b92
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 29896cb1-147d-4e15-a815-a2a55e7ce483
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 7a266977-5fea-4a23-84d0-e967650cc91e
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 4762b1ab-4dc9-4735-a9b3-6ab26a6c19d7
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 143dbfef-3376-4892-998a-a769acb78a05
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ 0fd1db6a-c51e-4e18-85aa-3064103b345f
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ 880f1579-fae9-4cb4-9fd5-b533203fa0e0
_rbmlo = rbm_seqs ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .< -300));

# ╔═╡ 14da7dad-da45-457f-ac4a-3a59dabf0898
_rbmhi = rbm_seqs ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .> -300));

# ╔═╡ d71ccb35-7881-4571-8c5c-bac572630c37
ΔR_sam = (
    nanmean(shape_data_rep0.shape_reactivities[:, :, conds_sam_rep0]; dim=3) .- 
    shape_data_rep0.shape_reactivities[:, :, only(conds_mg_rep0)]
);

# ╔═╡ c78f246c-f20d-4f6f-afa8-c12a316f488d
ΔR_sam_avg_seed = nanmean(ΔR_sam[:, seed_seqs]; dim=2)

# ╔═╡ 071338b6-2bb0-449d-b562-e23155436612
ΔR_sam_std_seed = nanstd(ΔR_sam[:, seed_seqs]; dim=2);

# ╔═╡ bd37015e-a662-4b7f-9694-8558306cd0d9
ΔR_sam_avg_full = nanmean(ΔR_sam[:, full_seqs]; dim=2)

# ╔═╡ 484ed85f-fe50-44e0-8e5d-b79c1b22317b
ΔR_sam_std_full = nanstd(ΔR_sam[:, full_seqs]; dim=2);

# ╔═╡ 301a213f-fd4b-4570-bfbb-d0b12d4f76e0
ΔR_sam_avg_rbmlo = nanmean(ΔR_sam[:, _rbmlo]; dim=2)

# ╔═╡ be52f566-e2b0-4cce-bfea-3b7e40bf077a
ΔR_sam_std_rbmlo = nanstd(ΔR_sam[:, _rbmlo]; dim=2);

# ╔═╡ 3dfad0ac-309b-49e7-9836-cd256fe59005
ΔR_sam_avg_rbmhi = nanmean(ΔR_sam[:, _rbmhi]; dim=2)

# ╔═╡ 91758ffb-1c13-4c8d-8ae5-fd66c0225fee
ΔR_sam_std_rbmhi = nanstd(ΔR_sam[:, _rbmhi]; dim=2);

# ╔═╡ 1d573ee1-eae1-4a96-b888-5b709101709e
ΔR_sam_avg_inf = nanmean(ΔR_sam[:, inf_seqs]; dim=2)

# ╔═╡ 91a97d59-c9fa-4141-a9e1-9cf9f0663405
ΔR_sam_std_inf = nanstd(ΔR_sam[:, inf_seqs]; dim=2);

# ╔═╡ 08a76f04-f5bd-4850-82e4-4c4fcea9dc9a
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, seed_seqs, conds_sam_rep0];

# ╔═╡ b0c3959c-fb33-4caf-bd16-96885a0cd302
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, seed_seqs, conds_sam_rep0];

# ╔═╡ fdb72e46-1410-488a-89b7-815c288277af
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ a9b44f72-eec6-4a23-a308-a4fafeac34e8
_sites = 3:107

# ╔═╡ 0276201f-b846-4573-b28d-e7e756ba0a50
md"""
# Figures
"""

# ╔═╡ 5ef6ad69-e574-4643-912b-db489114995e
let fig = Makie.Figure()
	ax = Makie.Axis(
		fig[1,1], width=300, height=300, xlabel="SHAPE reactivity", ylabel="frequency", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true
	)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), label="base paired", normalization=:pdf, bins=-2:0.05:6, linewidth=3, color=:teal)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), label="not paired", normalization=:pdf, bins=-2:0.05:6, linewidth=3, color=:orange)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	#Makie.axislegend(ax, framevisible=false, patchlabelgap=3, position=(-0.02, 1))
	Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r)
	
	_dummy_ax = Makie.Axis(fig[1,2], width=20, xgridvisible=false, ygridvisible=false)
	Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
	Makie.hidexdecorations!(_dummy_ax)
	Makie.hideydecorations!(_dummy_ax)
	
	ax = Makie.Axis(fig[1,3], width=300, height=300, xlabel="SHAPE reactivity", ylabel="frequency", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), label="b.p.", normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), label="n.p.", normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[pks, nat_seqs, conds_mg_rep0])), label="p.k.", normalization=:pdf, bins=-2:0.1:6, linewidth=3, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	ax = Makie.Axis(fig[1,4], width=300, height=300, xlabel="SHAPE reactivity", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[pks, nat_seqs, conds_sam_rep0])), label="pseudoknot", normalization=:pdf, bins=-2:0.1:6, linewidth=3, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	_xs = _sites
	
	ax = Makie.Axis(fig[2,:], width=900, height=150, xticks=5:10:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_seed - ΔR_sam_std_seed/2)[_xs], (ΔR_sam_avg_seed + ΔR_sam_std_seed/2)[_xs], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _xs, ΔR_sam_avg_seed[_xs], linewidth=1, color=:gray)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_seed[_xs], markersize=5, color=:black, label="Natural")
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_rbmlo - ΔR_sam_std_rbmlo/2)[_xs], (ΔR_sam_avg_rbmlo + ΔR_sam_std_rbmlo/2)[_xs], markersize=5, color=(:blue, 0.25))
	Makie.lines!(ax, _xs, ΔR_sam_avg_rbmlo[_xs], linewidth=1, color=:blue)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_rbmlo[_xs], markersize=5, color=:blue, label="RBM (RBMscore>300)")
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.xlims!(1, 108)
	Makie.hidespines!(ax, :t, :r, :b)
	Makie.hidexdecorations!(ax)
	
	ax = Makie.Axis(fig[3,:], width=900, height=150, xticks=5:10:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_seed - ΔR_sam_std_seed/2)[_xs], (ΔR_sam_avg_seed + ΔR_sam_std_seed/2)[_xs], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _xs, ΔR_sam_avg_seed[_xs], linewidth=1, color=:gray)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_seed[_xs], markersize=5, color=:black, label="Natural")
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_inf - ΔR_sam_std_inf/2)[_xs], (ΔR_sam_avg_inf + ΔR_sam_std_inf/2)[_xs], markersize=5, color=(:red, 0.25))
	Makie.lines!(ax, _xs, ΔR_sam_avg_inf[_xs], linewidth=1, color=:red)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_inf[_xs], markersize=5, color=:red, label="rCM")
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax, :t, :r)
	Makie.xlims!(1, 108)
	
	# Makie.Label(fig[1,1][1,1,Makie.TopLeft()], "A)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,2][1,1,Makie.TopLeft()], "B)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,3][1,1,Makie.TopLeft()], "C)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[2,:][1,1,Makie.TopLeft()], "D)", font=:bold, padding=(0,0,0,0))
	# Makie.Label(fig[3,:][1,1,Makie.TopLeft()], "E)", font=:bold, padding=(0,0,0,0))
	
	Makie.resize_to_layout!(fig)
	#Makie.save("/workspaces/SamApp.jl/notebooks/2024-03-14 New paper figures/Figures/SHAPE reactivities.pdf", fig)
	fig
end

# ╔═╡ 94019df5-2de3-4d18-aa73-4d6b1537b835
md"""
# Significance of correlation difference by bootstrapping

Following https://stats.stackexchange.com/questions/278751/how-do-i-determine-whether-two-correlations-are-significantly-different.
"""

# ╔═╡ 1eced2e8-6c9c-4482-8c5d-6f0212ca1bde
cor(ΔR_sam_avg_seed[_sites], ΔR_sam_avg_inf[_sites])

# ╔═╡ 946c2fb9-ed2d-40bd-8942-38c21d6d38e5
cor(ΔR_sam_avg_seed[_sites], ΔR_sam_avg_rbmlo[_sites])

# ╔═╡ 1c183be2-7c10-4ddb-b9df-5b8216de3768
let fig = Makie.Figure()
	N = 50000
	resampled_cor_inf = [cor(nanmean(ΔR_sam[_sites, rand(seed_seqs, length(seed_seqs))]; dim=2), nanmean(ΔR_sam[_sites, rand(inf_seqs, length(inf_seqs))]; dim=2)) for _ = 1:N]
	resampled_cor_rbm = [cor(nanmean(ΔR_sam[_sites, rand(seed_seqs, length(seed_seqs))]; dim=2), nanmean(ΔR_sam[_sites, rand(_rbmlo, length(_rbmlo))]; dim=2)) for _ = 1:N]

	bins = -0.2:0.01:1
	
	ax = Makie.Axis(fig[1,1]; width=300, height=300, xgridvisible=false, ygridvisible=false, xlabel="Bootstrapped correlation")
	Makie.hist!(ax, resampled_cor_inf; normalization=:pdf, color=(:red, 0.5), linewidth=2, bins)
	Makie.hist!(ax, resampled_cor_rbm; normalization=:pdf, color=(:blue, 0.5), linewidth=2, bins)
	Makie.stephist!(ax, resampled_cor_inf; normalization=:pdf, label="CM", color=:red, linewidth=2, bins)
	Makie.stephist!(ax, resampled_cor_rbm; normalization=:pdf, label="RBM", color=:blue, linewidth=2, bins)
	Makie.xlims!(ax, -0.2, 1)
	Makie.ylims!(ax, 0, 12)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	ax = Makie.Axis(fig[1,2]; width=300, height=300, xgridvisible=false, ygridvisible=false, xlabel="Bootstrapped correlation difference (RBM - CM)")
	Makie.hist!(ax, resampled_cor_rbm - resampled_cor_inf; normalization=:pdf, color=(:black, 0.5), linewidth=2, bins)
	Makie.stephist!(ax, resampled_cor_rbm - resampled_cor_inf; normalization=:pdf, color=:black, linewidth=2, bins)
	Makie.xlims!(ax, -0.2, 1)
	Makie.ylims!(ax, 0, 12)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ e9a343df-5c1c-4650-ae6d-482ed24af863
let N = 1000000
	resampled_cor_inf = [cor(nanmean(ΔR_sam[_sites, rand(seed_seqs, length(seed_seqs))]; dim=2), nanmean(ΔR_sam[_sites, rand(inf_seqs, length(inf_seqs))]; dim=2)) for _ = 1:N]
	resampled_cor_rbm = [cor(nanmean(ΔR_sam[_sites, rand(seed_seqs, length(seed_seqs))]; dim=2), nanmean(ΔR_sam[_sites, rand(_rbmlo, length(_rbmlo))]; dim=2)) for _ = 1:N]
	mean(resampled_cor_rbm .< resampled_cor_inf)
end

# ╔═╡ 65dd700e-8bf0-4c2e-83f3-a2ef6b592255
1/1000000

# ╔═╡ d492a8ee-ae2a-40a0-bf15-6c466097cd78
md"""
# Figure with larger std band
"""

# ╔═╡ a67aee24-3efb-4d69-bf42-d4103b2487f2
let fig = Makie.Figure()
	ax = Makie.Axis(
		fig[1,1], width=300, height=300, xlabel="SHAPE reactivity", ylabel="frequency", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true
	)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), label="base paired", normalization=:pdf, bins=-2:0.05:6, linewidth=3, color=:teal)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), label="not paired", normalization=:pdf, bins=-2:0.05:6, linewidth=3, color=:orange)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	#Makie.axislegend(ax, framevisible=false, patchlabelgap=3, position=(-0.02, 1))
	Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r)
	
	_dummy_ax = Makie.Axis(fig[1,2], width=20, xgridvisible=false, ygridvisible=false)
	Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
	Makie.hidexdecorations!(_dummy_ax)
	Makie.hideydecorations!(_dummy_ax)
	
	ax = Makie.Axis(fig[1,3], width=300, height=300, xlabel="SHAPE reactivity", ylabel="frequency", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), label="b.p.", normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), label="n.p.", normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[pks, nat_seqs, conds_mg_rep0])), label="p.k.", normalization=:pdf, bins=-2:0.1:6, linewidth=3, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	ax = Makie.Axis(fig[1,4], width=300, height=300, xlabel="SHAPE reactivity", xgridvisible=false, ygridvisible=false, xticks=-2:2:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0])), normalization=:pdf, bins=-2:0.05:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, vec(shape_data_rep0.shape_reactivities[pks, nat_seqs, conds_sam_rep0])), label="pseudoknot", normalization=:pdf, bins=-2:0.1:6, linewidth=3, color=:black)
	Makie.xlims!(-2.2, 6)
	Makie.ylims!(-0.07, 2)
	Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	_xs = 3:107
	
	ax = Makie.Axis(fig[2,:], width=900, height=150, xticks=5:10:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_seed - ΔR_sam_std_seed)[_xs], (ΔR_sam_avg_seed + ΔR_sam_std_seed)[_xs], markersize=5, color=(:gray, 0.2))
	Makie.lines!(ax, _xs, ΔR_sam_avg_seed[_xs], linewidth=1, color=:gray)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_seed[_xs], markersize=5, color=:black, label="Natural")
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_rbmlo - ΔR_sam_std_rbmlo)[_xs], (ΔR_sam_avg_rbmlo + ΔR_sam_std_rbmlo)[_xs], markersize=5, color=(:blue, 0.2))
	Makie.lines!(ax, _xs, ΔR_sam_avg_rbmlo[_xs], linewidth=1, color=:blue)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_rbmlo[_xs], markersize=5, color=:blue, label="RBM (RBMscore>300)")
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.xlims!(1, 108)
	Makie.hidespines!(ax, :t, :r, :b)
	Makie.hidexdecorations!(ax)
	
	ax = Makie.Axis(fig[3,:], width=900, height=150, xticks=5:10:108, yticks=-2:1:1, xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_seed - ΔR_sam_std_seed)[_xs], (ΔR_sam_avg_seed + ΔR_sam_std_seed)[_xs], markersize=5, color=(:gray, 0.2))
	Makie.lines!(ax, _xs, ΔR_sam_avg_seed[_xs], linewidth=1, color=:gray)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_seed[_xs], markersize=5, color=:black, label="Natural")
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_inf - ΔR_sam_std_inf)[_xs], (ΔR_sam_avg_inf + ΔR_sam_std_inf)[_xs], markersize=5, color=(:red, 0.2))
	Makie.lines!(ax, _xs, ΔR_sam_avg_inf[_xs], linewidth=1, color=:red)
	Makie.scatter!(ax, _xs, ΔR_sam_avg_inf[_xs], markersize=5, color=:red, label="rCM")
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax, :t, :r)
	Makie.xlims!(1, 108)
	
	# Makie.Label(fig[1,1][1,1,Makie.TopLeft()], "A)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,2][1,1,Makie.TopLeft()], "B)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,3][1,1,Makie.TopLeft()], "C)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[2,:][1,1,Makie.TopLeft()], "D)", font=:bold, padding=(0,0,0,0))
	# Makie.Label(fig[3,:][1,1,Makie.TopLeft()], "E)", font=:bold, padding=(0,0,0,0))
	
	Makie.resize_to_layout!(fig)
	#Makie.save("/workspaces/SamApp.jl/notebooks/2024-03-14 New paper figures/Figures/SHAPE reactivities.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═d443b9aa-7449-4447-983c-18263e61be12
# ╠═91f4edd8-290d-4270-83c1-f7c6281e9f68
# ╠═59e3b0b4-c955-4985-bc79-bbc9220539a6
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
# ╠═e743ce81-29e8-4bac-b22e-79e123db3326
# ╠═9299b1d4-a573-4808-9e26-396c2e04b8cc
# ╠═530dbe76-e3a2-4658-ab86-daf6588cf372
# ╠═0c2a4cd9-a42a-42fc-b85b-46e5987a77cc
# ╠═8da6f416-1a00-4cc8-a9c0-cdd76ec364b0
# ╠═d705aa70-ad0e-47aa-a31b-714984d3b65c
# ╠═18452b0f-3c62-492c-a41b-6a9043ce86b1
# ╠═4206fb3f-99c4-4e4b-bbc9-a97891a2f321
# ╠═dcf0701e-3fcc-4199-8285-803305ee38b2
# ╠═d1387ef6-6497-4f0b-873c-4de5f8aa0715
# ╠═0297ebf2-0822-4ae0-b4e1-73b8ca1c4452
# ╠═c011571f-74e7-4bae-be5e-918ad678397c
# ╠═e39876f7-5399-4e4a-a1d5-101355741f28
# ╠═d527106f-7e5b-49a5-930a-fd847a782542
# ╠═3ac302e6-b7c6-474d-af2b-c12420ebd0be
# ╠═64cb3148-1d0a-4f7d-8859-ded1e2fc5283
# ╠═f9b70f45-5ce5-460a-b7bd-46a9cca2c8f3
# ╠═4e123003-aede-4e18-956f-7f87a85d5083
# ╠═6bdf0eb7-511a-4900-bc51-66e960cd391c
# ╠═6994d41d-0fdf-4dd6-bc8e-122d745ad34c
# ╠═44042f22-0ae9-44a4-ad22-7d316e976b92
# ╠═29896cb1-147d-4e15-a815-a2a55e7ce483
# ╠═7a266977-5fea-4a23-84d0-e967650cc91e
# ╠═4762b1ab-4dc9-4735-a9b3-6ab26a6c19d7
# ╠═143dbfef-3376-4892-998a-a769acb78a05
# ╠═0fd1db6a-c51e-4e18-85aa-3064103b345f
# ╠═880f1579-fae9-4cb4-9fd5-b533203fa0e0
# ╠═14da7dad-da45-457f-ac4a-3a59dabf0898
# ╠═d71ccb35-7881-4571-8c5c-bac572630c37
# ╠═c78f246c-f20d-4f6f-afa8-c12a316f488d
# ╠═071338b6-2bb0-449d-b562-e23155436612
# ╠═bd37015e-a662-4b7f-9694-8558306cd0d9
# ╠═484ed85f-fe50-44e0-8e5d-b79c1b22317b
# ╠═301a213f-fd4b-4570-bfbb-d0b12d4f76e0
# ╠═be52f566-e2b0-4cce-bfea-3b7e40bf077a
# ╠═3dfad0ac-309b-49e7-9836-cd256fe59005
# ╠═91758ffb-1c13-4c8d-8ae5-fd66c0225fee
# ╠═1d573ee1-eae1-4a96-b888-5b709101709e
# ╠═91a97d59-c9fa-4141-a9e1-9cf9f0663405
# ╠═08a76f04-f5bd-4850-82e4-4c4fcea9dc9a
# ╠═b0c3959c-fb33-4caf-bd16-96885a0cd302
# ╠═fdb72e46-1410-488a-89b7-815c288277af
# ╠═a9b44f72-eec6-4a23-a308-a4fafeac34e8
# ╠═0276201f-b846-4573-b28d-e7e756ba0a50
# ╠═5ef6ad69-e574-4643-912b-db489114995e
# ╠═94019df5-2de3-4d18-aa73-4d6b1537b835
# ╠═1eced2e8-6c9c-4482-8c5d-6f0212ca1bde
# ╠═946c2fb9-ed2d-40bd-8942-38c21d6d38e5
# ╠═1c183be2-7c10-4ddb-b9df-5b8216de3768
# ╠═e9a343df-5c1c-4650-ae6d-482ed24af863
# ╠═65dd700e-8bf0-4c2e-83f3-a2ef6b592255
# ╠═d492a8ee-ae2a-40a0-bf15-6c466097cd78
# ╠═a67aee24-3efb-4d69-bf42-d4103b2487f2

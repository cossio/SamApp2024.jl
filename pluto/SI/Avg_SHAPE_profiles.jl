### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ c66b138e-05dc-4f93-ba69-1d7a2e8e2bc5
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ d0862cc4-369b-46f9-a612-c93832c73f58
using BioSequences: LongRNA

# ╔═╡ 5834c057-3087-4850-bcf9-c83c3754622d
using Distributions: Gamma

# ╔═╡ b4a50bf4-dc61-4279-860e-149555c1ef1b
using Makie: @L_str

# ╔═╡ 98157aeb-819d-41ce-b453-1a7535a13fd7
using NaNStatistics: nanmean

# ╔═╡ 8d5835b7-950c-4c16-82e9-53376e9df802
using NaNStatistics: nanstd

# ╔═╡ 9640f397-ce77-4ce2-b963-5e7297b00019
using NaNStatistics: nancor

# ╔═╡ c111079d-2fe2-4abf-894b-e5a0567e0351
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ d0da5b04-9fd8-48f0-ab9d-829d2daacbc8
using Statistics: cor

# ╔═╡ 4ddd5cb8-93c9-4dcd-ad69-01926f2c56a7
using Statistics: mean

# ╔═╡ 3b512609-68e8-443a-a5f7-2126c56da6bc
using NaNStatistics: nansum

# ╔═╡ 68115a4a-e980-4194-aafa-069616330d43
md"""
# Imports
"""

# ╔═╡ 986dc1d3-023f-4389-9a42-c91264a15c24
import PlutoUI

# ╔═╡ 61f57084-e80c-4ffc-9b8f-9a128984151e
import CairoMakie

# ╔═╡ 279ca90f-35d4-4702-8206-695128358782
import Makie

# ╔═╡ d9e94bf5-82f3-4a99-a304-7b6e432698c2
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 57d5fca7-ff24-40a6-86b6-3c86a8be338f
import SamApp2024

# ╔═╡ 6c3e0faa-2a8f-467f-b9dd-e78174125bf0
import StatsBase

# ╔═╡ cbc38a01-15db-494c-9838-bbb9b60d1dba
import Infernal

# ╔═╡ 16e81b9f-a494-4a84-b861-fcdb4442f67a
import Rfam

# ╔═╡ b8d19662-fa27-4610-8091-c85ccac45be0
import FASTX

# ╔═╡ 7cf5352f-284b-442c-bd3a-bbd6877034c9
PlutoUI.TableOfContents()

# ╔═╡ ac1376ac-2fe3-42cc-8f2e-3fd8766a09ea
md"""
# Load data
"""

# ╔═╡ cc53f278-b947-4f29-898d-93e41d583a0f
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ c48bd796-3ba5-4a87-9c9a-83559b99718d
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 3bf379e7-e0ed-4b71-8826-47783635ef40
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 3301c847-8709-4aed-b62f-d87a441e89e6
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 54436436-1cad-4708-92ff-83d9498a947b
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 52a73bc3-943c-4af3-9053-401aca0d5221
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 5ab91f1e-4deb-4727-9588-3c417b559211
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 35d18d70-30a6-4b27-9977-a622aacd5cc7
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ e2e76c04-06d4-4260-8a4b-7a6fca3667f3
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 33c9124c-c057-4c6e-b110-7a2b3ec37155
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 3d6ad8f8-9509-465a-9af8-dcb420055724
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 2e09ec7f-778b-4b58-bbef-f383b3957106
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ e03eb8a3-e1c7-4bf7-a775-da70d8a4a99c
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ 064da099-08b5-4fb1-853c-1ad7cf0902bc
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ ce97d10c-57be-4e3c-b0fa-7dc0200d85ae
_rbmlo = rbm_seqs ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .< -300));

# ╔═╡ f9da2ba1-5e14-4b66-a8d5-35dd293de6b5
_rbmhi = rbm_seqs ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .> -300));

# ╔═╡ 662dad5e-309e-4a51-8aa9-7b8e9a68b933
ΔR_sam = (
    nanmean(shape_data_rep0.shape_reactivities[:, :, conds_sam_rep0]; dim=3) .- 
    shape_data_rep0.shape_reactivities[:, :, only(conds_mg_rep0)]
);

# ╔═╡ e9e4caf5-53a4-45c0-b60a-8f30c86a1f0a
ΔR_mg = (
    shape_data_rep0.shape_reactivities[:, :, only(conds_mg_rep0)] .-
    shape_data_rep0.shape_reactivities[:, :, only(conds_30C_rep0)]
);

# ╔═╡ ea7d249a-2f97-4257-a384-23bcf898e1d1
ΔR_sam_avg_seed = nanmean(ΔR_sam[:, seed_seqs]; dim=2)

# ╔═╡ d46610f4-5880-4ed9-a11b-f46af4a350ec
ΔR_sam_std_seed = nanstd(ΔR_sam[:, seed_seqs]; dim=2)

# ╔═╡ c6d25a44-acae-4ac1-85d8-9c771ade2ad2
ΔR_mg_avg_seed = nanmean(ΔR_mg[:, seed_seqs]; dim=2)

# ╔═╡ 521c0736-50b6-402e-bd66-993f573f5bd2
ΔR_mg_std_seed = nanstd(ΔR_mg[:, seed_seqs]; dim=2)

# ╔═╡ 2bf75879-eb41-4811-af27-2b75d5a03b94
ΔR_sam_avg_full = nanmean(ΔR_sam[:, full_seqs]; dim=2)

# ╔═╡ 50630028-9ef9-4a99-9960-239e78905ea1
ΔR_sam_std_full = nanstd(ΔR_sam[:, full_seqs]; dim=2);

# ╔═╡ 2422f14c-0b02-4059-8fd2-72bd6110ac9c
ΔR_sam_avg_rbmlo = nanmean(ΔR_sam[:, _rbmlo]; dim=2)

# ╔═╡ ab59563c-24aa-4343-aabd-bc83af634d47
ΔR_sam_std_rbmlo = nanstd(ΔR_sam[:, _rbmlo]; dim=2)

# ╔═╡ 3a341520-7ec3-47c5-bff0-69dcde6befa3
ΔR_mg_avg_rbmlo = nanmean(ΔR_mg[:, _rbmlo]; dim=2)

# ╔═╡ a0759269-97b3-46e8-a116-08c92044094a
ΔR_mg_std_rbmlo = nanstd(ΔR_mg[:, _rbmlo]; dim=2)

# ╔═╡ 85529c26-e375-4fd3-a0ac-2c29214d2932
ΔR_sam_avg_rbmhi = nanmean(ΔR_sam[:, _rbmhi]; dim=2)

# ╔═╡ 28e44344-004f-4d49-9dc4-736ca1232f0b
ΔR_sam_std_rbmhi = nanstd(ΔR_sam[:, _rbmhi]; dim=2)

# ╔═╡ 7ae461b8-e707-42cd-b5fa-45416de4959a
ΔR_mg_avg_rbmhi = nanmean(ΔR_mg[:, _rbmhi]; dim=2)

# ╔═╡ 5f55fb6d-1d25-4d08-a221-f1e965c32e18
ΔR_mg_std_rbmhi = nanstd(ΔR_mg[:, _rbmhi]; dim=2)

# ╔═╡ 2e1dcf63-3045-4b3c-aec6-b0226a0dd509
ΔR_sam_avg_inf = nanmean(ΔR_sam[:, inf_seqs]; dim=2)

# ╔═╡ 4338e4c9-b19a-45c4-b183-7c0bdb722849
ΔR_sam_std_inf = nanstd(ΔR_sam[:, inf_seqs]; dim=2)

# ╔═╡ 3f5ec14c-d887-41f4-8729-ff8eb7f9b66d
ΔR_mg_avg_inf = nanmean(ΔR_mg[:, inf_seqs]; dim=2)

# ╔═╡ 66300085-651d-4a5a-86c9-5366e4d2388c
ΔR_mg_std_inf = nanstd(ΔR_mg[:, inf_seqs]; dim=2)

# ╔═╡ 1e1f186f-0d28-4dcf-9fee-adbafe9317cb
_sites = 3:107 # do not plot sites at the edges because they contain NAN / Inf

# ╔═╡ ccef9ed1-5921-4db7-b31f-6c4f7158810b
md"""
# Figure: Response to Mg
"""

# ╔═╡ 8aefd350-7e20-494f-9422-d1351d729dff
let fig = Makie.Figure()

	ax = Makie.Axis(fig[1,1]; width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_seed - ΔR_mg_std_seed/2)[_sites], (ΔR_mg_avg_seed + ΔR_mg_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_rbmlo - ΔR_mg_std_rbmlo/2)[_sites], (ΔR_mg_avg_rbmlo + ΔR_mg_std_rbmlo/2)[_sites], markersize=5, color=(:blue, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_rbmlo[_sites], linewidth=1, color=:blue)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_rbmlo[_sites], markersize=5, color=:blue)
	
	
	ax = Makie.Axis(fig[2,1], width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_seed - ΔR_mg_std_seed/2)[_sites], (ΔR_mg_avg_seed + ΔR_mg_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_rbmhi - ΔR_mg_std_rbmhi/2)[_sites], (ΔR_mg_avg_rbmhi + ΔR_mg_std_rbmhi/2)[_sites], markersize=5, color=(:red, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_rbmhi[_sites], linewidth=1, color=:red)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_rbmhi[_sites], markersize=5, color=:red)
	
	
	ax = Makie.Axis(fig[3,1], width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_seed - ΔR_mg_std_seed/2)[_sites], (ΔR_mg_avg_seed + ΔR_mg_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_mg_avg_inf - ΔR_mg_std_inf/2)[_sites], (ΔR_mg_avg_inf + ΔR_mg_std_inf/2)[_sites], markersize=5, color=(:red, 0.25))
	Makie.lines!(ax, _sites, ΔR_mg_avg_inf[_sites], linewidth=1, color=:red)
	Makie.scatter!(ax, _sites, ΔR_mg_avg_inf[_sites], markersize=5, color=:red)

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "B)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[3,1][1, 1, Makie.TopLeft()], "C)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	Makie.save("/DATA/cossio/SAM/2024/SamApp2024.jl/pluto/SI/Figures/2024-07-30 Avg SHAPE profiles - Mg.pdf", fig)
	fig
end

# ╔═╡ 95ecd643-da11-4b50-b7dc-00dfff1c54ce
nancor(ΔR_mg_avg_seed[_sites], ΔR_mg_avg_rbmlo[_sites])

# ╔═╡ b802efbc-1246-4437-a743-8b9e359df09f
nancor(ΔR_mg_avg_seed[_sites], ΔR_mg_avg_inf[_sites])

# ╔═╡ bf4201e8-ea61-49dd-9746-06f6d283e931
nancor(ΔR_mg_avg_seed[_sites], ΔR_mg_avg_rbmhi[_sites])

# ╔═╡ 462694bb-2ed1-4b0b-9a1b-c3f5e046ea32
md"""
# Figure: Response to SAM
"""

# ╔═╡ 273aa7d7-dea1-409d-9e1a-e845b3cadaec
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_seed - ΔR_sam_std_seed/2)[_sites], (ΔR_sam_avg_seed + ΔR_sam_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_rbmlo - ΔR_sam_std_rbmlo/2)[_sites], (ΔR_sam_avg_rbmlo + ΔR_sam_std_rbmlo/2)[_sites], markersize=5, color=(:blue, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_rbmlo[_sites], linewidth=1, color=:blue)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_rbmlo[_sites], markersize=5, color=:blue)
	
	
	ax = Makie.Axis(fig[2,1], width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_seed - ΔR_sam_std_seed/2)[_sites], (ΔR_sam_avg_seed + ΔR_sam_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_rbmhi - ΔR_sam_std_rbmhi/2)[_sites], (ΔR_sam_avg_rbmhi + ΔR_sam_std_rbmhi/2)[_sites], markersize=5, color=(:red, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_rbmhi[_sites], linewidth=1, color=:red)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_rbmhi[_sites], markersize=5, color=:red)
	
	
	ax = Makie.Axis(fig[3,1], width=700, height=100, xticks=1:5:108, yticks=-2:1:1, xlabel="site", ylabel="diff. reactivity", xgridvisible=false, ygridvisible=false)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_seed - ΔR_sam_std_seed/2)[_sites], (ΔR_sam_avg_seed + ΔR_sam_std_seed/2)[_sites], markersize=5, color=(:gray, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_seed[_sites], linewidth=1, color=:gray)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_seed[_sites], markersize=5, color=:black)
	
	Makie.band!(ax, _sites, (ΔR_sam_avg_inf - ΔR_sam_std_inf/2)[_sites], (ΔR_sam_avg_inf + ΔR_sam_std_inf/2)[_sites], markersize=5, color=(:red, 0.25))
	Makie.lines!(ax, _sites, ΔR_sam_avg_inf[_sites], linewidth=1, color=:red)
	Makie.scatter!(ax, _sites, ΔR_sam_avg_inf[_sites], markersize=5, color=:red)

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "B)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[3,1][1, 1, Makie.TopLeft()], "C)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	Makie.save("/DATA/cossio/SAM/2024/SamApp2024.jl/pluto/SI/Figures/2024-07-30 Avg SHAPE profiles - SAM.pdf", fig)
	fig
end

# ╔═╡ 211a2d6f-bb96-4a9d-83c8-3e653d85ba47
nancor(ΔR_sam_avg_seed[_sites], ΔR_sam_avg_rbmlo[_sites])

# ╔═╡ d1e686b6-2e9b-4c75-9b6a-51d97e65c23f
nancor(ΔR_sam_avg_seed[_sites], ΔR_sam_avg_inf[_sites])

# ╔═╡ 7b236bcd-6127-4ff5-bb32-d18bd0b197b1
nancor(ΔR_sam_avg_seed[_sites], ΔR_sam_avg_rbmhi[_sites])

# ╔═╡ Cell order:
# ╠═68115a4a-e980-4194-aafa-069616330d43
# ╠═c66b138e-05dc-4f93-ba69-1d7a2e8e2bc5
# ╠═986dc1d3-023f-4389-9a42-c91264a15c24
# ╠═61f57084-e80c-4ffc-9b8f-9a128984151e
# ╠═279ca90f-35d4-4702-8206-695128358782
# ╠═d9e94bf5-82f3-4a99-a304-7b6e432698c2
# ╠═57d5fca7-ff24-40a6-86b6-3c86a8be338f
# ╠═6c3e0faa-2a8f-467f-b9dd-e78174125bf0
# ╠═cbc38a01-15db-494c-9838-bbb9b60d1dba
# ╠═16e81b9f-a494-4a84-b861-fcdb4442f67a
# ╠═b8d19662-fa27-4610-8091-c85ccac45be0
# ╠═d0862cc4-369b-46f9-a612-c93832c73f58
# ╠═5834c057-3087-4850-bcf9-c83c3754622d
# ╠═b4a50bf4-dc61-4279-860e-149555c1ef1b
# ╠═98157aeb-819d-41ce-b453-1a7535a13fd7
# ╠═8d5835b7-950c-4c16-82e9-53376e9df802
# ╠═9640f397-ce77-4ce2-b963-5e7297b00019
# ╠═c111079d-2fe2-4abf-894b-e5a0567e0351
# ╠═d0da5b04-9fd8-48f0-ab9d-829d2daacbc8
# ╠═4ddd5cb8-93c9-4dcd-ad69-01926f2c56a7
# ╠═3b512609-68e8-443a-a5f7-2126c56da6bc
# ╠═7cf5352f-284b-442c-bd3a-bbd6877034c9
# ╠═ac1376ac-2fe3-42cc-8f2e-3fd8766a09ea
# ╠═cc53f278-b947-4f29-898d-93e41d583a0f
# ╠═c48bd796-3ba5-4a87-9c9a-83559b99718d
# ╠═3bf379e7-e0ed-4b71-8826-47783635ef40
# ╠═3301c847-8709-4aed-b62f-d87a441e89e6
# ╠═54436436-1cad-4708-92ff-83d9498a947b
# ╠═52a73bc3-943c-4af3-9053-401aca0d5221
# ╠═5ab91f1e-4deb-4727-9588-3c417b559211
# ╠═35d18d70-30a6-4b27-9977-a622aacd5cc7
# ╠═e2e76c04-06d4-4260-8a4b-7a6fca3667f3
# ╠═33c9124c-c057-4c6e-b110-7a2b3ec37155
# ╠═3d6ad8f8-9509-465a-9af8-dcb420055724
# ╠═2e09ec7f-778b-4b58-bbef-f383b3957106
# ╠═e03eb8a3-e1c7-4bf7-a775-da70d8a4a99c
# ╠═064da099-08b5-4fb1-853c-1ad7cf0902bc
# ╠═ce97d10c-57be-4e3c-b0fa-7dc0200d85ae
# ╠═f9da2ba1-5e14-4b66-a8d5-35dd293de6b5
# ╠═662dad5e-309e-4a51-8aa9-7b8e9a68b933
# ╠═e9e4caf5-53a4-45c0-b60a-8f30c86a1f0a
# ╠═ea7d249a-2f97-4257-a384-23bcf898e1d1
# ╠═d46610f4-5880-4ed9-a11b-f46af4a350ec
# ╠═c6d25a44-acae-4ac1-85d8-9c771ade2ad2
# ╠═521c0736-50b6-402e-bd66-993f573f5bd2
# ╠═2bf75879-eb41-4811-af27-2b75d5a03b94
# ╠═50630028-9ef9-4a99-9960-239e78905ea1
# ╠═2422f14c-0b02-4059-8fd2-72bd6110ac9c
# ╠═ab59563c-24aa-4343-aabd-bc83af634d47
# ╠═3a341520-7ec3-47c5-bff0-69dcde6befa3
# ╠═a0759269-97b3-46e8-a116-08c92044094a
# ╠═85529c26-e375-4fd3-a0ac-2c29214d2932
# ╠═28e44344-004f-4d49-9dc4-736ca1232f0b
# ╠═7ae461b8-e707-42cd-b5fa-45416de4959a
# ╠═5f55fb6d-1d25-4d08-a221-f1e965c32e18
# ╠═2e1dcf63-3045-4b3c-aec6-b0226a0dd509
# ╠═4338e4c9-b19a-45c4-b183-7c0bdb722849
# ╠═3f5ec14c-d887-41f4-8729-ff8eb7f9b66d
# ╠═66300085-651d-4a5a-86c9-5366e4d2388c
# ╠═1e1f186f-0d28-4dcf-9fee-adbafe9317cb
# ╠═ccef9ed1-5921-4db7-b31f-6c4f7158810b
# ╠═8aefd350-7e20-494f-9422-d1351d729dff
# ╠═95ecd643-da11-4b50-b7dc-00dfff1c54ce
# ╠═b802efbc-1246-4437-a743-8b9e359df09f
# ╠═bf4201e8-ea61-49dd-9746-06f6d283e931
# ╠═462694bb-2ed1-4b0b-9a1b-c3f5e046ea32
# ╠═273aa7d7-dea1-409d-9e1a-e845b3cadaec
# ╠═211a2d6f-bb96-4a9d-83c8-3e653d85ba47
# ╠═d1e686b6-2e9b-4c75-9b6a-51d97e65c23f
# ╠═7b236bcd-6127-4ff5-bb32-d18bd0b197b1

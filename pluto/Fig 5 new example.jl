### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8e34985c-368c-4a4f-99c0-5cc16d195974
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 4afbfcac-2e59-44d0-a850-5c49c7167720
using BioSequences: LongRNA

# ╔═╡ 75374067-8150-4721-b83f-79fae3ac12b2
using DataFrames: DataFrame

# ╔═╡ dcce26c2-893a-4963-9852-2a7880167da2
using Makie: @L_str

# ╔═╡ e0bca8ed-50f4-482d-896a-edb2eeb154fe
using NaNStatistics: nansum

# ╔═╡ ba434302-f435-4dd8-82c0-c5b67028c848
using Random: bitrand

# ╔═╡ 1fde47a6-3818-4c28-8934-f6a4e1c880f8
using Statistics: cor

# ╔═╡ 0bddc5ea-29de-4f27-bb14-e1e889a18993
using Statistics: mean

# ╔═╡ ee13efeb-8568-4583-97fe-65ede489e356
using StatsBase: countmap

# ╔═╡ 1f9e4648-3b0f-4ce8-acc1-8aa4a51e7242
md"# Imports"

# ╔═╡ 5e10a3a4-e5bb-441c-8ab5-6b1db64d9ebb
import Makie

# ╔═╡ c01cb31b-d0c9-4f42-82f6-43f60bb3d67c
import CairoMakie

# ╔═╡ d3ecc35c-ac26-4b86-a021-f5cf284de8cc
import FASTX

# ╔═╡ d7847339-7164-46f0-88e8-5b8193de577a
import Infernal

# ╔═╡ 4e597e09-0311-47d2-9f26-8ec6a1218128
import SamApp2024

# ╔═╡ 5002f57d-c32c-42d9-a350-417e85615943
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 20af2d04-5621-458c-b129-501f6a0da126
import Rfam

# ╔═╡ 398d2110-7faa-4878-bf53-f61d6e375f05
import PlutoUI

# ╔═╡ a94f1ab3-c6d6-427d-8836-d3f1e35f738d
import Unitful

# ╔═╡ d806a548-11a5-41a3-8591-1f52f46ee8c8
import ViennaRNA

# ╔═╡ fce1dfd3-89b0-4ff1-a0f9-c71109a683f8
import StatsBase

# ╔═╡ 6b943684-49d3-4eee-b5d6-aeb7fe4fef1f
import KernelDensity

# ╔═╡ b759c7bb-dc94-4dc5-bc83-6d9e29e93748
PlutoUI.TableOfContents()

# ╔═╡ 46eddbfb-70af-4bf4-b448-9f21067c6871
md"# Load data"

# ╔═╡ 8eb1f3ce-fbc9-461f-9e2a-ac8707702f9a
@show Rfam.get_rfam_directory() Rfam.get_rfam_version();

# ╔═╡ b9653255-2c4d-4758-9eb9-667937bde77d
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20240730_with_pdb();

# ╔═╡ 8fa8e888-3ce8-4726-83c3-7cafe3db9f14
# All merged data, for the reactivity profiles plots
shape_data_all_merged = SamApp2024.load_shapemapper_data_pierre_demux_20240801_with_pdb_repls_merged();

# ╔═╡ eecaf592-337e-4483-acf8-43cf7a68abd0
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ fc3cbbbc-d4f3-4789-894d-76921085e6fa
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 19279ab8-afdd-4ecb-bf46-41e7f0ce956e
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 8056ed26-0adb-4600-9a51-c0273e1ef322
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 5bcfc9a7-147f-4e3d-b1b9-89cb878fe9b0
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ cb10c61d-2066-4d01-a0d5-fe49e6a2e9b8
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 82dd4082-50cc-4eff-9407-cd7808d64688
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ f65cd678-a8c5-495d-b104-de5161cabb5e
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ f1a565af-6d83-427a-be5b-46792243f631
conds_SAM_all_merged = map(identity, indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_allrep", "SAMAP_1M7_1SAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 9a12a8c9-e86d-4ddd-a6c2-ffa419e06ea2
conds_Mg_all_merged = map(identity, indexin(["SAMAP_1M7_noSAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 7110c60b-1c56-4678-bad8-1cd1c6de4034
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 3e703ec9-f5ef-4fe2-ac60-9b6898249304
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ f9afb22d-0990-4616-8d03-7a0f75cf3bd0
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ b0163fae-c4f5-43dc-ad63-a8932d3ef3a4
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 646fa540-843f-4738-966e-f0dd69e7241b
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ d93a6822-20e0-417f-a93e-983ad647a750
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 46ff2334-74b3-45e0-8c1e-2d4e8a80de7b
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 26359219-b16b-4a73-bc48-7ce5deea4b3a
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ a4e7ef5c-7462-496c-86bc-9fdb3ad3611d
_thresh = log(5)

# ╔═╡ 57e2e479-a340-4345-9b52-8f11587344b5
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ ea96c786-cd96-4337-bb0a-20cbddde86d7
# structural motifs
struct_bands = [
    (; x0=0.5, xf=8.5, color="blue", alpha=0.1), # P1
    (; x0=100.5, xf=108.5, color="blue", alpha=0.1), # P1
    (; x0=11.5, xf=16.5, color="green", alpha=0.1), # P2
    (; x0=20.5, xf=23.5, color="green", alpha=0.1), # P2
    (; x0=28.5, xf=31.5, color="green", alpha=0.1), # P2
    (; x0=37.5, xf=42.5, color="green", alpha=0.1), # P2
    (; x0=42.5, xf=46.5, color="orange", alpha=0.1), # P3
    (; x0=47.5, xf=53.5, color="orange", alpha=0.1), # P3
    (; x0=60.5, xf=64.5, color="orange", alpha=0.1), # P3
    (; x0=66.5, xf=72.5, color="orange", alpha=0.1), # P3
    (; x0=80.5, xf=86.5, color="teal", alpha=0.1), # P4
    (; x0=91.5, xf=97.5, color="teal", alpha=0.1), # P4
    (; x0=24.5, xf=28.5, color="red", alpha=0.1), # Pk
    (; x0=76.5, xf=80.5, color="red", alpha=0.1), # Pk
];

# ╔═╡ 6fd2910a-bd6e-48c6-acec-9c9a4d3a74b9
md"# Secondary structure"

# ╔═╡ faa6a32d-6b02-4c85-aaaa-a38855ae119a
wuss = SamApp2024.rfam_ss("RF00162"; inserts=false)

# ╔═╡ 1aa04802-6a89-47cb-8994-aa9b40319fb5
ss = SamApp2024.clean_wuss(wuss)

# ╔═╡ 01c3a737-9faa-4b5a-8f5a-47ff3ad121f1
length(ss)

# ╔═╡ 0615574c-94c5-4679-bec3-7b952f58dd05
shape_data_all_merged.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), conds_SAM_all_merged[1]]

# ╔═╡ 9dd3cf23-2841-47df-94ba-861e911bb49a
shape_data_all_merged.aligned_sequences[only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10"))]

# ╔═╡ 09ed5c06-8ab8-4ea9-b0c5-98b256838071
length("GUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAC")

# ╔═╡ 1e5a5d35-3b84-4c42-af1c-310375ba449a
md"# Example 1: 4KQY"

# ╔═╡ 23e51bbf-8db0-470f-ba50-d6999480d922
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10"))
	@show shape_data_045.aptamer_ids[n_ex]
	
	width = 700
	height = 100
	xticks = 5:5:108

	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -1.5, 1)
	Makie.ylims!(ax_react, -0.5, 6)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example.pdf", fig)
	fig
end

# ╔═╡ ea8e20bc-fd8e-43cd-8ef5-6a653f28f53f
aptamers_df = SamApp2024.probed_aptamers_table_20221027()

# ╔═╡ 9ecff5e8-4b91-4338-aaaa-1200684d3765
pwd()

# ╔═╡ ff713d18-cb3a-4707-b694-a7acc6f0d1d2
pdb10_index_in_aptamers_df = only(findall(aptamers_df.name .== "SAMAP-PDB10"))

# ╔═╡ 9f491a8e-fee1-409c-b2bd-0fbc0be1270c
replace(aptamers_df.sequence[pdb10_index_in_aptamers_df], 'T' => 'U')

# ╔═╡ d12fb5cb-1315-4284-bb3b-57deb741cb87
length("GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAAC")

# ╔═╡ 4826b171-a3b5-4411-8357-bf2fc13c2f01
md"# Example 2"

# ╔═╡ 99def13a-3700-45c6-9127-64754f234957
conds_SAM_all_merged

# ╔═╡ e7ef002b-2567-434d-8e0f-b707ccd389d3
shape_data_all_merged.conditions

# ╔═╡ af800276-6e2f-4a14-b6f4-d565b0de6f45
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "APSAMN7"))
	@show shape_data_045.aptamer_ids[n_ex]
	
	width = 700
	height = 100
	xticks = 5:5:108

	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -1.5, 1)
	Makie.ylims!(ax_react, -0.5, 6)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example2.pdf", fig)
	fig
end

# ╔═╡ f88073e3-144f-459a-8910-7ad5a0c2340d
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "APSAMN7"))
	@show shape_data_045.aptamer_ids[n_ex]
	
	width = 700
	height = 100
	xticks = 5:5:108

	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[2]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -1.5, 1)
	Makie.ylims!(ax_react, -0.5, 6)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example2.pdf", fig)
	fig
end

# ╔═╡ 1b7a906e-f240-4c53-b7f0-719ce088908f
conds_SAM_all_merged

# ╔═╡ bbd8df1d-ce78-45cc-90a9-2eb1ba112f7b
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "APSAMN7"))
	@show shape_data_045.aptamer_ids[n_ex]
	
	width = 700
	height = 100
	xticks = 5:5:108

	_R_sam = shape_data_rep0.shape_reactivities[:, n_ex, conds_sam_rep0[1]]
	_R_mg = shape_data_rep0.shape_reactivities[:, n_ex, only(conds_mg_rep0)]
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -1.5, 1)
	Makie.ylims!(ax_react, -0.5, 6)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example2.pdf", fig)
	fig
end

# ╔═╡ cbd07af7-0688-467b-bf06-6e51f0cd9693
aptamers_df.description[only(findall(aptamers_df.name .== "APSAMN7"))]

# ╔═╡ 41102bd6-a2dc-46c6-9510-50f35ade2c12
md"# Example 3: 2GIS"

# ╔═╡ f028cfab-f18a-43d3-86f6-809833b13d81
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0"))
	@show shape_data_045.aptamer_ids[n_ex]
	
	width = 700
	height = 100
	xticks = 5:5:108

	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	#Makie.scatter!(ax_diff, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -1.5, 1)
	Makie.ylims!(ax_react, -0.5, 6)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA/cossio/SAM/2024/SamApp2024.jl/pluto/Figures/Fig5new_SHAPE_example_2GIS.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═1f9e4648-3b0f-4ce8-acc1-8aa4a51e7242
# ╠═8e34985c-368c-4a4f-99c0-5cc16d195974
# ╠═5e10a3a4-e5bb-441c-8ab5-6b1db64d9ebb
# ╠═c01cb31b-d0c9-4f42-82f6-43f60bb3d67c
# ╠═d3ecc35c-ac26-4b86-a021-f5cf284de8cc
# ╠═d7847339-7164-46f0-88e8-5b8193de577a
# ╠═4e597e09-0311-47d2-9f26-8ec6a1218128
# ╠═5002f57d-c32c-42d9-a350-417e85615943
# ╠═20af2d04-5621-458c-b129-501f6a0da126
# ╠═398d2110-7faa-4878-bf53-f61d6e375f05
# ╠═a94f1ab3-c6d6-427d-8836-d3f1e35f738d
# ╠═d806a548-11a5-41a3-8591-1f52f46ee8c8
# ╠═fce1dfd3-89b0-4ff1-a0f9-c71109a683f8
# ╠═6b943684-49d3-4eee-b5d6-aeb7fe4fef1f
# ╠═4afbfcac-2e59-44d0-a850-5c49c7167720
# ╠═75374067-8150-4721-b83f-79fae3ac12b2
# ╠═dcce26c2-893a-4963-9852-2a7880167da2
# ╠═e0bca8ed-50f4-482d-896a-edb2eeb154fe
# ╠═ba434302-f435-4dd8-82c0-c5b67028c848
# ╠═1fde47a6-3818-4c28-8934-f6a4e1c880f8
# ╠═0bddc5ea-29de-4f27-bb14-e1e889a18993
# ╠═ee13efeb-8568-4583-97fe-65ede489e356
# ╠═b759c7bb-dc94-4dc5-bc83-6d9e29e93748
# ╠═46eddbfb-70af-4bf4-b448-9f21067c6871
# ╠═8eb1f3ce-fbc9-461f-9e2a-ac8707702f9a
# ╠═b9653255-2c4d-4758-9eb9-667937bde77d
# ╠═8fa8e888-3ce8-4726-83c3-7cafe3db9f14
# ╠═eecaf592-337e-4483-acf8-43cf7a68abd0
# ╠═fc3cbbbc-d4f3-4789-894d-76921085e6fa
# ╠═19279ab8-afdd-4ecb-bf46-41e7f0ce956e
# ╠═8056ed26-0adb-4600-9a51-c0273e1ef322
# ╠═5bcfc9a7-147f-4e3d-b1b9-89cb878fe9b0
# ╠═cb10c61d-2066-4d01-a0d5-fe49e6a2e9b8
# ╠═82dd4082-50cc-4eff-9407-cd7808d64688
# ╠═f65cd678-a8c5-495d-b104-de5161cabb5e
# ╠═f1a565af-6d83-427a-be5b-46792243f631
# ╠═9a12a8c9-e86d-4ddd-a6c2-ffa419e06ea2
# ╠═7110c60b-1c56-4678-bad8-1cd1c6de4034
# ╠═3e703ec9-f5ef-4fe2-ac60-9b6898249304
# ╠═f9afb22d-0990-4616-8d03-7a0f75cf3bd0
# ╠═b0163fae-c4f5-43dc-ad63-a8932d3ef3a4
# ╠═646fa540-843f-4738-966e-f0dd69e7241b
# ╠═d93a6822-20e0-417f-a93e-983ad647a750
# ╠═46ff2334-74b3-45e0-8c1e-2d4e8a80de7b
# ╠═26359219-b16b-4a73-bc48-7ce5deea4b3a
# ╠═a4e7ef5c-7462-496c-86bc-9fdb3ad3611d
# ╠═57e2e479-a340-4345-9b52-8f11587344b5
# ╠═ea96c786-cd96-4337-bb0a-20cbddde86d7
# ╠═6fd2910a-bd6e-48c6-acec-9c9a4d3a74b9
# ╠═faa6a32d-6b02-4c85-aaaa-a38855ae119a
# ╠═1aa04802-6a89-47cb-8994-aa9b40319fb5
# ╠═01c3a737-9faa-4b5a-8f5a-47ff3ad121f1
# ╠═0615574c-94c5-4679-bec3-7b952f58dd05
# ╠═9dd3cf23-2841-47df-94ba-861e911bb49a
# ╠═09ed5c06-8ab8-4ea9-b0c5-98b256838071
# ╠═1e5a5d35-3b84-4c42-af1c-310375ba449a
# ╠═23e51bbf-8db0-470f-ba50-d6999480d922
# ╠═ea8e20bc-fd8e-43cd-8ef5-6a653f28f53f
# ╠═9ecff5e8-4b91-4338-aaaa-1200684d3765
# ╠═ff713d18-cb3a-4707-b694-a7acc6f0d1d2
# ╠═9f491a8e-fee1-409c-b2bd-0fbc0be1270c
# ╠═d12fb5cb-1315-4284-bb3b-57deb741cb87
# ╠═4826b171-a3b5-4411-8357-bf2fc13c2f01
# ╠═99def13a-3700-45c6-9127-64754f234957
# ╠═e7ef002b-2567-434d-8e0f-b707ccd389d3
# ╠═af800276-6e2f-4a14-b6f4-d565b0de6f45
# ╠═f88073e3-144f-459a-8910-7ad5a0c2340d
# ╠═1b7a906e-f240-4c53-b7f0-719ce088908f
# ╠═bbd8df1d-ce78-45cc-90a9-2eb1ba112f7b
# ╠═cbd07af7-0688-467b-bf06-6e51f0cd9693
# ╠═41102bd6-a2dc-46c6-9510-50f35ade2c12
# ╠═f028cfab-f18a-43d3-86f6-809833b13d81

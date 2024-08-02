### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 5b37ecc0-9ce5-4e2d-933a-ac9f9c062d0d
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 152889a2-8b96-49a5-8eea-6ed487385438
using BioSequences: LongRNA

# ╔═╡ 005094b2-d36f-4759-933a-2783e4e4f7d9
using DataFrames: DataFrame

# ╔═╡ 311eab7d-c733-472c-a120-58b245acb81c
using Makie: @L_str

# ╔═╡ ef73626e-278b-4a02-9d3d-9cda50c3beae
using NaNStatistics: nansum

# ╔═╡ 9d1d521b-562f-49ea-9c61-7eff2a778045
using Random: bitrand

# ╔═╡ 59f0c590-97ea-4539-8a49-d178cfa27c07
using Statistics: cor

# ╔═╡ baaa2949-7f41-44da-913f-9083748e8466
using Statistics: mean

# ╔═╡ 05720d06-3cfc-4db9-9c9a-184b9b4eeeb7
using StatsBase: countmap

# ╔═╡ b6a715f3-153e-445c-a61e-87382d15a0ba
md"# Imports"

# ╔═╡ 7c92a3c6-481b-4359-bd24-6c908efda46d
import Makie

# ╔═╡ de471d8a-18de-4a6b-86f1-184128c83733
import CairoMakie

# ╔═╡ 2951323f-8e46-46e3-86ec-8aba305d204b
import FASTX

# ╔═╡ 47952a8c-39f5-4450-8850-f3a71e9d341b
import Infernal

# ╔═╡ 615f2f7a-7a0d-4360-a4b7-f3bb343e8cf0
import SamApp2024

# ╔═╡ 358aa00c-8371-4699-aafd-a4def8520b29
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ c5f5306e-7cb3-4e3e-81a0-3a109c460915
import Rfam

# ╔═╡ 58338cfb-c5b9-4263-8ff5-34c4abef2fed
import PlutoUI

# ╔═╡ 537eb2be-632b-4579-af6f-4c83e7abc1ae
import Unitful

# ╔═╡ ee8092aa-9eff-4d96-991b-b594f061a80d
import ViennaRNA

# ╔═╡ 61948391-f70a-45f7-bc65-4a9e367fa0ba
import StatsBase

# ╔═╡ 0f921b35-7342-4e61-8e24-545cf58cb031
import KernelDensity

# ╔═╡ 5603f363-96a1-40e5-8e38-e0d2552177c1
PlutoUI.TableOfContents()

# ╔═╡ 94e3a304-66ef-4825-860f-b6adc222372d
md"# Load data"

# ╔═╡ 7e01165e-41f1-4328-b284-950a065c7738
@show Rfam.get_rfam_directory() Rfam.get_rfam_version();

# ╔═╡ 1dd97772-6f38-4ffc-aecc-03b37a07c763
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20240730_with_pdb();

# ╔═╡ e9055ce5-23a1-4c0b-b929-32d7082a7d86
# All merged data, for the reactivity profiles plots
shape_data_all_merged = SamApp2024.load_shapemapper_data_pierre_demux_20240801_with_pdb_repls_merged();

# ╔═╡ 46d25430-f3fd-43a0-9962-427aaa8b7955
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ c0678b07-73ba-44fa-b5c2-0bb14356e2ed
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 6f10fefb-7c50-4251-9c37-14c8505b808a
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ d3000ccb-3069-464c-a6a8-d822eb3468e6
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 77813a6a-4ceb-4a18-bbd9-387e340389a8
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 1bd16e79-f300-4e03-9136-4157e00ac59e
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 5b807dda-6acb-45ae-aeb2-50f649540a5c
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 502c11cb-e45d-46af-a094-7ea0bf438b0a
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ ee73d339-3b15-461f-9040-5c1ed7d7a566
conds_SAM_all_merged = map(identity, indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_allrep", "SAMAP_1M7_1SAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ fcd284a0-50d3-42ad-ba57-89f568f637c0
conds_Mg_all_merged = map(identity, indexin(["SAMAP_1M7_noSAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 51bb4443-42f4-439f-8255-d1f2fd89477e
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ b156d28e-2459-46b3-9238-9e5ccf532e37
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ ed764596-8ef6-4ca9-9780-aea485c9ac4d
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ d74aac70-c56a-4596-b14e-f926c237b57e
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 4c4819dd-3f82-489f-ab64-4756dc5e0473
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ c99ea59e-015d-490c-96a2-45c7f669e46e
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ bc0664d7-370d-4f62-a4c0-5db8b4e61141
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 5233e4f7-9519-4f54-9462-ed86fed7967b
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ cb9b6efc-68e1-4774-974d-44bd5236eabe
md"# Natural and PDB examples (from Repl. 0)"

# ╔═╡ 0b64262b-0955-48fa-b21c-dc583b766f20
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0];

# ╔═╡ 71e94a83-69c7-4be4-a78f-f9e64c23b509
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0];

# ╔═╡ 2db2fd10-8e88-444a-9a47-44c20d58f328
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ 4c8e11c4-e7b2-450f-92d1-0bd02c914a4f
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 84b73033-032d-4166-9530-77bf2e7a9156
_thresh = log(1)

# ╔═╡ 1390357f-3524-4ec5-abbc-7976ddd988ae
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ 0b59b6a0-ce7d-4a07-8021-c278892629ef
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))

# ╔═╡ 02b44cc7-ab81-486f-a47e-aca031042003
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

# ╔═╡ 714afebd-0f15-4a64-9e59-4d63352432d2
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 40ced149-4220-4e65-a512-31b2aabc14c1
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ 68bbd0f3-a714-4dc9-b127-1d3331f5c1bb
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ cb5be1aa-8452-45c2-a19b-745979ea54ce
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

# ╔═╡ 5703b0d9-c713-4f78-a096-88ba30ab07c9
sum(_responds_sam_yes_rep0), sum(_responds_sam_nop_rep0), sum(_inconclusive_rep0)

# ╔═╡ 7dae3dbc-ac62-4f59-bdc9-1bdb1fa6e4ad
shape_data_045.aptamer_origin[end-1:end]

# ╔═╡ 04aa8a36-f8cc-4038-b79b-ff82740c1fdc
_responds_sam_yes_rep0[shape_data_045.aptamer_origin .== "PDB"]

# ╔═╡ d19b01a0-bd2b-4280-8e50-a7d86a81ccf8
only(_responds_sam_yes_rep0[shape_data_045.aptamer_names .== "APSAMN172"])

# ╔═╡ 996fa88c-512c-4dea-b01b-24db423999c7
only(_responds_sam_yes_rep0[shape_data_045.aptamer_names .== "APSAMN30"])

# ╔═╡ b4de3dcd-dc4d-4a24-a228-48a85dfa8163
only(_responds_sam_yes_rep0[shape_data_045.aptamer_names .== "APSAMN96"]), only(_responds_sam_nop_rep0[shape_data_045.aptamer_names .== "APSAMN96"])

# ╔═╡ 6eb31fd6-59f6-4311-87ed-94a59670a72b
only(_responds_sam_yes_rep0[shape_data_045.aptamer_names .== "APSAMS25"]), only(_responds_sam_nop_rep0[shape_data_045.aptamer_names .== "APSAMS25"])

# ╔═╡ a247dad7-844d-4f9a-ba71-b27a24d3c028
only(_responds_sam_yes_rep0[shape_data_045.aptamer_names .== "APSAMS10"]), only(_responds_sam_nop_rep0[shape_data_045.aptamer_names .== "APSAMS10"])

# ╔═╡ 82c4b4c7-48cc-48ce-86cd-8633e239541f
md"""
# PDB examples (from Repl. 45)

PDB examples are present only in Repl 45 (they have no data in repl. 0).
"""

# ╔═╡ a1eddfa2-2600-4c1c-af16-c179a1bf3c4b
bps_reactivities_rep45 = shape_data_rep45.shape_reactivities[bps, nat_seqs, conds_sam_rep45];

# ╔═╡ 060ad3b4-0ef4-4d36-a26e-1609b27635c6
nps_reactivities_rep45 = shape_data_rep45.shape_reactivities[nps, nat_seqs, conds_sam_rep45];

# ╔═╡ e4c2f662-f6a9-4b9d-b591-ffb8e1ceb6d3
all_reactivities_rep45 = shape_data_rep45.shape_reactivities[:, nat_seqs, conds_sam_rep45];

# ╔═╡ cabe9e9b-c9f4-4a36-bde9-12ee38e7ff9a
shape_stats_rep45 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep45,
    paired_reactivities = bps_reactivities_rep45,
    unpaired_reactivities = nps_reactivities_rep45,
    all_reactivities = all_reactivities_rep45,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ f810b7ef-c682-4cd8-847d-e8adfc54ac42
md"# Reactivity profiles of PDB examples (from Repl. 0)"

# ╔═╡ 2b47c8c3-33a8-43b9-be9f-37f3594846d1
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

# ╔═╡ ad38de09-1a74-413f-8676-3a3c261b7310
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10"))
	_width = 700
	_height = 100

	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]
	
	ax_react_1 = Makie.Axis(
		fig[1,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:4:8, xtrimspine=true, ytrimspine=true, title=shape_data_045.aptamer_ids[n_ex]
	)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_1, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_1, 1:108, _R_mg, color=:gray, step=:center, label="no SAM")
	Makie.stairs!(ax_react_1, 1:108, _R_sam, color=:purple, step=:center, label="with SAM")
	Makie.hidespines!(ax_react_1, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_1)
	#Makie.axislegend(ax_react_1, position=(0.0, -13), framevisible=false)
	
	ax_diff_1 = Makie.Axis(fig[2,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_1, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_1, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff_1, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff_1, 0, 109)
	Makie.hidespines!(ax_diff_1, :r, :t)
	#Makie.hidexdecorations!(ax_diff_1)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react_1, ax_diff_1)
	Makie.ylims!(ax_diff_1, -1.5, 1)
	Makie.ylims!(ax_react_1, -0.5, 8)
	
	Makie.xlims!(ax_react_1, 0.5, 108.5)
	Makie.xlims!(ax_diff_1,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 7d867cb8-23f0-45c0-b3b9-b377e5de73eb
let fig = Makie.Figure()
	n_ex = only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0"))
	_width = 700
	_height = 100

	_R_sam = shape_data_rep0.shape_reactivities[:, n_ex, conds_sam_rep0[3]]
	_R_mg = shape_data_rep0.shape_reactivities[:, n_ex, only(conds_mg_rep0)]
	
	ax_react_1 = Makie.Axis(
		fig[1,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:4:8, xtrimspine=true, ytrimspine=true, title=shape_data_045.aptamer_ids[n_ex]
	)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_1, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_1, 1:108, _R_mg, color=:gray, step=:center, label="no SAM")
	Makie.stairs!(ax_react_1, 1:108, _R_sam, color=:purple, step=:center, label="with SAM")
	Makie.hidespines!(ax_react_1, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_1)
	#Makie.axislegend(ax_react_1, position=(0.0, -13), framevisible=false)
	
	ax_diff_1 = Makie.Axis(fig[2,1]; valign=:bottom, width=_width, height=_height, xticks=5:10:108, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, yticks=-1:1, xtrimspine=true, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_1, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_1, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff_1, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff_1, 0, 109)
	Makie.hidespines!(ax_diff_1, :r, :t)
	#Makie.hidexdecorations!(ax_diff_1)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react_1, ax_diff_1)
	Makie.ylims!(ax_diff_1, -1.5, 1)
	Makie.ylims!(ax_react_1, -0.5, 8)
	
	Makie.xlims!(ax_react_1, 0.5, 108.5)
	Makie.xlims!(ax_diff_1,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ dc5ea24a-704f-49ee-b384-68c5168d4ef6
any(isfinite, shape_data_rep0.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0")), conds_sam_rep0]),
any(isfinite, shape_data_rep0.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0")), conds_mg_rep0]),
any(isfinite, shape_data_rep0.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), conds_sam_rep0]),
any(isfinite, shape_data_rep0.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), conds_mg_rep0])

# ╔═╡ 1943e448-20a3-413b-b96c-13950225d861
any(isfinite, shape_data_rep45.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0")), conds_sam_rep45]),
any(isfinite, shape_data_rep45.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0")), conds_mg_rep45]),
any(isfinite, shape_data_rep45.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), conds_sam_rep45]),
any(isfinite, shape_data_rep45.shape_reactivities[:, only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB10")), conds_mg_rep45])

# ╔═╡ Cell order:
# ╠═b6a715f3-153e-445c-a61e-87382d15a0ba
# ╠═5b37ecc0-9ce5-4e2d-933a-ac9f9c062d0d
# ╠═7c92a3c6-481b-4359-bd24-6c908efda46d
# ╠═de471d8a-18de-4a6b-86f1-184128c83733
# ╠═2951323f-8e46-46e3-86ec-8aba305d204b
# ╠═47952a8c-39f5-4450-8850-f3a71e9d341b
# ╠═615f2f7a-7a0d-4360-a4b7-f3bb343e8cf0
# ╠═358aa00c-8371-4699-aafd-a4def8520b29
# ╠═c5f5306e-7cb3-4e3e-81a0-3a109c460915
# ╠═58338cfb-c5b9-4263-8ff5-34c4abef2fed
# ╠═537eb2be-632b-4579-af6f-4c83e7abc1ae
# ╠═ee8092aa-9eff-4d96-991b-b594f061a80d
# ╠═61948391-f70a-45f7-bc65-4a9e367fa0ba
# ╠═0f921b35-7342-4e61-8e24-545cf58cb031
# ╠═152889a2-8b96-49a5-8eea-6ed487385438
# ╠═005094b2-d36f-4759-933a-2783e4e4f7d9
# ╠═311eab7d-c733-472c-a120-58b245acb81c
# ╠═ef73626e-278b-4a02-9d3d-9cda50c3beae
# ╠═9d1d521b-562f-49ea-9c61-7eff2a778045
# ╠═59f0c590-97ea-4539-8a49-d178cfa27c07
# ╠═baaa2949-7f41-44da-913f-9083748e8466
# ╠═05720d06-3cfc-4db9-9c9a-184b9b4eeeb7
# ╠═5603f363-96a1-40e5-8e38-e0d2552177c1
# ╠═94e3a304-66ef-4825-860f-b6adc222372d
# ╠═7e01165e-41f1-4328-b284-950a065c7738
# ╠═1dd97772-6f38-4ffc-aecc-03b37a07c763
# ╠═e9055ce5-23a1-4c0b-b929-32d7082a7d86
# ╠═46d25430-f3fd-43a0-9962-427aaa8b7955
# ╠═c0678b07-73ba-44fa-b5c2-0bb14356e2ed
# ╠═6f10fefb-7c50-4251-9c37-14c8505b808a
# ╠═d3000ccb-3069-464c-a6a8-d822eb3468e6
# ╠═77813a6a-4ceb-4a18-bbd9-387e340389a8
# ╠═1bd16e79-f300-4e03-9136-4157e00ac59e
# ╠═5b807dda-6acb-45ae-aeb2-50f649540a5c
# ╠═502c11cb-e45d-46af-a094-7ea0bf438b0a
# ╠═ee73d339-3b15-461f-9040-5c1ed7d7a566
# ╠═fcd284a0-50d3-42ad-ba57-89f568f637c0
# ╠═51bb4443-42f4-439f-8255-d1f2fd89477e
# ╠═b156d28e-2459-46b3-9238-9e5ccf532e37
# ╠═ed764596-8ef6-4ca9-9780-aea485c9ac4d
# ╠═d74aac70-c56a-4596-b14e-f926c237b57e
# ╠═4c4819dd-3f82-489f-ab64-4756dc5e0473
# ╠═c99ea59e-015d-490c-96a2-45c7f669e46e
# ╠═bc0664d7-370d-4f62-a4c0-5db8b4e61141
# ╠═5233e4f7-9519-4f54-9462-ed86fed7967b
# ╠═cb9b6efc-68e1-4774-974d-44bd5236eabe
# ╠═0b64262b-0955-48fa-b21c-dc583b766f20
# ╠═71e94a83-69c7-4be4-a78f-f9e64c23b509
# ╠═2db2fd10-8e88-444a-9a47-44c20d58f328
# ╠═4c8e11c4-e7b2-450f-92d1-0bd02c914a4f
# ╠═84b73033-032d-4166-9530-77bf2e7a9156
# ╠═1390357f-3524-4ec5-abbc-7976ddd988ae
# ╠═0b59b6a0-ce7d-4a07-8021-c278892629ef
# ╠═02b44cc7-ab81-486f-a47e-aca031042003
# ╠═714afebd-0f15-4a64-9e59-4d63352432d2
# ╠═40ced149-4220-4e65-a512-31b2aabc14c1
# ╠═68bbd0f3-a714-4dc9-b127-1d3331f5c1bb
# ╠═cb5be1aa-8452-45c2-a19b-745979ea54ce
# ╠═5703b0d9-c713-4f78-a096-88ba30ab07c9
# ╠═7dae3dbc-ac62-4f59-bdc9-1bdb1fa6e4ad
# ╠═04aa8a36-f8cc-4038-b79b-ff82740c1fdc
# ╠═d19b01a0-bd2b-4280-8e50-a7d86a81ccf8
# ╠═996fa88c-512c-4dea-b01b-24db423999c7
# ╠═b4de3dcd-dc4d-4a24-a228-48a85dfa8163
# ╠═6eb31fd6-59f6-4311-87ed-94a59670a72b
# ╠═a247dad7-844d-4f9a-ba71-b27a24d3c028
# ╠═82c4b4c7-48cc-48ce-86cd-8633e239541f
# ╠═a1eddfa2-2600-4c1c-af16-c179a1bf3c4b
# ╠═060ad3b4-0ef4-4d36-a26e-1609b27635c6
# ╠═e4c2f662-f6a9-4b9d-b591-ffb8e1ceb6d3
# ╠═cabe9e9b-c9f4-4a36-bde9-12ee38e7ff9a
# ╠═f810b7ef-c682-4cd8-847d-e8adfc54ac42
# ╠═2b47c8c3-33a8-43b9-be9f-37f3594846d1
# ╠═ad38de09-1a74-413f-8676-3a3c261b7310
# ╠═7d867cb8-23f0-45c0-b3b9-b377e5de73eb
# ╠═dc5ea24a-704f-49ee-b384-68c5168d4ef6
# ╠═1943e448-20a3-413b-b96c-13950225d861

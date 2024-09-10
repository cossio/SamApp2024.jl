### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 910ed40a-a45d-4abe-9719-9f09dff09e6b
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 42d180ab-1092-45ca-ac7f-174dc493f779
using BioSequences: LongRNA

# ╔═╡ c4232a5e-be0d-414c-adc7-e07acce8f779
using Distributions: Gamma

# ╔═╡ 5ec543b4-a285-49f6-bbda-bb71180e6601
using Makie: @L_str

# ╔═╡ 6bda50c6-7d04-4d3e-86a1-ec4de656c246
using NaNStatistics: nanmean

# ╔═╡ db91eaa3-317c-4bef-890b-1b1fa47f2935
using NaNStatistics: nanstd

# ╔═╡ e063dbe6-e7e6-459f-95c1-17696aa0388a
using NaNStatistics: nansum

# ╔═╡ b231f8f3-bad4-4763-8c7d-e335b40a34b6
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 9e8ea7ad-43fa-4ac4-b602-04a6b44cfc07
using Statistics: cor

# ╔═╡ 20eb3376-f702-43ec-92df-1634cf7a98fd
using Statistics: mean

# ╔═╡ 0b31b947-6c71-40a8-afa9-cc23eabed235
md"# Imports"

# ╔═╡ 995559e5-397c-4049-a024-fdabc0ee24a1
import PlutoUI

# ╔═╡ dfe55eab-b72b-43fc-99a8-c22f74f96f9b
import CairoMakie

# ╔═╡ 5e8f32d4-5283-41e7-b6ac-4dff426bd2ae
import Makie

# ╔═╡ a796fa03-f142-4ae2-824e-dd97fd2e746f
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 44002f4c-7bf9-4562-b79b-0d281dba401c
import SamApp2024

# ╔═╡ 6c5308c4-d360-4f54-802b-e9abbaded7b9
import StatsBase

# ╔═╡ 653db4db-12d2-4481-abe4-81af83fd1d10
PlutoUI.TableOfContents()

# ╔═╡ b268fbd9-2f3d-48a1-a2ed-ecb11f998b4b
md"# Load data (Repl.0)"

# ╔═╡ ab147f81-129c-4762-ab05-6843faacf6f6
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ c013080b-8a47-4f16-b3dc-205a3c77173d
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ d5833ba8-3a0b-4434-8ddd-2d498e8c81bc
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ c2091dea-5509-47a2-930a-455483d83bd3
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ fdb4d4ba-6a81-4e07-b889-4cd7fa4e5f00
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 7d5363e5-4a20-4055-9408-33efb3bd8fde
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ fd93ce37-da5c-4dc1-aba4-6128ebab6b2f
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 1e3f1c68-ac2f-4ba9-b62c-6b9880701eae
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 77aea6ba-4153-4309-8f83-9f2c35cc55ad
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ baa4d77c-2c14-400d-a2be-73feb797a665
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 5bba752f-4057-47c6-a45d-391fe08ca8df
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ e4442aea-22c0-4fd4-9c49-9f6eb5047322
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 53ae987f-e5ab-44c5-9fba-aa860ddd2bb8
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 5f3a94a0-c0b7-4175-99f1-012ed9057368
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 90f0cd0f-3396-4e6e-8059-2bd38ad10f6e
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 8cccd2d8-c3e9-43e9-9a12-abe0cdeb3264
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ b46a63bd-de0f-4a6e-bd3c-4ff06b15a505
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ 7d7e7e3e-91eb-4990-bf19-c221f3b33a84
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ ede16ee1-c49b-4485-ac2a-1147f0f57dc2
md"# Load data (500 seqs)"

# ╔═╡ d4c5b09e-f1df-471d-be61-114815dd78f1
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ 89e6fd9c-0bf8-476d-bf73-fee6053f6b8b
conds_sam_500 = [1,2];

# ╔═╡ 13108515-5e21-416b-b8cf-d5b22fd215b0
conds_mg_500 = [4];

# ╔═╡ 1fdf70cd-8de7-4465-a9f9-11d8bd5d1cc3
conds_30C_500 = [6];

# ╔═╡ dcfbf05d-bdca-4579-940e-36b951335d76
bps_reactivities_500 = shape_data_500.shape_reactivities[bps, :, conds_sam_500];

# ╔═╡ 294262b4-25af-4c83-b610-a6216b53b102
nps_reactivities_500 = shape_data_500.shape_reactivities[nps, :, conds_sam_500];

# ╔═╡ 44554cfb-2a5e-4938-bfb3-9ff4a7055802
all_reactivities_500 = shape_data_500.shape_reactivities[:, :, conds_sam_500];

# ╔═╡ 2af1e54c-af9c-4917-94b5-0b781d09cce3
shape_stats_500 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = bps_reactivities_500,
    unpaired_reactivities = nps_reactivities_500,
    all_reactivities = all_reactivities_500,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ bc18fd34-0ea1-4034-8945-96dbb20b4b89
md"# Compute protection scores"

# ╔═╡ 3a24b18e-7c8f-45e2-8e72-bb280e746571
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0];

# ╔═╡ 36d1bfc5-5cb8-4b3d-ace2-1dcd40fd21f3
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0];

# ╔═╡ cbf50471-70bd-4a61-855e-e7b7befe473d
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ 14f3514f-790a-4557-b57f-4c0b3a53953a
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ 5739e646-c62e-4577-b34d-bdea1a6c65e0
_thresh = log(5)

# ╔═╡ df0f1809-8209-4e12-a158-431f7acd0f91
function null_model_FP_shape_data(shape_data)
	shape_data_null = deepcopy(shape_data)
	for n = axes(shape_data.shape_M, 2)
		_is = rand(nps, 108)
		for c = conds_sam_rep0
			shape_data_null.shape_M[:, n, c] .= shape_data.shape_M[_is, n, c]
			shape_data_null.shape_U[:, n, c] .= shape_data.shape_U[_is, n, c]
			shape_data_null.shape_D[:, n, c] .= shape_data.shape_D[_is, n, c]

			shape_data_null.shape_M_stderr[:, n, c] .= shape_data.shape_M_stderr[_is, n, c]
			shape_data_null.shape_U_stderr[:, n, c] .= shape_data.shape_U_stderr[_is, n, c]
			shape_data_null.shape_D_stderr[:, n, c] .= shape_data.shape_D_stderr[_is, n, c]

			shape_data_null.shape_M_depth[:, n, c] .= shape_data.shape_M_depth[_is, n, c]
			shape_data_null.shape_U_depth[:, n, c] .= shape_data.shape_U_depth[_is, n, c]
			shape_data_null.shape_D_depth[:, n, c] .= shape_data.shape_D_depth[_is, n, c]

			shape_data_null.shape_reactivities[:, n, c] .= shape_data.shape_reactivities[_is, n, c]
		end
	end

	return shape_data_null
end

# ╔═╡ 08cb7fa3-421a-429f-bd8d-3fbcb09b6376
function null_model_FN_shape_data(shape_data)
	shape_data_null = deepcopy(shape_data)
	for n = axes(shape_data.shape_M, 2)
		_is = rand(bps, 108)
		for c = conds_sam_rep0
			shape_data_null.shape_M[:, n, c] .= shape_data.shape_M[_is, n, c]
			shape_data_null.shape_U[:, n, c] .= shape_data.shape_U[_is, n, c]
			shape_data_null.shape_D[:, n, c] .= shape_data.shape_D[_is, n, c]

			shape_data_null.shape_M_stderr[:, n, c] .= shape_data.shape_M_stderr[_is, n, c]
			shape_data_null.shape_U_stderr[:, n, c] .= shape_data.shape_U_stderr[_is, n, c]
			shape_data_null.shape_D_stderr[:, n, c] .= shape_data.shape_D_stderr[_is, n, c]

			shape_data_null.shape_M_depth[:, n, c] .= shape_data.shape_M_depth[_is, n, c]
			shape_data_null.shape_U_depth[:, n, c] .= shape_data.shape_U_depth[_is, n, c]
			shape_data_null.shape_D_depth[:, n, c] .= shape_data.shape_D_depth[_is, n, c]

			shape_data_null.shape_reactivities[:, n, c] .= shape_data.shape_reactivities[_is, n, c]
		end
	end

	return shape_data_null
end

# ╔═╡ 0bc97164-25ec-46d9-80fd-14b42b652fc0
md"# Make table"

# ╔═╡ 0d091096-900e-4f3e-b4f4-252f99523c7d


# ╔═╡ Cell order:
# ╠═0b31b947-6c71-40a8-afa9-cc23eabed235
# ╠═910ed40a-a45d-4abe-9719-9f09dff09e6b
# ╠═995559e5-397c-4049-a024-fdabc0ee24a1
# ╠═dfe55eab-b72b-43fc-99a8-c22f74f96f9b
# ╠═5e8f32d4-5283-41e7-b6ac-4dff426bd2ae
# ╠═a796fa03-f142-4ae2-824e-dd97fd2e746f
# ╠═44002f4c-7bf9-4562-b79b-0d281dba401c
# ╠═6c5308c4-d360-4f54-802b-e9abbaded7b9
# ╠═42d180ab-1092-45ca-ac7f-174dc493f779
# ╠═c4232a5e-be0d-414c-adc7-e07acce8f779
# ╠═5ec543b4-a285-49f6-bbda-bb71180e6601
# ╠═6bda50c6-7d04-4d3e-86a1-ec4de656c246
# ╠═db91eaa3-317c-4bef-890b-1b1fa47f2935
# ╠═e063dbe6-e7e6-459f-95c1-17696aa0388a
# ╠═b231f8f3-bad4-4763-8c7d-e335b40a34b6
# ╠═9e8ea7ad-43fa-4ac4-b602-04a6b44cfc07
# ╠═20eb3376-f702-43ec-92df-1634cf7a98fd
# ╠═653db4db-12d2-4481-abe4-81af83fd1d10
# ╠═b268fbd9-2f3d-48a1-a2ed-ecb11f998b4b
# ╠═ab147f81-129c-4762-ab05-6843faacf6f6
# ╠═c013080b-8a47-4f16-b3dc-205a3c77173d
# ╠═d5833ba8-3a0b-4434-8ddd-2d498e8c81bc
# ╠═c2091dea-5509-47a2-930a-455483d83bd3
# ╠═fdb4d4ba-6a81-4e07-b889-4cd7fa4e5f00
# ╠═7d5363e5-4a20-4055-9408-33efb3bd8fde
# ╠═fd93ce37-da5c-4dc1-aba4-6128ebab6b2f
# ╠═1e3f1c68-ac2f-4ba9-b62c-6b9880701eae
# ╠═77aea6ba-4153-4309-8f83-9f2c35cc55ad
# ╠═baa4d77c-2c14-400d-a2be-73feb797a665
# ╠═5bba752f-4057-47c6-a45d-391fe08ca8df
# ╠═e4442aea-22c0-4fd4-9c49-9f6eb5047322
# ╠═53ae987f-e5ab-44c5-9fba-aa860ddd2bb8
# ╠═5f3a94a0-c0b7-4175-99f1-012ed9057368
# ╠═90f0cd0f-3396-4e6e-8059-2bd38ad10f6e
# ╠═8cccd2d8-c3e9-43e9-9a12-abe0cdeb3264
# ╠═b46a63bd-de0f-4a6e-bd3c-4ff06b15a505
# ╠═7d7e7e3e-91eb-4990-bf19-c221f3b33a84
# ╠═ede16ee1-c49b-4485-ac2a-1147f0f57dc2
# ╠═d4c5b09e-f1df-471d-be61-114815dd78f1
# ╠═89e6fd9c-0bf8-476d-bf73-fee6053f6b8b
# ╠═13108515-5e21-416b-b8cf-d5b22fd215b0
# ╠═1fdf70cd-8de7-4465-a9f9-11d8bd5d1cc3
# ╠═dcfbf05d-bdca-4579-940e-36b951335d76
# ╠═294262b4-25af-4c83-b610-a6216b53b102
# ╠═44554cfb-2a5e-4938-bfb3-9ff4a7055802
# ╠═2af1e54c-af9c-4917-94b5-0b781d09cce3
# ╠═bc18fd34-0ea1-4034-8945-96dbb20b4b89
# ╠═3a24b18e-7c8f-45e2-8e72-bb280e746571
# ╠═36d1bfc5-5cb8-4b3d-ace2-1dcd40fd21f3
# ╠═cbf50471-70bd-4a61-855e-e7b7befe473d
# ╠═14f3514f-790a-4557-b57f-4c0b3a53953a
# ╠═5739e646-c62e-4577-b34d-bdea1a6c65e0
# ╠═df0f1809-8209-4e12-a158-431f7acd0f91
# ╠═08cb7fa3-421a-429f-bd8d-3fbcb09b6376
# ╠═0bc97164-25ec-46d9-80fd-14b42b652fc0
# ╠═0d091096-900e-4f3e-b4f4-252f99523c7d

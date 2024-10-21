### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 910ed40a-a45d-4abe-9719-9f09dff09e6b
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 42d180ab-1092-45ca-ac7f-174dc493f779
using BioSequences: LongRNA

# ╔═╡ 5ec543b4-a285-49f6-bbda-bb71180e6601
using Makie: @L_str

# ╔═╡ e063dbe6-e7e6-459f-95c1-17696aa0388a
using NaNStatistics: nansum

# ╔═╡ b231f8f3-bad4-4763-8c7d-e335b40a34b6
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ eeca4e4a-7dff-49a1-bb87-8514a19043d2
using DataFrames: DataFrame

# ╔═╡ 0b31b947-6c71-40a8-afa9-cc23eabed235
md"# Imports"

# ╔═╡ 995559e5-397c-4049-a024-fdabc0ee24a1
import PlutoUI

# ╔═╡ a796fa03-f142-4ae2-824e-dd97fd2e746f
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 44002f4c-7bf9-4562-b79b-0d281dba401c
import SamApp2024

# ╔═╡ b56b53a7-fe6d-4194-8a39-9b5121ed856d
import CSV

# ╔═╡ 653db4db-12d2-4481-abe4-81af83fd1d10
PlutoUI.TableOfContents()

# ╔═╡ f90147bf-7dde-4b7f-84b4-3ca7eea9fdd6
md"# General data"

# ╔═╡ e54bbdee-8d7b-48ac-a68a-0039b9b9ba36
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 2bc7ae6b-9eb3-444a-816b-83f689e29773
# Hallmark sites
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ 6e70dd1c-dec3-471d-8003-8f02349a1aba
# signifcance threshold for protection scores
_thresh = log(5)

# ╔═╡ b268fbd9-2f3d-48a1-a2ed-ecb11f998b4b
md"# Load data (Repl.0)"

# ╔═╡ ab147f81-129c-4762-ab05-6843faacf6f6
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ c013080b-8a47-4f16-b3dc-205a3c77173d
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ c2091dea-5509-47a2-930a-455483d83bd3
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ fdb4d4ba-6a81-4e07-b889-4cd7fa4e5f00
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 7d5363e5-4a20-4055-9408-33efb3bd8fde
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ baa4d77c-2c14-400d-a2be-73feb797a665
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 53ae987f-e5ab-44c5-9fba-aa860ddd2bb8
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 5f3a94a0-c0b7-4175-99f1-012ed9057368
inf_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 90f0cd0f-3396-4e6e-8059-2bd38ad10f6e
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 8cccd2d8-c3e9-43e9-9a12-abe0cdeb3264
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ b46a63bd-de0f-4a6e-bd3c-4ff06b15a505
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ 122ee22d-0adc-4ac1-bfc4-a3db3446dec9
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 7f52c28a-6813-4bfe-b303-2293dc8c15c9
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ bd2f02a4-98a0-469e-8922-695f9b50c832
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 49905d69-de0b-44af-86d1-9fd366498c0a
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ 6c5416af-9730-40b9-b70e-3aaa5048e705
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3));

# ╔═╡ ed61e8ec-7a49-4f16-9a90-b340e47cb643
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3));

# ╔═╡ db3a8c33-0e4e-4f71-ae15-9a41813cc365
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 83dcd13b-aa52-4dd6-b224-85eb65b10f4b
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ 7d7e7e3e-91eb-4990-bf19-c221f3b33a84
aptamer_rbm_energies_rep0 = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ ede16ee1-c49b-4485-ac2a-1147f0f57dc2
md"# Load data (500 seqs)"

# ╔═╡ d4c5b09e-f1df-471d-be61-114815dd78f1
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ 89e6fd9c-0bf8-476d-bf73-fee6053f6b8b
conds_sam_500, conds_mg_500, conds_30C_500 = [1,2], [4], [6];

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

# ╔═╡ 4546ffc8-d928-4740-a8a7-b1e99e093288
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :,  conds_mg_500]; dim=(1,3));

# ╔═╡ b83bb091-66d5-4958-b930-76eb6e947f86
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ a6c5b767-0740-4cc1-93cf-51c2efa96309
_responds_sam_yes_500 = (x_mg_500 .< -_thresh) .& (x_sam_500 .> +_thresh);

# ╔═╡ 0ac21b72-6937-4942-a547-3dda848a39fe
_responds_sam_nop_500 = (x_mg_500 .> +_thresh) .| (x_sam_500 .< -_thresh);

# ╔═╡ df8b0931-bab7-44c3-a8da-eaac2fb55fb3
aptamer_rbm_energies_500 = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_500.aligned_sequences
];

# ╔═╡ 0bc97164-25ec-46d9-80fd-14b42b652fc0
md"# Make table"

# ╔═╡ 406c30f3-a087-411e-ac43-b73d826897cf
_responsive_sam_rep0 = ifelse.(_responds_sam_yes_rep0, "Responsive", ifelse.(_responds_sam_nop_rep0, "Non-responsive", "Inconclusive"))

# ╔═╡ dfddbf66-aa46-4220-9398-4456cf45028b
_responsive_sam_500 = ifelse.(_responds_sam_yes_500, "Responsive", ifelse.(_responds_sam_nop_500, "Non-responsive", "Inconclusive"))

# ╔═╡ 0d091096-900e-4f3e-b4f4-252f99523c7d
df = DataFrame(;
    aligned_sequences = [shape_data_rep0.aligned_sequences; shape_data_500.aligned_sequences],
    aptamer_origin = [shape_data_rep0.aptamer_origin; shape_data_500.aptamer_origin],
	experiment = [fill("Experiment_1", length(shape_data_rep0.aligned_sequences)); fill("Experiment_2", length(shape_data_500.aligned_sequences))],
    responsive = [_responsive_sam_rep0; _responsive_sam_500],
	RBM_score = [-aptamer_rbm_energies_rep0; -aptamer_rbm_energies_500],
	Protect_Score_Hallmark_Mg = [x_mg_rep0; x_mg_500],
	Protect_Score_Hallmark_SAM = [x_sam_rep0; x_sam_500]
)

# ╔═╡ d14b0358-1839-489e-85aa-7edded875a88
md"# Export table"

# ╔═╡ 15e163fe-8765-487e-8883-0348b0ecb455
CSV.write("Data/2024-09-11_Suppl_Table.csv", df)

# ╔═╡ Cell order:
# ╠═0b31b947-6c71-40a8-afa9-cc23eabed235
# ╠═910ed40a-a45d-4abe-9719-9f09dff09e6b
# ╠═995559e5-397c-4049-a024-fdabc0ee24a1
# ╠═a796fa03-f142-4ae2-824e-dd97fd2e746f
# ╠═44002f4c-7bf9-4562-b79b-0d281dba401c
# ╠═b56b53a7-fe6d-4194-8a39-9b5121ed856d
# ╠═42d180ab-1092-45ca-ac7f-174dc493f779
# ╠═5ec543b4-a285-49f6-bbda-bb71180e6601
# ╠═e063dbe6-e7e6-459f-95c1-17696aa0388a
# ╠═b231f8f3-bad4-4763-8c7d-e335b40a34b6
# ╠═eeca4e4a-7dff-49a1-bb87-8514a19043d2
# ╠═653db4db-12d2-4481-abe4-81af83fd1d10
# ╠═f90147bf-7dde-4b7f-84b4-3ca7eea9fdd6
# ╠═e54bbdee-8d7b-48ac-a68a-0039b9b9ba36
# ╠═2bc7ae6b-9eb3-444a-816b-83f689e29773
# ╠═6e70dd1c-dec3-471d-8003-8f02349a1aba
# ╠═b268fbd9-2f3d-48a1-a2ed-ecb11f998b4b
# ╠═ab147f81-129c-4762-ab05-6843faacf6f6
# ╠═c013080b-8a47-4f16-b3dc-205a3c77173d
# ╠═c2091dea-5509-47a2-930a-455483d83bd3
# ╠═fdb4d4ba-6a81-4e07-b889-4cd7fa4e5f00
# ╠═7d5363e5-4a20-4055-9408-33efb3bd8fde
# ╠═baa4d77c-2c14-400d-a2be-73feb797a665
# ╠═53ae987f-e5ab-44c5-9fba-aa860ddd2bb8
# ╠═5f3a94a0-c0b7-4175-99f1-012ed9057368
# ╠═90f0cd0f-3396-4e6e-8059-2bd38ad10f6e
# ╠═8cccd2d8-c3e9-43e9-9a12-abe0cdeb3264
# ╠═b46a63bd-de0f-4a6e-bd3c-4ff06b15a505
# ╠═122ee22d-0adc-4ac1-bfc4-a3db3446dec9
# ╠═7f52c28a-6813-4bfe-b303-2293dc8c15c9
# ╠═bd2f02a4-98a0-469e-8922-695f9b50c832
# ╠═49905d69-de0b-44af-86d1-9fd366498c0a
# ╠═6c5416af-9730-40b9-b70e-3aaa5048e705
# ╠═ed61e8ec-7a49-4f16-9a90-b340e47cb643
# ╠═db3a8c33-0e4e-4f71-ae15-9a41813cc365
# ╠═83dcd13b-aa52-4dd6-b224-85eb65b10f4b
# ╠═7d7e7e3e-91eb-4990-bf19-c221f3b33a84
# ╠═ede16ee1-c49b-4485-ac2a-1147f0f57dc2
# ╠═d4c5b09e-f1df-471d-be61-114815dd78f1
# ╠═89e6fd9c-0bf8-476d-bf73-fee6053f6b8b
# ╠═dcfbf05d-bdca-4579-940e-36b951335d76
# ╠═294262b4-25af-4c83-b610-a6216b53b102
# ╠═44554cfb-2a5e-4938-bfb3-9ff4a7055802
# ╠═2af1e54c-af9c-4917-94b5-0b781d09cce3
# ╠═4546ffc8-d928-4740-a8a7-b1e99e093288
# ╠═b83bb091-66d5-4958-b930-76eb6e947f86
# ╠═a6c5b767-0740-4cc1-93cf-51c2efa96309
# ╠═0ac21b72-6937-4942-a547-3dda848a39fe
# ╠═df8b0931-bab7-44c3-a8da-eaac2fb55fb3
# ╠═0bc97164-25ec-46d9-80fd-14b42b652fc0
# ╠═406c30f3-a087-411e-ac43-b73d826897cf
# ╠═dfddbf66-aa46-4220-9398-4456cf45028b
# ╠═0d091096-900e-4f3e-b4f4-252f99523c7d
# ╠═d14b0358-1839-489e-85aa-7edded875a88
# ╠═15e163fe-8765-487e-8883-0348b0ecb455

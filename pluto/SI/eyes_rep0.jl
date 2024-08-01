### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 935fc628-348a-11ef-38cf-758397ff9cf6
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 1f9db468-3570-42d0-9d69-de3a4d8fae2c
using BioSequences: LongRNA

# ╔═╡ 852b0931-fc64-4915-b816-3d714f13eff0
using DataFrames: DataFrame

# ╔═╡ 1d213d26-0f44-4607-b688-4c31f037b24f
using Distributions: Gamma

# ╔═╡ d0d0e370-44e2-4b1d-8968-e4ca8be37997
using Distributions: logpdf

# ╔═╡ eaf9fbbe-2b68-4afc-81f6-fd1b7b3e613a
using Distributions: pdf

# ╔═╡ e9d79be8-1360-4e8f-8d45-50ad6fd8979f
using Distributions: Poisson

# ╔═╡ e5be3a23-c5b0-4a1b-bbf9-1403dd6f53cf
using LinearAlgebra: Diagonal

# ╔═╡ 34f63831-ba17-4589-b863-3ef2483f4240
using LinearAlgebra: eigen

# ╔═╡ ba8b6c25-cae1-40dc-bd05-761a51c5f66c
using Makie: @L_str

# ╔═╡ 6fae1e7c-f3ee-4b0c-b542-c623b4dcee78
using NaNStatistics: nanmean, nansum

# ╔═╡ 6a20d879-aeb9-44a6-b831-294b38afd3ca
using Random: bitrand

# ╔═╡ 11585ed8-ee1d-483a-8bc2-fc8efb1aa5cb
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 42d62fa4-850a-47a9-8494-b09abe978132
using Statistics: cor

# ╔═╡ ced9ef56-8787-47f6-a665-1b9114c45ef6
using Statistics: mean

# ╔═╡ 658c3424-a2e2-4d28-a475-6dd34ea67b04
using Statistics: middle

# ╔═╡ 3f29b445-5d71-4b5a-a054-eb98ec7d24c4
using StatsBase: countmap

# ╔═╡ cc0b4412-9d86-4f71-b614-7c2a19385657
using Unitful: ustrip

# ╔═╡ d008c260-2165-4217-97a7-d4d47e3372df
md"# Imports"

# ╔═╡ 099f78f0-4877-4fa7-8907-ab152950ef62
import CairoMakie

# ╔═╡ a07bc9ea-169f-4baf-acf3-9734d0a7a922
import Makie

# ╔═╡ cf7b1a2e-9aff-4dca-a1fe-73f3dc1068aa
import SamApp2024

# ╔═╡ 65b859b6-cbef-4513-b5eb-79b2abd75b1b
import CSV

# ╔═╡ 09686bc0-2a00-46c9-96cf-f66b42267888
import FASTX

# ╔═╡ dc8b78bf-38d9-4093-9cf2-8abcb2b90727
import HDF5

# ╔═╡ 89861771-8d02-431f-875f-0c3b63189463
import Infernal

# ╔═╡ 591be198-8444-476d-8448-34e5e23502b7
import KernelDensity

# ╔═╡ e93c97a1-e0fa-4565-ab34-e4a1d819481f
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 6998ff93-78e1-4b61-8d42-d301bd6e9ae3
import Rfam

# ╔═╡ 81a16243-3915-4135-ad89-5e2cc4b9b9d3
import StatsBase

# ╔═╡ 18572c57-560d-4503-8e68-a55bb8595e00
import ViennaRNA

# ╔═╡ af740c9c-bfb1-4a5a-afa1-badf804be1ec
import ViennaRNA_jll

# ╔═╡ ae85207c-862e-4221-a983-37044952e9c1
import PlutoUI

# ╔═╡ ca9ee2cd-c8b4-49ec-bac9-e2288b4fd877
PlutoUI.TableOfContents()

# ╔═╡ b51759a0-df1d-4b4e-af11-51618bb31574
md"# Load data"

# ╔═╡ 5be4163b-374c-4901-8376-7c1fe6ac0066
@show Rfam.get_rfam_directory();

# ╔═╡ 6a735f9e-2f44-4bcb-8828-58d3983ceec0
@show Rfam.get_rfam_version();

# ╔═╡ dfda0eab-fc96-4c4e-8ff2-89d3e4fd0cf6
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 215c7d66-760d-46fe-a74d-fa84fd0eeca3
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ f27542a4-c369-4d45-bcca-3bf235cd8d39
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 17f5d388-ae37-4203-8429-f016bd57d12e
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 60275ad9-e808-4b07-8dad-0ff0c3c90b87
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ eb513cdd-f513-405e-a43f-44d929d43498
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 09dcc509-a01a-4158-bbeb-e19f161a5793
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 1e9e40ab-a425-44cf-a3ba-3298fe56c894
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 4e56953c-c760-482a-a0f5-ae0c8c737056
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 7f37a57a-52d2-4f5b-b0a6-f1aab4c2a179
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 932fc140-fa8b-4e6d-8a17-b7db2886a24e
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ ecb48a4c-9242-4fc8-acb5-e02b1ac185cb
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 58dd48a3-6d97-4573-9e86-efec48914c50
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 1b4e7280-e537-4433-bbea-91125a3fdcce
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 7486bb16-522c-4e9f-bd72-c6e658bc6226
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 16f73e28-0b10-4d66-aefa-2f59d0d77cb2
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ a12bfa61-e00b-44b4-b257-20514d3a6043
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ 9c18f5d0-7cd7-417f-aa67-0a43874fbbca
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0];

# ╔═╡ a5ac34f9-984d-48b5-a8ad-6db78a2f8e9d
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0];

# ╔═╡ f447b9de-41cd-40d7-bcfa-cea1d4c200fa
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ 8b702321-fd78-4e72-a40d-fb149b36a8c0
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 660ddb1b-537d-4a4a-87c1-12cee75cfe33
_thresh = log(5)

# ╔═╡ 224fdc0a-859e-4e0e-8291-d1d48d24970c
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ 58e4e854-93e9-422f-b799-95d896c8249c
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))

# ╔═╡ e963f6e3-835e-4a70-b8d2-b6aa1b0c1da1
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

# ╔═╡ e803f158-379f-4e85-bc06-c8dad2534e37
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ b82f538b-c4a7-4164-a583-8297de242b30
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ 1119960b-b5db-4e64-b24c-93c46e331ea9
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ 4eae5eaa-37d8-43ff-be2a-c2d3ccf9a01d
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

# ╔═╡ f949a0c1-064d-4f61-9b3b-92a93a65f289
_responds_sam_yes_rep0

# ╔═╡ 27fec5a4-13c6-461f-9f90-ec136fb25c3b
_selection = ["APSAMN172", "APSAMN75", "APSAMN96", "APSAMS25", "APSAMS10"]

# ╔═╡ 47634614-5886-4de6-8608-da57354ce289
[_responds_sam_yes_rep0[only(findall(==(_s), shape_stats_rep0.aptamer_names))] for _s = _selection]

# ╔═╡ 48d6dda0-9885-4653-85eb-ce42e6a1da4f
[_responds_sam_nop_rep0[only(findall(==(_s), shape_stats_rep0.aptamer_names))] for _s = _selection]

# ╔═╡ 67d9cc7e-c984-4069-91be-e471500ca68d
[_inconclusive_rep0[only(findall(==(_s), shape_stats_rep0.aptamer_names))] for _s = _selection]

# ╔═╡ Cell order:
# ╠═d008c260-2165-4217-97a7-d4d47e3372df
# ╠═935fc628-348a-11ef-38cf-758397ff9cf6
# ╠═099f78f0-4877-4fa7-8907-ab152950ef62
# ╠═a07bc9ea-169f-4baf-acf3-9734d0a7a922
# ╠═cf7b1a2e-9aff-4dca-a1fe-73f3dc1068aa
# ╠═65b859b6-cbef-4513-b5eb-79b2abd75b1b
# ╠═09686bc0-2a00-46c9-96cf-f66b42267888
# ╠═dc8b78bf-38d9-4093-9cf2-8abcb2b90727
# ╠═89861771-8d02-431f-875f-0c3b63189463
# ╠═591be198-8444-476d-8448-34e5e23502b7
# ╠═e93c97a1-e0fa-4565-ab34-e4a1d819481f
# ╠═6998ff93-78e1-4b61-8d42-d301bd6e9ae3
# ╠═81a16243-3915-4135-ad89-5e2cc4b9b9d3
# ╠═18572c57-560d-4503-8e68-a55bb8595e00
# ╠═af740c9c-bfb1-4a5a-afa1-badf804be1ec
# ╠═ae85207c-862e-4221-a983-37044952e9c1
# ╠═1f9db468-3570-42d0-9d69-de3a4d8fae2c
# ╠═852b0931-fc64-4915-b816-3d714f13eff0
# ╠═1d213d26-0f44-4607-b688-4c31f037b24f
# ╠═d0d0e370-44e2-4b1d-8968-e4ca8be37997
# ╠═eaf9fbbe-2b68-4afc-81f6-fd1b7b3e613a
# ╠═e9d79be8-1360-4e8f-8d45-50ad6fd8979f
# ╠═e5be3a23-c5b0-4a1b-bbf9-1403dd6f53cf
# ╠═34f63831-ba17-4589-b863-3ef2483f4240
# ╠═ba8b6c25-cae1-40dc-bd05-761a51c5f66c
# ╠═6fae1e7c-f3ee-4b0c-b542-c623b4dcee78
# ╠═6a20d879-aeb9-44a6-b831-294b38afd3ca
# ╠═11585ed8-ee1d-483a-8bc2-fc8efb1aa5cb
# ╠═42d62fa4-850a-47a9-8494-b09abe978132
# ╠═ced9ef56-8787-47f6-a665-1b9114c45ef6
# ╠═658c3424-a2e2-4d28-a475-6dd34ea67b04
# ╠═3f29b445-5d71-4b5a-a054-eb98ec7d24c4
# ╠═cc0b4412-9d86-4f71-b614-7c2a19385657
# ╠═ca9ee2cd-c8b4-49ec-bac9-e2288b4fd877
# ╠═b51759a0-df1d-4b4e-af11-51618bb31574
# ╠═5be4163b-374c-4901-8376-7c1fe6ac0066
# ╠═6a735f9e-2f44-4bcb-8828-58d3983ceec0
# ╠═dfda0eab-fc96-4c4e-8ff2-89d3e4fd0cf6
# ╠═215c7d66-760d-46fe-a74d-fa84fd0eeca3
# ╠═f27542a4-c369-4d45-bcca-3bf235cd8d39
# ╠═17f5d388-ae37-4203-8429-f016bd57d12e
# ╠═60275ad9-e808-4b07-8dad-0ff0c3c90b87
# ╠═eb513cdd-f513-405e-a43f-44d929d43498
# ╠═09dcc509-a01a-4158-bbeb-e19f161a5793
# ╠═1e9e40ab-a425-44cf-a3ba-3298fe56c894
# ╠═4e56953c-c760-482a-a0f5-ae0c8c737056
# ╠═7f37a57a-52d2-4f5b-b0a6-f1aab4c2a179
# ╠═932fc140-fa8b-4e6d-8a17-b7db2886a24e
# ╠═ecb48a4c-9242-4fc8-acb5-e02b1ac185cb
# ╠═58dd48a3-6d97-4573-9e86-efec48914c50
# ╠═1b4e7280-e537-4433-bbea-91125a3fdcce
# ╠═7486bb16-522c-4e9f-bd72-c6e658bc6226
# ╠═16f73e28-0b10-4d66-aefa-2f59d0d77cb2
# ╠═a12bfa61-e00b-44b4-b257-20514d3a6043
# ╠═9c18f5d0-7cd7-417f-aa67-0a43874fbbca
# ╠═a5ac34f9-984d-48b5-a8ad-6db78a2f8e9d
# ╠═f447b9de-41cd-40d7-bcfa-cea1d4c200fa
# ╠═8b702321-fd78-4e72-a40d-fb149b36a8c0
# ╠═660ddb1b-537d-4a4a-87c1-12cee75cfe33
# ╠═224fdc0a-859e-4e0e-8291-d1d48d24970c
# ╠═58e4e854-93e9-422f-b799-95d896c8249c
# ╠═e963f6e3-835e-4a70-b8d2-b6aa1b0c1da1
# ╠═e803f158-379f-4e85-bc06-c8dad2534e37
# ╠═b82f538b-c4a7-4164-a583-8297de242b30
# ╠═1119960b-b5db-4e64-b24c-93c46e331ea9
# ╠═4eae5eaa-37d8-43ff-be2a-c2d3ccf9a01d
# ╠═f949a0c1-064d-4f61-9b3b-92a93a65f289
# ╠═27fec5a4-13c6-461f-9f90-ec136fb25c3b
# ╠═47634614-5886-4de6-8608-da57354ce289
# ╠═48d6dda0-9885-4653-85eb-ce42e6a1da4f
# ╠═67d9cc7e-c984-4069-91be-e471500ca68d

### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f08d8cfc-711e-47c7-a5a9-f70e3c0cb1ac
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 35b53fe4-c79d-48b0-bbe6-b8c6d469c3c0
using BioSequences: LongRNA

# ╔═╡ b738e180-90d8-419b-95bd-c383e4b19012
using DataFrames: DataFrame

# ╔═╡ 4c6704f4-8802-4f5c-80c5-25d8520fcc25
using Makie: @L_str

# ╔═╡ 7f4e47e0-980e-440c-9e81-ea25a0149f34
using NaNStatistics: nansum

# ╔═╡ 8ff5bd07-2672-4565-b14f-deadc0a3b2d3
using Random: bitrand

# ╔═╡ c4f91ffb-cc6e-4347-804a-887ca45be249
using Statistics: cor

# ╔═╡ a718fd95-2963-4b35-bddd-d3b53a2faf52
using Statistics: mean

# ╔═╡ 4244cffe-81f2-41ff-950d-97e973947a32
using StatsBase: countmap

# ╔═╡ d42dae48-fdc2-11ef-191c-b7b143d538a5
md"# Imports"

# ╔═╡ dd8fcacd-6feb-439f-8495-861731ef26db
import Makie

# ╔═╡ a80399d9-96e5-45a1-8d0c-0ab364781e3e
import CairoMakie

# ╔═╡ 3553a742-fbfb-4553-abe0-17a4044f6eae
import FASTX

# ╔═╡ 059b139e-3d06-41d5-a89a-6990db4ab9fa
import Infernal

# ╔═╡ e71e2863-3a8e-4854-812e-ec14b5458350
import SamApp2024

# ╔═╡ 57681d3b-7d23-4b9c-ab1c-d24b72845654
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 4aa377c9-05fc-45df-9ab2-0500952b41f7
import Rfam

# ╔═╡ 59e38478-1323-4047-b591-c0a788e898c8
import PlutoUI

# ╔═╡ 41fb8639-2472-4341-98a1-b8f11475e31e
import Unitful

# ╔═╡ 030731fb-dd82-4d0f-a9d7-026af8378d20
import ViennaRNA

# ╔═╡ 1c2478bc-bdfc-4596-a053-b5736cbe6504
import StatsBase

# ╔═╡ 6bc48dfb-9009-4aee-8f47-3623b8834ee7
import KernelDensity

# ╔═╡ 6155ff9f-af5f-4c1a-8e6e-965f11efd285
PlutoUI.TableOfContents()

# ╔═╡ ed1f893e-29e0-446f-bbd4-dc5516614aa1
md"# General data"

# ╔═╡ 1a341857-b3e4-4b7e-8172-13c399b0fca3
_sites = SamApp2024.hallmark_sites_20230507

# ╔═╡ e0367f59-2f02-4a09-8018-e71317695f3b
md"# DMS data"

# ╔═╡ af8c9beb-f844-4238-8967-2dbff72ac27c
dms_df = SamApp2024.load_dms_data_sequences_table_20250303_with_aligned_sequences()

# ╔═╡ aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
dms_data = SamApp2024.load_dms_data_20250303()

# ╔═╡ e00ec7b5-4f06-4200-9bd0-9e486e4322bc
dms_data_primers = dms_df.primer[[only(findall(dms_df.name .== name)) for name = dms_data.aptamer_names]]

# ╔═╡ 4875da0e-9898-4498-8b29-cb60282e2cc2
unique(dms_data_primers)

# ╔═╡ 2b18233e-d2ae-485e-9ad7-17e20a26fee7
num_sites, num_seqs, num_conds = size(dms_data.shape_reactivities)

# ╔═╡ 59a5fee0-59f1-49da-9940-8827c5a23b99
bps, nps, pks = SamApp2024.RF00162_sites_paired()

# ╔═╡ 025f0c01-0e0d-4a50-b9cd-becd304a7a42
all_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=1:num_sites for n=1:num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ 11a25b40-e29b-4f48-8be7-2cc91feb8491
bps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=bps for n=1:num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ cd5b25ff-6b98-41f3-af8d-f07665f4547f
nps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=nps for n=1:num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ c69fbd19-69f4-4290-9a92-7bd5d57f79ea
dms_stats = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = dms_data,
    paired_reactivities = bps_dms_reactivities,
    unpaired_reactivities = nps_dms_reactivities,
    all_reactivities = all_dms_reactivities,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ d7ba4366-918e-4b62-8756-691d07f10bf8
dms_data.conditions

# ╔═╡ 8fe64552-4504-4d55-ab9e-275cdbc51293
x_mg_dms = [
	begin
		if ismissing(dms_data.aligned_sequence[n])
			NaN
		else
			reactivities = [dms_stats.shape_log_odds[i,n,2] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
			# if isempty(reactivities)
			# 	NaN
			# else
			# 	nansum(reactivities)
			# end
			nansum(reactivities)
		end
	end for n=1:400
];

# ╔═╡ fe8dbdeb-5c0b-41a4-b32e-ef742fde4146
x_sam_dms = [
	begin
		if ismissing(dms_data.aligned_sequence[n])
			NaN
		else
			reactivities = [dms_stats.shape_log_odds[i,n,1] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
			# if isempty(reactivities)
			# 	NaN
			# else
			# 	nansum(reactivities)
			# end
			nansum(reactivities)
		end
	end for n=1:400
];

# ╔═╡ 7959344a-1e02-45ef-9ef3-71e892e5379c
_thresh_dms = log(5)

# ╔═╡ c5ffed2c-c1f3-41a8-82c1-a5842fba5ce4
_responds_sam_yes_dms = (x_mg_dms .< -_thresh_dms) .& (x_sam_dms .> +_thresh_dms);

# ╔═╡ 010c18f0-8aea-4c91-b4d1-300117449683
_responds_sam_nop_dms = (x_mg_dms .> +_thresh_dms) .| (x_sam_dms .< -_thresh_dms);

# ╔═╡ f22df52d-cb6b-4fe8-9a16-bc1488a0122f
_responds_sam_inconcl_dms = ((!).(_responds_sam_yes_dms)) .& ((!).(_responds_sam_nop_dms));

# ╔═╡ c4fe44f9-fd26-4450-9479-afc596abe4ee
md"# Load SHAPE data (rep.0)"

# ╔═╡ f1f0126d-e212-47bc-8acc-736bb5905c0a
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ b7fe3fe6-a590-414e-8cbf-05daee93bab1
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 110aeb3e-2eab-44cf-ba3a-238075218ae4
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 8caabc4e-2ac7-4ee8-a2c3-d162078904f5
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ e9d1d644-e4db-4a8a-b1d2-f8ffb843215b
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 1c2a3b46-ac85-49d9-8789-b899e4ea1009
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 3f572a36-7a0a-48f6-9868-6b1c113c614a
inf_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 2e2d454d-55e4-4bd7-bdd3-823b1515a3c1
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ e2d9671d-175f-44b5-a1e0-6a041f3bc91a
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ cb7fc980-abf0-44e7-86e3-3c85b0435f64
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ 0414c6b2-105b-40fd-8bdf-d171de1fb60a
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 777b984e-febb-4ac3-b1c4-3487ab250ba5
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 0cdbd1e1-ab5d-434a-901d-cbf5042c8d6e
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 964fa067-a76d-4798-89a4-1b443a0ed48d
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ b38219f2-2f5b-46d9-97af-8e52a04fa886
_thresh_rep0 = log(5)

# ╔═╡ 48c1b560-bfa3-4e6f-9a62-482e56a3228f
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))

# ╔═╡ ca87405c-3a9f-4b62-8a36-d1db46edc13e
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

# ╔═╡ 61613274-8e04-46ba-90dc-aa5b4e0c661f
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh_rep0) .& (x_sam_rep0 .> +_thresh_rep0);

# ╔═╡ 70b61f1d-c2a7-4f3b-bbfd-aabcca807470
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh_rep0) .| (x_sam_rep0 .< -_thresh_rep0);

# ╔═╡ 9471fa84-538f-4b54-8e70-97eb279f0f7e
_responds_sam_inconcl_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ c6aaed41-40ae-43ec-94c0-e70679251316
md"# Load SHAPE data (500 aptamers)"

# ╔═╡ 022735fe-bae2-4b21-a786-3d4e373eb24f
conds_sam_500 = [1,2];

# ╔═╡ 804b0de7-57b2-4ef1-8fe5-af03466207f1
conds_mg_500 = [4];

# ╔═╡ 16341f0f-4e41-446c-9095-dd90f21538da
conds_30C_500 = [6];

# ╔═╡ 13d43c06-6a78-4f81-b735-02b6dc6c1662
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ 9e1fca3b-3733-4fc2-b00e-beb8c318713d
bps_reactivities_500 = shape_data_500.shape_reactivities[bps, :, conds_sam_500];

# ╔═╡ ba975787-84e4-479b-ad49-f7e54b90c5ce
nps_reactivities_500 = shape_data_500.shape_reactivities[nps, :, conds_sam_500];

# ╔═╡ 709eeec8-b9d0-4a63-91d1-9d6aec643ce5
all_reactivities_500 = shape_data_500.shape_reactivities[:, :, conds_sam_500];

# ╔═╡ fad3d549-76e9-4ad2-8c3c-83f68270ecd1
shape_stats_500 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = bps_reactivities_500,
    unpaired_reactivities = nps_reactivities_500,
    all_reactivities = all_reactivities_500,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ 261e6ff5-1434-4ea5-9294-293f5d903236
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :,  conds_mg_500]; dim=(1,3));

# ╔═╡ 8de6a4d9-df94-498d-9a25-674631eec2e9
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ a2f33ae5-6b7d-41ea-800e-be5d70058b2e
_thresh_500 = log(5)

# ╔═╡ fa17604e-683e-4a18-8e88-fbad8cd70fc4
_responds_sam_yes_500 = (x_mg_500 .< -_thresh_500) .& (x_sam_500 .> +_thresh_500);

# ╔═╡ 59cf3f69-683f-4a81-a60f-e5205d9b046c
_responds_sam_nop_500 = (x_mg_500 .> +_thresh_500) .| (x_sam_500 .< -_thresh_500);

# ╔═╡ f515605a-f8e9-4c86-9388-ef5565520074
_responds_sam_inconcl_500 = ((!).(_responds_sam_yes_500)) .& ((!).(_responds_sam_nop_500));

# ╔═╡ cf94281c-f97d-4966-bd3f-ddad7b8c21d0
md"# Comparison DMS vs. SHAPE Repl.0"

# ╔═╡ 10e1166d-5e07-45aa-8adc-1b92f370d26d
_dms_rep0_indices = indexin(dms_data.aligned_sequence, shape_data_rep0.aligned_sequences)

# ╔═╡ 74fffba4-b0d1-4f15-b430-81e1d677d2f4
_responds_sam_yes_rep0_dms_seqs = [isnothing(n) ? nothing : _responds_sam_yes_rep0[n] for n = _dms_rep0_indices]

# ╔═╡ 29729f40-dbc3-4d28-946e-cee8ec62bd81
_responds_sam_nop_rep0_dms_seqs = [isnothing(n) ? nothing : _responds_sam_nop_rep0[n] for n = _dms_rep0_indices]

# ╔═╡ 84be17bf-eafa-404b-a264-b64e6206a725
_responds_sam_inconcl_rep0_dms_seqs = [isnothing(n) ? nothing : _responds_sam_inconcl_rep0[n] for n = _dms_rep0_indices]

# ╔═╡ 9597c913-72e6-41f2-9072-983b03a278d5
# DMS yes and SHAPE (rep.0) yes
count(_responds_sam_yes_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_yes_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ bfb20f7d-e892-4254-a2b6-8700d0958072
# DMS yes and SHAPE (rep.0) nop
count(_responds_sam_yes_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 54b11c1f-271b-461c-9662-e838c56d4e87
# DMS yes and SHAPE (rep.0) inconclusive
count(_responds_sam_yes_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 72e0ed21-0f65-4ffa-9cb1-02e45bb7a474
# DMS nop and SHAPE (rep.0) yes
count(_responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_yes_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ d25c2bc7-930c-4133-ab2a-2ef29c0b6e82
# DMS nop and SHAPE (rep.0) nop
count(_responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ fda5a42a-509e-41e6-98e6-598261190276
# DMS nop and SHAPE (rep.0) inconclusive
count(_responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 33ca29be-cb6c-47dd-82f5-43c4596bb84e
# DMS inconclusive and SHAPE (rep.0) yes
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_yes_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 25b42825-6964-4ab8-b043-3740112b7dfe
# DMS inconclusive and SHAPE (rep.0) nop
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 0a698578-2c83-42b5-b7fa-a80fb0887bba
# DMS inconclusive and SHAPE (rep.0) inconclusive
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 755c1c17-20b4-4ee9-96ef-201473651a58
# How many RBM are confirmed ?
count(_responds_sam_yes_dms[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 16d309f1-480a-4b60-90cf-a9ef9bbe67a1
# How many RBM are confirmed ?
count(_responds_sam_nop_dms[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 5e6d4cc3-e1a9-4535-aced-835d77303d83
# How many RBM are confirmed ?
count(_responds_sam_inconcl_dms[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)] .&& _responds_sam_inconcl_rep0_dms_seqs[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)])

# ╔═╡ 154f87e7-3199-4232-980c-f6f29b1db1fe
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 5660fe93-b5b5-4bdf-9f80-918815c766bb
# RBM only
[
	begin
		count(dms[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)] .&& rep0[rbm_seqs_rep0 ∩ findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# primer 1
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[findall(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[1])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ c7a94bc0-c9d2-45e5-ac3a-75e25e5c288f
# primer 2
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[findall(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[2])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 2b0d56a5-de91-45af-83aa-5ca3de9b8427
# primer 3
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[findall(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[3])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 6a825260-2887-41d4-9e20-fd66ec0e2e8c
# primer 4
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[findall(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[4])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_dms_seqs, _responds_sam_nop_rep0_dms_seqs, _responds_sam_inconcl_rep0_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 6f7f502f-3a03-4f94-ba92-75fc71427564
unique(dms_data_primers)

# ╔═╡ b4acfb01-3865-4a34-8b25-3e474a306da4
length(dms_data_primers)

# ╔═╡ 9163a506-bc9f-4877-978f-ef854b858be1
length(_dms_rep0_indices)

# ╔═╡ 6aa691a1-6e62-4077-83f6-e168810f4ec2
md"# Comparison DMS vs. SHAPE 500"

# ╔═╡ 8a604d45-5890-4578-8a82-01065faf1a1c
_dms_500_indices = indexin(dms_data.aligned_sequence, map(string, shape_data_500.aligned_sequences))

# ╔═╡ cc2bd6ed-6e09-453d-91ea-552aebe1ca1d
_responds_sam_yes_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_yes_500[n] for n = _dms_500_indices]

# ╔═╡ deb8707e-3e53-4f95-ac51-e6ebf5bffe9c
_responds_sam_nop_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_nop_500[n] for n = _dms_500_indices]

# ╔═╡ f366dabf-e8b7-4417-9a1e-99412ec6c582
_responds_sam_inconcl_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_inconcl_500[n] for n = _dms_500_indices]

# ╔═╡ d56ea5f4-59a3-4f96-87d3-e435303e7cf2
# DMS yes and SHAPE (500) yes
count(_responds_sam_yes_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_yes_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ 4089033a-e18e-4be3-b3e9-b22d23de0656
# DMS yes and SHAPE (500) nop
count(_responds_sam_yes_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_nop_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ c6a5fbf0-971c-48e0-9bee-3512a788d11c
# DMS yes and SHAPE (500) inconclusive
count(_responds_sam_yes_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_inconcl_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ e9e5fe1d-2331-453e-bca8-e587f62da74d
# DMS nop and SHAPE (500) yes
count(_responds_sam_nop_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_yes_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ c19502ff-1dfb-477c-8069-5a2c272740fb
# DMS nop and SHAPE (500) nop
count(_responds_sam_nop_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_nop_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ d33911ab-0a8d-4d9f-b1f3-d4eda14ec64c
# DMS nop and SHAPE (500) inconclusive
count(_responds_sam_nop_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_inconcl_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ 2273495b-7150-4ffb-adb0-c26de1762c25
# DMS inconclusive and SHAPE (500) yes
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_yes_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ 7fd2cf30-9848-4f61-b641-82572f83591b
# DMS inconclusive and SHAPE (500) nop
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_nop_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ 0b8aec74-cf6e-4b18-8929-8ac6faeddce2
# DMS inconclusive and SHAPE (500) inconclusive
count(_responds_sam_inconcl_dms[findall(!isnothing, _dms_500_indices)] .&& _responds_sam_inconcl_500_dms_seqs[findall(!isnothing, _dms_500_indices)])

# ╔═╡ 6ce57a7b-0e16-43f8-a840-d2ab4e9e03d3
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500_dms_seqs, _responds_sam_nop_500_dms_seqs, _responds_sam_inconcl_500_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 1fda3f8c-4477-4039-8636-14b2ae0a5df3
# primer 1
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[findall(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[1])[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500_dms_seqs, _responds_sam_nop_500_dms_seqs, _responds_sam_inconcl_500_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 61738118-55ec-4c08-85f8-de345c103ee6
# primer 2
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[findall(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[2])[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500_dms_seqs, _responds_sam_nop_500_dms_seqs, _responds_sam_inconcl_500_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ d2553930-1e08-4ece-9ae7-23d1256bd2c7
# primer 3
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[findall(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[3])[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500_dms_seqs, _responds_sam_nop_500_dms_seqs, _responds_sam_inconcl_500_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ d8e143e4-5ef4-489c-9107-1cbdd1d78028
# primer 4
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[findall(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== unique(dms_data_primers)[4])[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500_dms_seqs, _responds_sam_nop_500_dms_seqs, _responds_sam_inconcl_500_dms_seqs), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ Cell order:
# ╠═d42dae48-fdc2-11ef-191c-b7b143d538a5
# ╠═f08d8cfc-711e-47c7-a5a9-f70e3c0cb1ac
# ╠═dd8fcacd-6feb-439f-8495-861731ef26db
# ╠═a80399d9-96e5-45a1-8d0c-0ab364781e3e
# ╠═3553a742-fbfb-4553-abe0-17a4044f6eae
# ╠═059b139e-3d06-41d5-a89a-6990db4ab9fa
# ╠═e71e2863-3a8e-4854-812e-ec14b5458350
# ╠═57681d3b-7d23-4b9c-ab1c-d24b72845654
# ╠═4aa377c9-05fc-45df-9ab2-0500952b41f7
# ╠═59e38478-1323-4047-b591-c0a788e898c8
# ╠═41fb8639-2472-4341-98a1-b8f11475e31e
# ╠═030731fb-dd82-4d0f-a9d7-026af8378d20
# ╠═1c2478bc-bdfc-4596-a053-b5736cbe6504
# ╠═6bc48dfb-9009-4aee-8f47-3623b8834ee7
# ╠═35b53fe4-c79d-48b0-bbe6-b8c6d469c3c0
# ╠═b738e180-90d8-419b-95bd-c383e4b19012
# ╠═4c6704f4-8802-4f5c-80c5-25d8520fcc25
# ╠═7f4e47e0-980e-440c-9e81-ea25a0149f34
# ╠═8ff5bd07-2672-4565-b14f-deadc0a3b2d3
# ╠═c4f91ffb-cc6e-4347-804a-887ca45be249
# ╠═a718fd95-2963-4b35-bddd-d3b53a2faf52
# ╠═4244cffe-81f2-41ff-950d-97e973947a32
# ╠═6155ff9f-af5f-4c1a-8e6e-965f11efd285
# ╠═ed1f893e-29e0-446f-bbd4-dc5516614aa1
# ╠═1a341857-b3e4-4b7e-8172-13c399b0fca3
# ╠═e0367f59-2f02-4a09-8018-e71317695f3b
# ╠═af8c9beb-f844-4238-8967-2dbff72ac27c
# ╠═aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
# ╠═e00ec7b5-4f06-4200-9bd0-9e486e4322bc
# ╠═4875da0e-9898-4498-8b29-cb60282e2cc2
# ╠═2b18233e-d2ae-485e-9ad7-17e20a26fee7
# ╠═59a5fee0-59f1-49da-9940-8827c5a23b99
# ╠═025f0c01-0e0d-4a50-b9cd-becd304a7a42
# ╠═11a25b40-e29b-4f48-8be7-2cc91feb8491
# ╠═cd5b25ff-6b98-41f3-af8d-f07665f4547f
# ╠═c69fbd19-69f4-4290-9a92-7bd5d57f79ea
# ╠═d7ba4366-918e-4b62-8756-691d07f10bf8
# ╠═8fe64552-4504-4d55-ab9e-275cdbc51293
# ╠═fe8dbdeb-5c0b-41a4-b32e-ef742fde4146
# ╠═7959344a-1e02-45ef-9ef3-71e892e5379c
# ╠═c5ffed2c-c1f3-41a8-82c1-a5842fba5ce4
# ╠═010c18f0-8aea-4c91-b4d1-300117449683
# ╠═f22df52d-cb6b-4fe8-9a16-bc1488a0122f
# ╠═c4fe44f9-fd26-4450-9479-afc596abe4ee
# ╠═f1f0126d-e212-47bc-8acc-736bb5905c0a
# ╠═b7fe3fe6-a590-414e-8cbf-05daee93bab1
# ╠═110aeb3e-2eab-44cf-ba3a-238075218ae4
# ╠═8caabc4e-2ac7-4ee8-a2c3-d162078904f5
# ╠═e9d1d644-e4db-4a8a-b1d2-f8ffb843215b
# ╠═1c2a3b46-ac85-49d9-8789-b899e4ea1009
# ╠═3f572a36-7a0a-48f6-9868-6b1c113c614a
# ╠═2e2d454d-55e4-4bd7-bdd3-823b1515a3c1
# ╠═e2d9671d-175f-44b5-a1e0-6a041f3bc91a
# ╠═cb7fc980-abf0-44e7-86e3-3c85b0435f64
# ╠═0414c6b2-105b-40fd-8bdf-d171de1fb60a
# ╠═777b984e-febb-4ac3-b1c4-3487ab250ba5
# ╠═0cdbd1e1-ab5d-434a-901d-cbf5042c8d6e
# ╠═964fa067-a76d-4798-89a4-1b443a0ed48d
# ╠═b38219f2-2f5b-46d9-97af-8e52a04fa886
# ╠═48c1b560-bfa3-4e6f-9a62-482e56a3228f
# ╠═ca87405c-3a9f-4b62-8a36-d1db46edc13e
# ╠═61613274-8e04-46ba-90dc-aa5b4e0c661f
# ╠═70b61f1d-c2a7-4f3b-bbfd-aabcca807470
# ╠═9471fa84-538f-4b54-8e70-97eb279f0f7e
# ╠═c6aaed41-40ae-43ec-94c0-e70679251316
# ╠═022735fe-bae2-4b21-a786-3d4e373eb24f
# ╠═804b0de7-57b2-4ef1-8fe5-af03466207f1
# ╠═16341f0f-4e41-446c-9095-dd90f21538da
# ╠═13d43c06-6a78-4f81-b735-02b6dc6c1662
# ╠═9e1fca3b-3733-4fc2-b00e-beb8c318713d
# ╠═ba975787-84e4-479b-ad49-f7e54b90c5ce
# ╠═709eeec8-b9d0-4a63-91d1-9d6aec643ce5
# ╠═fad3d549-76e9-4ad2-8c3c-83f68270ecd1
# ╠═261e6ff5-1434-4ea5-9294-293f5d903236
# ╠═8de6a4d9-df94-498d-9a25-674631eec2e9
# ╠═a2f33ae5-6b7d-41ea-800e-be5d70058b2e
# ╠═fa17604e-683e-4a18-8e88-fbad8cd70fc4
# ╠═59cf3f69-683f-4a81-a60f-e5205d9b046c
# ╠═f515605a-f8e9-4c86-9388-ef5565520074
# ╠═cf94281c-f97d-4966-bd3f-ddad7b8c21d0
# ╠═10e1166d-5e07-45aa-8adc-1b92f370d26d
# ╠═74fffba4-b0d1-4f15-b430-81e1d677d2f4
# ╠═29729f40-dbc3-4d28-946e-cee8ec62bd81
# ╠═84be17bf-eafa-404b-a264-b64e6206a725
# ╠═9597c913-72e6-41f2-9072-983b03a278d5
# ╠═bfb20f7d-e892-4254-a2b6-8700d0958072
# ╠═54b11c1f-271b-461c-9662-e838c56d4e87
# ╠═72e0ed21-0f65-4ffa-9cb1-02e45bb7a474
# ╠═d25c2bc7-930c-4133-ab2a-2ef29c0b6e82
# ╠═fda5a42a-509e-41e6-98e6-598261190276
# ╠═33ca29be-cb6c-47dd-82f5-43c4596bb84e
# ╠═25b42825-6964-4ab8-b043-3740112b7dfe
# ╠═0a698578-2c83-42b5-b7fa-a80fb0887bba
# ╠═755c1c17-20b4-4ee9-96ef-201473651a58
# ╠═16d309f1-480a-4b60-90cf-a9ef9bbe67a1
# ╠═5e6d4cc3-e1a9-4535-aced-835d77303d83
# ╠═154f87e7-3199-4232-980c-f6f29b1db1fe
# ╠═5660fe93-b5b5-4bdf-9f80-918815c766bb
# ╠═40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# ╠═c7a94bc0-c9d2-45e5-ac3a-75e25e5c288f
# ╠═2b0d56a5-de91-45af-83aa-5ca3de9b8427
# ╠═6a825260-2887-41d4-9e20-fd66ec0e2e8c
# ╠═6f7f502f-3a03-4f94-ba92-75fc71427564
# ╠═b4acfb01-3865-4a34-8b25-3e474a306da4
# ╠═9163a506-bc9f-4877-978f-ef854b858be1
# ╠═6aa691a1-6e62-4077-83f6-e168810f4ec2
# ╠═8a604d45-5890-4578-8a82-01065faf1a1c
# ╠═cc2bd6ed-6e09-453d-91ea-552aebe1ca1d
# ╠═deb8707e-3e53-4f95-ac51-e6ebf5bffe9c
# ╠═f366dabf-e8b7-4417-9a1e-99412ec6c582
# ╠═d56ea5f4-59a3-4f96-87d3-e435303e7cf2
# ╠═4089033a-e18e-4be3-b3e9-b22d23de0656
# ╠═c6a5fbf0-971c-48e0-9bee-3512a788d11c
# ╠═e9e5fe1d-2331-453e-bca8-e587f62da74d
# ╠═c19502ff-1dfb-477c-8069-5a2c272740fb
# ╠═d33911ab-0a8d-4d9f-b1f3-d4eda14ec64c
# ╠═2273495b-7150-4ffb-adb0-c26de1762c25
# ╠═7fd2cf30-9848-4f61-b641-82572f83591b
# ╠═0b8aec74-cf6e-4b18-8929-8ac6faeddce2
# ╠═6ce57a7b-0e16-43f8-a840-d2ab4e9e03d3
# ╠═1fda3f8c-4477-4039-8636-14b2ae0a5df3
# ╠═61738118-55ec-4c08-85f8-de345c103ee6
# ╠═d2553930-1e08-4ece-9ae7-23d1256bd2c7
# ╠═d8e143e4-5ef4-489c-9107-1cbdd1d78028

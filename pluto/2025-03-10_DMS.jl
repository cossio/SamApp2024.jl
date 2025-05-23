### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ 6b5e7730-952d-4533-b6a7-7caa14d82fb0
using NaNStatistics: nanmean, nanstd

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
md"# Load DMS data"

# ╔═╡ af8c9beb-f844-4238-8967-2dbff72ac27c
dms_df = SamApp2024.load_dms_data_sequences_table_20250303_with_aligned_sequences()

# ╔═╡ aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
dms_data = SamApp2024.load_dms_data_20250303()

# ╔═╡ e00ec7b5-4f06-4200-9bd0-9e486e4322bc
dms_data_primers = dms_df.primer[[only(findall(dms_df.name .== name)) for name = dms_data.aptamer_names]]

# ╔═╡ 9763e16f-1fc6-4c79-8fed-176f847d9898
length(dms_data_primers)

# ╔═╡ ffa01cd2-cc36-403d-86e8-c3a7e68ac97b
length(dms_data.aptamer_names)

# ╔═╡ 2c0519f1-367d-4fdd-a346-ad550c2e80c6
dms_data 

# ╔═╡ 8fc06bf7-a380-4698-8971-e96df77764b9
length(dms_df.name)

# ╔═╡ df8d1d78-6cd2-40a5-a739-221ed6d03df8
dms_data.aptamer_names ⊆ dms_df.name

# ╔═╡ 5663156d-3342-4cac-9446-1fee39923ce4
dms_data.aptamer_names == dms_df.name[1:400]

# ╔═╡ 4875da0e-9898-4498-8b29-cb60282e2cc2
unique(dms_data_primers)

# ╔═╡ 50c9a9c5-0253-4cc5-aaa5-d58472e8920c
for p = unique(dms_data_primers)
	println(p, "\t", count(!ismissing(dms_df.aligned_sequence) .&& dms_df.primer .== p))
end

# ╔═╡ 5f64c1d1-b429-4946-a122-a02fcb7ff1f6
dms_df.aligned_sequence

# ╔═╡ 8ab280ee-8a1a-495a-ac75-dee2ca199759
dms_df

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

# ╔═╡ 02133d97-7e8a-4d8b-b212-c9ea18fcb9e9
dms_AC_count = [
	begin
		if ismissing(dms_data.aligned_sequence[n])
			NaN
		else
			length([i for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')])
		end
	end for n=1:400
]

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

# ╔═╡ 21f4d2dd-ad45-485f-bcc9-21a735f31c57
unique(shape_data_500.aptamer_origin)

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
# indices of DMS probed sequences in Rep0 (or nothing if it is not in Rep0)
_dms_rep0_indices = indexin(dms_data.aptamer_names, shape_data_rep0.aptamer_names)

# ╔═╡ c40c54d5-b448-41a1-b038-9e82cb804368
unique(shape_data_045.aptamer_origin)

# ╔═╡ 154f87e7-3199-4232-980c-f6f29b1db1fe
[ # all
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 80306132-cc00-4b7c-a012-3890aa00793f
[ # RBM only
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_rbm")[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 0d79ddcc-f6cc-405a-8a1c-a6d3336c1476
[ # natural
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 1d20ceef-c8ef-420d-9c31-e04db5d90a9b
[ # CM
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_inf")[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ eb047b2f-f822-4e03-94bd-aabab68f9610
unique(dms_data_primers)

# ╔═╡ 529b46ff-b9cd-45c9-b202-6b2faeacb40c
SamApp2024.primers20250607()

# ╔═╡ 40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# per primers
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== SamApp2024.primers20250607()[5])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 6aa691a1-6e62-4077-83f6-e168810f4ec2
md"# Comparison DMS vs. SHAPE 500"

# ╔═╡ 8a604d45-5890-4578-8a82-01065faf1a1c
# indices of DMS probed sequences in 500 (or nothing if it is not in 500)
_dms_500_indices = indexin(dms_data.aligned_sequence, map(string, shape_data_500.aligned_sequences[1:450]))

# ╔═╡ f62fd44f-b0b8-42a7-b77f-5e5ab259293b
string.(shape_data_500.aligned_sequences[filter(!isnothing, _dms_500_indices)]) == dms_data.aligned_sequence[findall(!isnothing, _dms_500_indices)]

# ╔═╡ 500a63b1-b4c8-4624-9ddb-107808167ef8
unique(shape_data_500.aptamer_origin)

# ╔═╡ cc2bd6ed-6e09-453d-91ea-552aebe1ca1d
_responds_sam_yes_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_yes_500[n] for n = _dms_500_indices]

# ╔═╡ deb8707e-3e53-4f95-ac51-e6ebf5bffe9c
_responds_sam_nop_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_nop_500[n] for n = _dms_500_indices]

# ╔═╡ f366dabf-e8b7-4417-9a1e-99412ec6c582
_responds_sam_inconcl_500_dms_seqs = [isnothing(n) ? nothing : _responds_sam_inconcl_500[n] for n = _dms_500_indices]

# ╔═╡ 6ce57a7b-0e16-43f8-a840-d2ab4e9e03d3
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[filter(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 4aecd441-932b-44cb-a11b-3866f29f45d0
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (shape_data_500.aptamer_origin .== "infernal")[filter(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 1fda3f8c-4477-4039-8636-14b2ae0a5df3
# per primer
[
	begin
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== SamApp2024.primers20250607()[1])[findall(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ cf9a262a-42ac-40a2-b818-0a020974a6b7
SamApp2024.primers20250607()[5]

# ╔═╡ e7b53f31-137d-43cb-8d80-bb9ddf420828
md"# Read depths"

# ╔═╡ 7c7c7ef4-c613-4868-ac35-57a498765ef6
nanmean(dms_data.shape_M_depth; dim=(2,3))

# ╔═╡ 66f55016-2e10-4590-8d65-f24609f25f48
unique(dms_data_primers)

# ╔═╡ ee9ff992-672e-46fa-87e1-c9c8adb87845
primers = SamApp2024.primers20250607()

# ╔═╡ 5f372f6c-9c6d-44b1-8b29-681b4383d820
SamApp2024.primers20250607()[1] ∈ unique(dms_data_primers)

# ╔═╡ 86651e08-1f33-4fc1-8780-888714d9c86c
let fig = Makie.Figure()
	groups = "primer" .* string.(1:5)
	colors = [:orange, :purple, :red, :blue, :teal]

	_width = 1000
	_height = 200
	
	ax = Makie.Axis(fig[1,1]; width=_width, height=_height, xticks=5:5:108, title="Modified condition", xlabel="Position         ", ylabel="Read depth")
	for (primer, gr, color) = zip(SamApp2024.primers20250607(), groups, colors)
		μ = nanmean(dms_data.shape_M_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		σ = nanstd(dms_data.shape_M_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 120)
	Makie.axislegend(ax)

	ax = Makie.Axis(fig[2,1]; width=_width, height=_height, xticks=5:5:108, title="Unmodified condition", xlabel="Position         ", ylabel="Read depth")
	for (primer, gr, color) = zip(SamApp2024.primers20250607(), groups, colors)
		μ = nanmean(dms_data.shape_U_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		σ = nanstd(dms_data.shape_U_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 120)
	Makie.axislegend(ax)
	
	ax = Makie.Axis(fig[3,1]; width=_width, height=_height, xticks=5:5:108, title="Denatured condition", xlabel="Position         ", ylabel="Read depth")
	for (primer, gr, color) = zip(SamApp2024.primers20250607(), groups, colors)
		μ = nanmean(dms_data.shape_D_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		σ = nanstd(dms_data.shape_D_depth[:, dms_data_primers .== primer, :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 120)
	Makie.axislegend(ax)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ a8c4c868-8d9f-4027-a265-700fbb19b895
let fig = Makie.Figure()
	groups = "primer" .* string.(1:5)
	colors = [:orange, :purple, :red, :blue, :teal]

	ax_succ = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate", xticks=[0, 2e4, 4e4])
	ax_incl = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate", xticks=[0, 2e4, 4e4])

	for (primer, color) = zip(SamApp2024.primers20250607(), colors)
		Makie.scatter!(ax_succ, [nanmean(dms_data.shape_M_depth[:, dms_data_primers .== primer, :])], [mean(_responds_sam_yes_dms[dms_data_primers .== primer])]; color, markersize=15)
		Makie.scatter!(ax_incl, [nanmean(dms_data.shape_M_depth[:, dms_data_primers .== primer, :])], [mean(_responds_sam_inconcl_dms[dms_data_primers .== primer])]; color, markersize=15)
	end

	for ax = (ax_succ,ax_incl)
		Makie.xlims!(ax, 0, 6e4)
		Makie.ylims!(ax, 0.05, 0.35)
	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 03b6232d-295c-4e37-be8b-72ada0e1a05e
length(dms_data_primers)

# ╔═╡ 86771e5d-a4b2-424b-8e72-2939148a4451
unique(shape_data_500.aptamer_origin)

# ╔═╡ 532fb8ca-296a-4afc-b0ba-a9a40e18b03e
length(dms_data.aptamer_names)

# ╔═╡ b2fb16b0-93a8-450f-bde3-7d19a2ab6399
length([n for n = dms_data.aptamer_names if startswith(n, "APSAMN")])

# ╔═╡ 7de7f351-2b52-43de-980c-947304970843
length([n for n = dms_data.aptamer_names if startswith(n, "APSAMS")])

# ╔═╡ 150b2b5e-6f77-420e-b4e5-c1c683044c4c
aptamer_names_500 = "APSAM-S2-" .* lpad.(0:499, 3, '0')

# ╔═╡ eba33d20-66dd-4cc4-b980-cd7013026c89
dms_aptamer_origin = [
	begin
		if n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_rbm"]
			"rbm"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_seed70"]
			"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_full30"]
				"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_inf"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "rbm"]
			"rbm"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "infernal"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "Infrared"]
			"Infrared"
		end
	end for n = dms_data.aptamer_names
]

# ╔═╡ 0e21f914-aec2-416d-a8fa-293767a3d12c
let fig = Makie.Figure()
	groups = "primer" .* string.(1:5)
	colors = [:orange, :purple, :red, :blue, :teal]

	for (n, seq_type) = enumerate(unique(dms_aptamer_origin))
		ax_succ = Makie.Axis(fig[1,n]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate", xticks=[0, 3e4, 6e4], title=seq_type)
		ax_incl = Makie.Axis(fig[2,n]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate", xticks=[0, 3e4, 6e4])

		for (n, (primer, color)) = enumerate(zip(SamApp2024.primers20250607(), colors))
			read_depth = [nanmean(dms_data.shape_M_depth[:, (dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type), :])]

			println("No. sequences for primer $n $primer $color, $seq_type: " , count((dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type)))
			
			Makie.scatter!(ax_succ, read_depth, [mean(_responds_sam_yes_dms[(dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type)])]; 
						   color, markersize=15)
			Makie.scatter!(ax_incl, read_depth, [mean(_responds_sam_inconcl_dms[(dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type)])]; 
						   color, markersize=15)
		end

		for ax = (ax_succ, ax_incl)
			Makie.xlims!(ax, 0, 8e4)
			Makie.ylims!(ax, -0.05, 0.55)
		end
	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 9026821c-a1eb-4d45-96b4-28384701d1ad
let fig = Makie.Figure()
	groups = "primer" .* string.(1:5)
	colors = [:orange, :purple, :red, :blue, :teal]

	for (n, seq_type) = enumerate(unique(dms_aptamer_origin))
		ax_succ = Makie.Axis(fig[1,n]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate / (1 - Incl. rate)", xticks=[0, 3e4, 6e4], title=seq_type)
		ax_incl = Makie.Axis(fig[2,n]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate", xticks=[0, 3e4, 6e4])

		for (primer, color) = zip(SamApp2024.primers20250607(), colors)
			read_depth = [nanmean(dms_data.shape_M_depth[:, (dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type), :])]
			response_rate = [mean(_responds_sam_yes_dms[(dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type)])]
			inconclusive_rate = [mean(_responds_sam_inconcl_dms[(dms_data_primers .== primer) .&& (dms_aptamer_origin .== seq_type)])]
			
			Makie.scatter!(ax_succ, read_depth, response_rate ./ (1 .- inconclusive_rate); color, markersize=15)
			Makie.scatter!(ax_incl, read_depth, inconclusive_rate; color, markersize=15)
		end

		for ax = (ax_succ,ax_incl)
			Makie.xlims!(ax, 0, 8e4)
			Makie.ylims!(ax, -0.05, 0.6)
		end
	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 9bf2fc7d-abd7-45a7-bd0d-455b86d7c4bb
let fig = Makie.Figure()
	colors = [:orange, :purple, :red, :blue, :teal]
	
	ax = Makie.Axis(fig[1,1]; width=230, height=230, xlabel="Read depth (M)", ylabel="Inconclusives / Total", xticks=[1e4, 3e4, 5e4])
	for (m, (primer, color)) = enumerate(zip(SamApp2024.primers20250607(), colors))
		read_depth = nanmean(dms_data.shape_M_depth[:, (dms_data_primers .== primer) .&& (dms_aptamer_origin .!= "Infrared"), :])
		inconclusive_rate = mean(_responds_sam_inconcl_dms[(dms_data_primers .== primer) .&& (dms_aptamer_origin .!= "Infrared")])
		Makie.scatter!(ax, read_depth, inconclusive_rate; color, label="Primer $m", markersize=15)
	end
	Makie.xlims!(ax, 0, 7e4)
	Makie.ylims!(ax, -0.05, 0.45)

	fig[1,2] = Makie.Legend(fig, ax, "Primers", framevisible = false)

	for (col, (origin, title)) = enumerate(zip(["natural", "rbm", "infernal"], ["Natural", "RBM", "CM"]))		
		ax_resp = Makie.Axis(fig[1,3][1,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. / Total", title, xticks=[1e4, 5e4])
		ax_resp_norm = Makie.Axis(fig[1,3][2,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. rate", xticks=[1e4, 5e4])
		
		for (m, (primer, color)) = enumerate(zip(SamApp2024.primers20250607(), colors))
			_flag = (dms_data_primers .== primer) .&& (dms_aptamer_origin .== origin)

			println("Primer $m, $primer, $origin ", count(_flag))
			
			if any(_flag)
				read_depth = nanmean(dms_data.shape_M_depth[:, _flag, :])
				inconclusive_rate = mean(_responds_sam_inconcl_dms[_flag])
				response_rate = mean(_responds_sam_yes_dms[_flag])
				Makie.scatter!(ax_resp, read_depth, response_rate; color, markersize=15)
				Makie.scatter!(ax_resp_norm, read_depth, response_rate / (1 - inconclusive_rate); color, markersize=15)
			end
		end
		
		for ax = (ax_resp, ax_resp_norm)
			Makie.xlims!(ax, 0, 7e4)
			Makie.ylims!(ax, -0.05, 0.6)
		end
	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 09496c9f-059c-4c24-9bb9-0dbc4e95a55b
length(dms_aptamer_origin)

# ╔═╡ 25f2170b-3fd2-4b3b-be42-8c877cae748f
unique(dms_aptamer_origin)

# ╔═╡ 4b46eae3-a8de-4eaf-8207-3aa89a6594fa
unique(dms_aptamer_origin)

# ╔═╡ 0779d301-a8f7-4cf0-a98f-c7fdd18b4430
length([n for n = dms_data.aptamer_names if n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_rbm"]])

# ╔═╡ 3e99bdff-31cc-47c3-8ca5-3bf21dfcdc22
length([n for n = dms_data.aptamer_names if n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "rbm"]])

# ╔═╡ bee7d3a1-2727-4314-8113-ad770800980b
length([n for n = dms_data.aptamer_names if n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "infernal"]])

# ╔═╡ d426efcc-4a7a-433c-930e-740a37a03e1d
length([n for n = dms_data.aptamer_names if n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "Infrared"]])

# ╔═╡ c928d8df-19c2-4411-9cbf-1183fd141a96
length([n for n = dms_data.aptamer_names if startswith(n, "APSAM-S2")])

# ╔═╡ 969913ea-2c77-4fff-8846-f7bc80a77629
152+100+148

# ╔═╡ ba53fbdc-97f5-444e-945a-3ff807d1b1ee
dms_data.aptamer_names[end]

# ╔═╡ b7018bdc-8240-4a2e-bc07-5d1d1053ae3d
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], df_groups.response_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], df_groups.response_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], df_groups.response_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.response_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5,9]], df_groups.response_rate[[5,9]]; color=:teal, label="Primer 5", markersize=15)

	incl_rate = (df_groups.total_sequences - df_groups.responsives) ./ df_groups.total_sequences
	
	ax = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], incl_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], incl_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], incl_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], incl_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5,9]], incl_rate[[5,9]]; color=:teal, label="Primer 5", markersize=15)
	
	#Makie.ylims!(ax, 0.05, 0.5)
	#Makie.xlims!(ax, 0, 5e4)
	fig[1,3] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

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
# ╠═6b5e7730-952d-4533-b6a7-7caa14d82fb0
# ╠═6155ff9f-af5f-4c1a-8e6e-965f11efd285
# ╠═ed1f893e-29e0-446f-bbd4-dc5516614aa1
# ╠═1a341857-b3e4-4b7e-8172-13c399b0fca3
# ╠═e0367f59-2f02-4a09-8018-e71317695f3b
# ╠═af8c9beb-f844-4238-8967-2dbff72ac27c
# ╠═aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
# ╠═e00ec7b5-4f06-4200-9bd0-9e486e4322bc
# ╠═9763e16f-1fc6-4c79-8fed-176f847d9898
# ╠═ffa01cd2-cc36-403d-86e8-c3a7e68ac97b
# ╠═2c0519f1-367d-4fdd-a346-ad550c2e80c6
# ╠═8fc06bf7-a380-4698-8971-e96df77764b9
# ╠═df8d1d78-6cd2-40a5-a739-221ed6d03df8
# ╠═5663156d-3342-4cac-9446-1fee39923ce4
# ╠═4875da0e-9898-4498-8b29-cb60282e2cc2
# ╠═50c9a9c5-0253-4cc5-aaa5-d58472e8920c
# ╠═5f64c1d1-b429-4946-a122-a02fcb7ff1f6
# ╠═8ab280ee-8a1a-495a-ac75-dee2ca199759
# ╠═2b18233e-d2ae-485e-9ad7-17e20a26fee7
# ╠═59a5fee0-59f1-49da-9940-8827c5a23b99
# ╠═025f0c01-0e0d-4a50-b9cd-becd304a7a42
# ╠═11a25b40-e29b-4f48-8be7-2cc91feb8491
# ╠═cd5b25ff-6b98-41f3-af8d-f07665f4547f
# ╠═c69fbd19-69f4-4290-9a92-7bd5d57f79ea
# ╠═d7ba4366-918e-4b62-8756-691d07f10bf8
# ╠═8fe64552-4504-4d55-ab9e-275cdbc51293
# ╠═fe8dbdeb-5c0b-41a4-b32e-ef742fde4146
# ╠═02133d97-7e8a-4d8b-b212-c9ea18fcb9e9
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
# ╠═21f4d2dd-ad45-485f-bcc9-21a735f31c57
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
# ╠═c40c54d5-b448-41a1-b038-9e82cb804368
# ╠═154f87e7-3199-4232-980c-f6f29b1db1fe
# ╠═80306132-cc00-4b7c-a012-3890aa00793f
# ╠═0d79ddcc-f6cc-405a-8a1c-a6d3336c1476
# ╠═1d20ceef-c8ef-420d-9c31-e04db5d90a9b
# ╠═eb047b2f-f822-4e03-94bd-aabab68f9610
# ╠═529b46ff-b9cd-45c9-b202-6b2faeacb40c
# ╠═40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# ╠═6aa691a1-6e62-4077-83f6-e168810f4ec2
# ╠═8a604d45-5890-4578-8a82-01065faf1a1c
# ╠═f62fd44f-b0b8-42a7-b77f-5e5ab259293b
# ╠═500a63b1-b4c8-4624-9ddb-107808167ef8
# ╠═cc2bd6ed-6e09-453d-91ea-552aebe1ca1d
# ╠═deb8707e-3e53-4f95-ac51-e6ebf5bffe9c
# ╠═f366dabf-e8b7-4417-9a1e-99412ec6c582
# ╠═6ce57a7b-0e16-43f8-a840-d2ab4e9e03d3
# ╠═4aecd441-932b-44cb-a11b-3866f29f45d0
# ╠═1fda3f8c-4477-4039-8636-14b2ae0a5df3
# ╠═cf9a262a-42ac-40a2-b818-0a020974a6b7
# ╠═e7b53f31-137d-43cb-8d80-bb9ddf420828
# ╠═7c7c7ef4-c613-4868-ac35-57a498765ef6
# ╠═66f55016-2e10-4590-8d65-f24609f25f48
# ╠═ee9ff992-672e-46fa-87e1-c9c8adb87845
# ╠═5f372f6c-9c6d-44b1-8b29-681b4383d820
# ╠═86651e08-1f33-4fc1-8780-888714d9c86c
# ╠═a8c4c868-8d9f-4027-a265-700fbb19b895
# ╠═0e21f914-aec2-416d-a8fa-293767a3d12c
# ╠═9026821c-a1eb-4d45-96b4-28384701d1ad
# ╠═9bf2fc7d-abd7-45a7-bd0d-455b86d7c4bb
# ╠═03b6232d-295c-4e37-be8b-72ada0e1a05e
# ╠═09496c9f-059c-4c24-9bb9-0dbc4e95a55b
# ╠═25f2170b-3fd2-4b3b-be42-8c877cae748f
# ╠═86771e5d-a4b2-424b-8e72-2939148a4451
# ╠═532fb8ca-296a-4afc-b0ba-a9a40e18b03e
# ╠═b2fb16b0-93a8-450f-bde3-7d19a2ab6399
# ╠═7de7f351-2b52-43de-980c-947304970843
# ╠═150b2b5e-6f77-420e-b4e5-c1c683044c4c
# ╠═eba33d20-66dd-4cc4-b980-cd7013026c89
# ╠═4b46eae3-a8de-4eaf-8207-3aa89a6594fa
# ╠═0779d301-a8f7-4cf0-a98f-c7fdd18b4430
# ╠═3e99bdff-31cc-47c3-8ca5-3bf21dfcdc22
# ╠═bee7d3a1-2727-4314-8113-ad770800980b
# ╠═d426efcc-4a7a-433c-930e-740a37a03e1d
# ╠═c928d8df-19c2-4411-9cbf-1183fd141a96
# ╠═969913ea-2c77-4fff-8846-f7bc80a77629
# ╠═ba53fbdc-97f5-444e-945a-3ff807d1b1ee
# ╠═b7018bdc-8240-4a2e-bc07-5d1d1053ae3d

### A Pluto.jl notebook ###
# v0.20.10

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
_sites = SamApp2024.hallmark_sites_20230507 # ∪ [9, 24]

# ╔═╡ 3eba092c-11a7-44c4-b712-17b64688f0ad
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

# ╔═╡ c064d796-a9eb-4951-84a0-856e4c875e22
bps, nps, pks = SamApp2024.RF00162_sites_paired();

# ╔═╡ e0367f59-2f02-4a09-8018-e71317695f3b
md"# Load DMS data"

# ╔═╡ af8c9beb-f844-4238-8967-2dbff72ac27c
dms_df = SamApp2024.load_dms_data_sequences_table_20250303_with_aligned_sequences();

# ╔═╡ aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
dms_data = SamApp2024.load_dms_data_20250303();

# ╔═╡ e00ec7b5-4f06-4200-9bd0-9e486e4322bc
dms_data_primers = dms_df.primer[[only(findall(dms_df.name .== name)) for name = dms_data.aptamer_names]];

# ╔═╡ 2b18233e-d2ae-485e-9ad7-17e20a26fee7
dms_num_sites, dms_num_seqs, dms_num_conds = size(dms_data.shape_reactivities);

# ╔═╡ 025f0c01-0e0d-4a50-b9cd-becd304a7a42
all_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=1:dms_num_sites for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ 11a25b40-e29b-4f48-8be7-2cc91feb8491
bps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=bps for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ cd5b25ff-6b98-41f3-af8d-f07665f4547f
nps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=nps for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

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
x_mg_dms = [ismissing(dms_data.aligned_sequence[n]) ? NaN : nansum([dms_stats.shape_log_odds[i,n,2] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:400];

# ╔═╡ fe8dbdeb-5c0b-41a4-b32e-ef742fde4146
x_sam_dms = [ismissing(dms_data.aligned_sequence[n]) ? NaN : nansum([dms_stats.shape_log_odds[i,n,1] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:400];

# ╔═╡ 02133d97-7e8a-4d8b-b212-c9ea18fcb9e9
dms_AC_count = [ismissing(dms_data.aligned_sequence[n]) ? NaN : length([i for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:400];

# ╔═╡ 7959344a-1e02-45ef-9ef3-71e892e5379c
_thresh_dms = log(5)

# ╔═╡ c5ffed2c-c1f3-41a8-82c1-a5842fba5ce4
_responds_sam_yes_dms = (x_mg_dms .< -_thresh_dms) .& (x_sam_dms .> +_thresh_dms);

# ╔═╡ 010c18f0-8aea-4c91-b4d1-300117449683
_responds_sam_nop_dms = (x_mg_dms .> +_thresh_dms) .| (x_sam_dms .< -_thresh_dms);

# ╔═╡ f22df52d-cb6b-4fe8-9a16-bc1488a0122f
_responds_sam_inconcl_dms = ((!).(_responds_sam_yes_dms)) .& ((!).(_responds_sam_nop_dms));

# ╔═╡ 1844dcba-9000-4d9a-abb0-866c39976417
dms_4kqy_data = SamApp2024.load_dms_data_20250609_pdb_4kqy();

# ╔═╡ b78e2d5c-6b1a-4aea-8b72-e926aa85457d
dms_stats_4kqy_dms = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = dms_4kqy_data,
    paired_reactivities = bps_dms_reactivities,
    unpaired_reactivities = nps_dms_reactivities,
    all_reactivities = all_dms_reactivities,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ 4fa39043-b67e-4a7d-8a6e-e97a9acfcde5
x_mg_dms_4kqy = [ismissing(dms_4kqy_data.aligned_sequence[n]) ? NaN : nansum([dms_stats_4kqy_dms.shape_log_odds[i,n,2] for i = _sites if dms_4kqy_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:1];

# ╔═╡ 1a09c421-da49-48f8-ba08-3e9bb65e47cb
x_sam_dms_4kqy = [ismissing(dms_4kqy_data.aligned_sequence[n]) ? NaN : nansum([dms_stats_4kqy_dms.shape_log_odds[i,n,1] for i = _sites if dms_4kqy_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:1];

# ╔═╡ 3bea40c7-1c63-4bbd-a1fa-e28409a0c0bf
_responds_sam_yes_dms_4kqy = (x_mg_dms_4kqy .< -_thresh_dms) .& (x_sam_dms_4kqy .> +_thresh_dms)

# ╔═╡ c4fe44f9-fd26-4450-9479-afc596abe4ee
md"# Load SHAPE data (rep.0)"

# ╔═╡ f1f0126d-e212-47bc-8acc-736bb5905c0a
# load SHAPE data
#shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 5d9f93d5-007c-4d1f-9b9b-a913a7f4c00a
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20240730_with_pdb();

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

# ╔═╡ 4076ed6d-6f5f-4207-ac76-e045aa0b1c75
# All merged data, for the reactivity profiles plots
shape_data_all_merged = SamApp2024.load_shapemapper_data_pierre_demux_20240801_with_pdb_repls_merged();

# ╔═╡ 306072f8-c317-40f3-aa43-5d23d24fbcc7
conds_SAM_all_merged = map(identity, indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_allrep", "SAMAP_1M7_1SAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 3affcd4d-e17d-456f-b24d-d66236cf1c9d
conds_Mg_all_merged = map(identity, indexin(["SAMAP_1M7_noSAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ 9b74ace3-7eac-461f-b6ad-fcd4513e42f4
rep0_aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ 69de6e55-f5c6-4199-bb53-4469f1426d33
md"# SHAPE repl.0 with or without A,C"

# ╔═╡ 9e889b57-de0e-4632-94d3-030a5ee58bd0
bps_reactivities_rep0_AC = [shape_data_rep0.shape_reactivities[i,n,c] for i=bps for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ 2d7bd92e-dcb8-4987-a73d-e91ca9532a57
nps_reactivities_rep0_AC = [shape_data_rep0.shape_reactivities[i,n,c] for i=nps for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ ea8af6b3-b427-4930-b1d3-d11809c2d6c3
all_reactivities_rep0_AC = [shape_data_rep0.shape_reactivities[i,n,c] for i=1:dms_num_sites for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

# ╔═╡ a82a1fbc-a450-4a61-987d-1ae1a40a95e6
shape_stats_rep0_AC = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0_AC,
    unpaired_reactivities = nps_reactivities_rep0_AC,
    all_reactivities = all_reactivities_rep0_AC,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 72b9f9c4-64e8-4d32-836d-954843e5d860
x_mg_rep0_AC = [
	begin
		if ismissing(shape_data_rep0.aligned_sequences[n])
			NaN
		else
			nansum([shape_stats_rep0_AC.shape_log_odds[i,n,c] for i = _sites for c=conds_mg_rep0 if shape_data_rep0.aligned_sequences[n][i] ∈ ('A', 'C')])
		end
	end for n=axes(shape_stats_rep0_AC.shape_log_odds, 2)
];

# ╔═╡ cf30e3f0-b478-43ae-88d4-433a2e5c54f5
x_sam_rep0_AC = [
	begin
		if ismissing(shape_data_rep0.aligned_sequences[n])
			NaN
		else
			nansum([shape_stats_rep0_AC.shape_log_odds[i,n,c] for i=_sites for c=conds_sam_rep0 if shape_data_rep0.aligned_sequences[n][i] ∈ ('A', 'C')])
		end
	end for n=axes(shape_stats_rep0_AC.shape_log_odds, 2)
];

# ╔═╡ b1f6d289-3fb9-4d03-a4e8-f65c118f536e
_responds_sam_yes_rep0_AC = (x_mg_rep0_AC .< -_thresh_rep0) .& (x_sam_rep0_AC .> +_thresh_rep0);

# ╔═╡ eafc4efa-6af5-4b02-8e2a-490cca3f9c7c
_responds_sam_nop_rep0_AC = (x_mg_rep0_AC .> +_thresh_rep0) .| (x_sam_rep0_AC .< -_thresh_rep0);

# ╔═╡ 7d04a758-95ab-46c1-ab13-9414afc5e560
_responds_sam_inconcl_rep0_AC = ((!).(_responds_sam_yes_rep0_AC)) .& ((!).(_responds_sam_nop_rep0_AC));

# ╔═╡ d3c68b7d-5f7c-4066-85cd-ba27b23d18bc
bps_reactivities_rep0_UG = [shape_data_rep0.shape_reactivities[i,n,c] for i=bps for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('U', 'G')]

# ╔═╡ 469b2593-9b77-478c-ab8f-5c6a04ffca19
nps_reactivities_rep0_UG = [shape_data_rep0.shape_reactivities[i,n,c] for i=nps for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('U', 'G')]

# ╔═╡ 4d03ee78-d457-4b49-a35e-4a49e351a254
all_reactivities_rep0_UG = [shape_data_rep0.shape_reactivities[i,n,c] for i=1:dms_num_sites for n=nat_seqs_rep0 for c=conds_sam_rep0 if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('U', 'G')]

# ╔═╡ 2794dc19-fb40-4e53-ac86-2b3494859810
shape_stats_rep0_UG = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0_UG,
    unpaired_reactivities = nps_reactivities_rep0_UG,
    all_reactivities = all_reactivities_rep0_UG,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ e19791c3-34e9-47bb-a818-f59bc961f005
x_mg_rep0_UG = [
	begin
		if ismissing(shape_data_rep0.aligned_sequences[n])
			NaN
		else
			nansum([shape_stats_rep0_UG.shape_log_odds[i,n,c] for i = _sites for c=conds_mg_rep0 if shape_data_rep0.aligned_sequences[n][i] ∈ ('U', 'G')])
		end
	end for n=axes(shape_stats_rep0_UG.shape_log_odds, 2)
];

# ╔═╡ 6f06ce95-c96e-4566-96b6-2c02db6fd12d
x_sam_rep0_UG = [
	begin
		if ismissing(shape_data_rep0.aligned_sequences[n])
			NaN
		else
			nansum([shape_stats_rep0_UG.shape_log_odds[i,n,c] for i=_sites for c=conds_sam_rep0 if shape_data_rep0.aligned_sequences[n][i] ∈ ('U', 'G')])
		end
	end for n=axes(shape_stats_rep0_UG.shape_log_odds, 2)
];

# ╔═╡ 449ae230-8736-4590-99a9-15be03527fb6
_responds_sam_yes_rep0_UG = (x_mg_rep0_UG .< -_thresh_rep0) .& (x_sam_rep0_UG .> +_thresh_rep0);

# ╔═╡ bb1bbf19-5bce-4ce9-9bd5-b2440db7e216
_responds_sam_nop_rep0_UG = (x_mg_rep0_UG .> +_thresh_rep0) .| (x_sam_rep0_UG .< -_thresh_rep0);

# ╔═╡ d2330e77-03ee-4d49-b717-1c0b005c34d1
_responds_sam_inconcl_rep0_UG = ((!).(_responds_sam_yes_rep0_UG)) .& ((!).(_responds_sam_nop_rep0_UG));

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

# ╔═╡ 58c59835-fcce-47dd-aa1c-1d6dbb6c4616
aptamer_rbm_energies_500 = [
    RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq = shape_data_500.aligned_sequences
];

# ╔═╡ cf94281c-f97d-4966-bd3f-ddad7b8c21d0
md"# Comparison DMS vs. SHAPE Repl.0"

# ╔═╡ 10e1166d-5e07-45aa-8adc-1b92f370d26d
# indices of DMS probed sequences in Rep0 (or nothing if it is not in Rep0)
_dms_rep0_indices = indexin(dms_data.aptamer_names, shape_data_rep0.aptamer_names)

# ╔═╡ ef3554ab-8bed-40ba-9f3e-800159d899a4
length(findall(startswith("APSAMN"), dms_data.aptamer_names))

# ╔═╡ edc476ef-30aa-4fc7-bf26-041b768512e7
findall(_responds_sam_yes_dms) ∩ findall(startswith("APSAMN"), dms_data.aptamer_names) # natural only

# ╔═╡ 1d37e22b-57b1-4ffb-a31b-229eb15eb4c1
let fig = Makie.Figure()
	ex_fig_5 = "APSAMN7"
	width = 700
	height = 100
	xticks = 5:5:108

	n_ex = only(findall(shape_data_045.aptamer_names .== ex_fig_5))
	@show shape_data_045.aptamer_ids[n_ex]
	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex, only(conds_Mg_all_merged)]

	n_ex_DMS = only(findall(dms_data.aptamer_names .== ex_fig_5))
	_R_sam_DMS = dms_data.shape_reactivities[:, n_ex_DMS, 1]
	_R_mg_DMS = dms_data.shape_reactivities[:, n_ex_DMS, 2]
	seq = dms_data.aligned_sequence[n_ex_DMS]

	@assert dms_data.aligned_sequence[n_ex_DMS] == shape_data_all_merged.aligned_sequences[n_ex]

	AC_positions = [i for (i, a) = enumerate(seq) if a == 'A' || a == 'C']
	UG_positions = [i for (i, a) = enumerate(seq) if a == 'U' || a == 'G']

	_R_sam_DMS[UG_positions] .= NaN
	_R_mg_DMS[UG_positions] .= NaN

	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react. (SHAPE)", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact. (SHAPE)", xgridvisible=false, ygridvisible=false, yticks=[-4,-2,0], xtrimspine=true, ytrimspine=true
	)

	ax_react_DMS = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react. (DMS)", xgridvisible=false, ygridvisible=false, yticks=0:2:8, xtrimspine=true, ytrimspine=true, yaxisposition=:right, yticklabelcolor=:orange, ylabelcolor=:orange
	)
	ax_diff_DMS = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact. (DMS)", xgridvisible=false, ygridvisible=false, yticks=-4:0.5:1, xtrimspine=true, ytrimspine=true, yaxisposition=:right, yticklabelcolor=:orange, ylabelcolor=:orange
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, color=:purple, label="with SAM")

	Makie.scatter!(ax_react_DMS, 1:108, _R_mg_DMS; color=:gray, markersize=5)
	Makie.scatter!(ax_react_DMS, 1:108, _R_sam_DMS; color=:orange, markersize=5)

	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)
	
	Makie.barplot!(ax_diff, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -3.8one.(_sites), markersize=7, color=:black, marker=:utriangle)

	Makie.scatter!(ax_diff_DMS, 1:108, _R_sam_DMS - _R_mg_DMS; color=ifelse.(_R_sam_DMS - _R_mg_DMS .< 0, :orange, :gray), markersize=5)

	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_diff_DMS, :l, :t, :b)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidespines!(ax_react_DMS, :l, :t, :b)
	Makie.hidexdecorations!(ax_react)
	Makie.hidexdecorations!(ax_react_DMS)

	Makie.linkxaxes!(ax_react, ax_diff, ax_react_DMS, ax_diff_DMS)
	
	Makie.ylims!(ax_diff, -4, 0.7)
	Makie.ylims!(ax_diff_DMS, -1, 0.7/4)
	Makie.ylims!(ax_react, -0.5, 6)
	Makie.ylims!(ax_react_DMS, -0.25, 3)

	for ax = (ax_react, ax_diff, ax_react_DMS, ax_diff_DMS)
		Makie.xlims!(ax, 0.5, 108.5)
	end

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example2.pdf", fig)
	fig
end

# ╔═╡ eb60b8a8-905e-43b3-a3a7-b940bd563d3c
let fig = Makie.Figure()
	ex_fig_5 = "APSAMN7"
	n_ex = only(findall(dms_data.aptamer_names .== ex_fig_5))
	#n_ex = only(findall(dms_data.aptamer_names .== "APSAMN62"))t
	#n_ex = 47

	good_coverage_flag = [mean(isnan, dms_data.shape_reactivities[:, n, :]) for n = axes(dms_data.shape_reactivities, 2)] .< 0.2
	@show count(good_coverage_flag)
	admissible_examples = findall(_responds_sam_yes_dms) ∩ findall(startswith("APSAMN"), dms_data.aptamer_names) ∩ findall(good_coverage_flag) # natural only
	@show length(admissible_examples)
	
	#n_ex = rand(admissible_examples) 
	@show n_ex
	
	# _R_sam = dms_data.shape_reactivities[:, n_ex, 1]
	# _R_mg = dms_data.shape_reactivities[:, n_ex, 2]
	# seq = dms_data.aligned_sequence[n_ex]

	_R_sam = dms_4kqy_data.shape_reactivities[:, 1, 1]
	_R_mg = dms_4kqy_data.shape_reactivities[:, 1, 2]
	seq = dms_4kqy_data.aligned_sequence[1]

	width = 700
	height = 100
	xticks = 5:5:108

	AC_positions = Int[]
	UG_positions = Int[]
	for i = 1:108
		if seq[i] == 'A' || seq[i] == 'C'
			push!(AC_positions, i)
		elseif seq[i] == 'U' || seq[i] == 'G'
			push!(UG_positions, i)
		end
	end

	_R_sam[UG_positions] .= NaN
	_R_mg[UG_positions] .= NaN
	
	ax_react = Makie.Axis(
		fig[1,1]; valign=:bottom, width, height, xticks, ylabel="react.", xgridvisible=false, ygridvisible=false, yticks=0:1:8, xtrimspine=true, ytrimspine=true
	)
	ax_diff = Makie.Axis(
		fig[2,1]; valign=:bottom, width, height, xticks, xlabel="site", ylabel="Δreact.", xgridvisible=false, ygridvisible=false, 
		yticks=-2:0.5:0, xtrimspine=true, ytrimspine=true
	)

	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react, x0, xf; color=(color, alpha))
		Makie.vspan!(ax_diff, x0, xf; color=(color, alpha))
	end
	
	Makie.stairs!(ax_react, 1:108, _R_mg; step=:center, linewidth=2, color=:gray, label="no SAM")
	Makie.stairs!(ax_react, 1:108, _R_sam; step=:center, linewidth=2, color=:purple, label="with SAM")
	#Makie.axislegend(ax_react, position=(0.0, -13), framevisible=false)
	#Makie.hidespines!(ax_react_1, :t, :r, :b)
	#Makie.hidexdecorations!(ax_react_1)

	#ΔR = replace(_R_sam - _R_mg, NaN => 0)
	ΔR = _R_sam - _R_mg
	
	Makie.barplot!(ax_diff, 1:108, ΔR; color=ifelse.(ΔR .< 0, :green, :gray))
	Makie.scatter!(ax_diff, _sites, -0.75one.(_sites), markersize=7, color=:black, marker=:utriangle)

	# Makie.scatter!(ax_react, AC_positions, 1.5one.(AC_positions), markersize=7, color=:blue, marker=:utriangle)
	# Makie.scatter!(ax_react, UG_positions, 1.5one.(UG_positions), markersize=7, color=:red, marker=:utriangle)
	
	Makie.xlims!(ax_diff, 0, 109)
	
	Makie.hidespines!(ax_diff, :r, :t)
	Makie.hidespines!(ax_react, :r, :t, :b)
	Makie.hidexdecorations!(ax_react)
	#Makie.scatter!(ax_diff_1, _sites, -0.2one.(_sites), color=:blue, markersize=5)

	Makie.linkxaxes!(ax_react, ax_diff)
	Makie.ylims!(ax_diff, -0.8, 0.2)
	Makie.ylims!(ax_react, -0.5, 1.5)
	
	Makie.xlims!(ax_react, 0.5, 108.5)
	Makie.xlims!(ax_diff,  0.5, 108.5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Fig5new_SHAPE_example.pdf", fig)
	fig
end

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

# ╔═╡ bba2554c-f415-473a-b442-a3ee19adb431
[ # RBM low energy only
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_rbm")[filter(!isnothing, _dms_rep0_indices)] .&& 
			((!ismissing).(rep0_aptamer_rbm_energies) .&& -rep0_aptamer_rbm_energies .> 300)[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ f1f4f486-8867-4c24-bee1-08ac1c927c96
findall(ismissing, (-rep0_aptamer_rbm_energies .> 300)[filter(!isnothing, _dms_rep0_indices)])

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

# ╔═╡ 40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# per primers
[
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (dms_data_primers .== SamApp2024.primers20250607()[5])[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ e727a7a7-5c40-43d0-9686-058483705c0e
md"## Comparison DMS+SHAPE on repl.0"

# ╔═╡ eca43047-bfeb-4e0a-b272-aa87167dfbbf
x_mg_dms_plus_shape_rep0 = x_mg_dms[findall(!isnothing, _dms_rep0_indices)] + x_mg_rep0[filter(!isnothing, _dms_rep0_indices)];

# ╔═╡ 882e33ac-2df0-4c02-9379-7896ca64aeea
x_sam_dms_plus_shape_rep0 = x_sam_dms[findall(!isnothing, _dms_rep0_indices)] + x_sam_rep0[filter(!isnothing, _dms_rep0_indices)];

# ╔═╡ a64dbc0f-b9c8-4b59-a33a-92c6d732b792
_responds_sam_yes_dms_plus_shape_rep0 = (x_mg_dms_plus_shape_rep0 .< -_thresh_rep0) .& (x_sam_dms_plus_shape_rep0 .> +_thresh_rep0);

# ╔═╡ b0998da2-0b6f-469e-a2dc-7e73e809f1c2
_responds_sam_nop_dms_plus_shape_rep0 = (x_mg_dms_plus_shape_rep0 .> +_thresh_rep0) .| (x_sam_dms_plus_shape_rep0 .< -_thresh_rep0);

# ╔═╡ 55708acd-bbfd-4b2b-b533-d5b1d86ecfbb
_responds_sam_inconcl_dms_plus_shape_rep0 = ((!).(_responds_sam_yes_dms_plus_shape_rep0)) .& ((!).(_responds_sam_nop_dms_plus_shape_rep0));

# ╔═╡ 07902e15-3f55-4355-a154-bf712aad41fe
[ # natural
	begin
		count(dms_plus_shape .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms_plus_shape = (_responds_sam_yes_dms_plus_shape_rep0, _responds_sam_nop_dms_plus_shape_rep0, _responds_sam_inconcl_dms_plus_shape_rep0)
]

# ╔═╡ fcba9ea4-e623-4ebe-ba13-39bdf48984db
[ # RBM only
	begin
		count(dms_plus_shape .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_rbm")[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms_plus_shape = (_responds_sam_yes_dms_plus_shape_rep0, _responds_sam_nop_dms_plus_shape_rep0, _responds_sam_inconcl_dms_plus_shape_rep0)
]

# ╔═╡ a7e32428-c6e2-4bf1-bb99-73da2831e202
[ # RBM low energy only
	begin
		count(dms_plus_shape .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_rbm")[filter(!isnothing, _dms_rep0_indices)] .&& 
			((!ismissing).(rep0_aptamer_rbm_energies) .&& -rep0_aptamer_rbm_energies .> 300)[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms_plus_shape = (_responds_sam_yes_dms_plus_shape_rep0, _responds_sam_nop_dms_plus_shape_rep0, _responds_sam_inconcl_dms_plus_shape_rep0)
]

# ╔═╡ c8149c54-d95a-4e2c-94ee-dc823219e381
[ # CM
	begin
		count(dms_plus_shape .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .== "RF00162_syn_inf")[filter(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0, _responds_sam_nop_rep0, _responds_sam_inconcl_rep0), dms_plus_shape = (_responds_sam_yes_dms_plus_shape_rep0, _responds_sam_nop_dms_plus_shape_rep0, _responds_sam_inconcl_dms_plus_shape_rep0)
]

# ╔═╡ 62e04c51-bfbd-4f3b-a057-ca0d75481ca8
md"## Comparison DMS vs SHAPE with A,C only"

# ╔═╡ c86f4278-c152-4b5b-9660-3c3134049e46
[ # natural (SHAPE A,C)
	begin
		count(
			dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)]
		)
	end for rep0 = (_responds_sam_yes_rep0_AC, _responds_sam_nop_rep0_AC, _responds_sam_inconcl_rep0_AC), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 34159350-a262-4121-b600-98d3e0a2f1b2
[ # natural (SHAPE A,C)
	begin
		count(
			dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)].&&
			_responds_sam_yes_rep0[filter(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)]
		)
	end for rep0 = (_responds_sam_yes_rep0_AC, _responds_sam_nop_rep0_AC, _responds_sam_inconcl_rep0_AC), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ fe0ed922-c2b8-4f25-bd35-211c4f70e02e
[ # natural (SHAPE U,G)
	begin
		count(
			dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)]
		)
	end for rep0 = (_responds_sam_yes_rep0_UG, _responds_sam_nop_rep0_UG, _responds_sam_inconcl_rep0_UG), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 5bff2e5e-4aab-45df-a088-0d97cccc4def
[ # natural (SHAPE A,C)
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)] .&&
			_responds_sam_yes_rep0[filter(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_AC, _responds_sam_nop_rep0_AC, _responds_sam_inconcl_rep0_AC), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ a8090cba-8d47-4371-95b8-3aacdac29594
[ # natural (SHAPE U, G)
	begin
		count(dms[findall(!isnothing, _dms_rep0_indices)] .&& rep0[filter(!isnothing, _dms_rep0_indices)] .&& (shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)] .&&_responds_sam_yes_rep0[filter(!isnothing, _dms_rep0_indices)] .&& 		_responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)])
	end for rep0 = (_responds_sam_yes_rep0_UG, _responds_sam_nop_rep0_UG, _responds_sam_inconcl_rep0_UG), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 75902f05-30c0-4e0e-b696-eb1e9e96cd2e
[ # natural (SHAPE A,C vs. U,G)
	begin
		count(rep0_AC[filter(!isnothing, _dms_rep0_indices)] .&& rep0_UG[filter(!isnothing, _dms_rep0_indices)] .&& 
			(shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)] .&& 
				_responds_sam_yes_rep0[filter(!isnothing, _dms_rep0_indices)] .&& _responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)])
	end for rep0_AC = (_responds_sam_yes_rep0_AC, _responds_sam_nop_rep0_AC, _responds_sam_inconcl_rep0_AC), rep0_UG = (_responds_sam_yes_rep0_UG, _responds_sam_nop_rep0_UG, _responds_sam_inconcl_rep0_UG)
]

# ╔═╡ f696a001-9e73-4cef-8940-760718dbb4c7
[ # natural (SHAPE A,C vs. U,G)
	begin
		count(rep0_AC[filter(!isnothing, _dms_rep0_indices)] .&& rep0_UG[filter(!isnothing, _dms_rep0_indices)] .&& 
			(shape_data_045.aptamer_origin .∈ Ref(["RF00162_full30", "RF00162_seed70"]))[filter(!isnothing, _dms_rep0_indices)] #.&& 
				#_responds_sam_yes_rep0[filter(!isnothing, _dms_rep0_indices)] #.&& _responds_sam_nop_dms[findall(!isnothing, _dms_rep0_indices)]
	)
	end for rep0_AC = (_responds_sam_yes_rep0_AC, _responds_sam_nop_rep0_AC, _responds_sam_inconcl_rep0_AC), rep0_UG = (_responds_sam_yes_rep0_UG, _responds_sam_nop_rep0_UG, _responds_sam_inconcl_rep0_UG)
]

# ╔═╡ 6aa691a1-6e62-4077-83f6-e168810f4ec2
md"# Comparison DMS vs. SHAPE 500"

# ╔═╡ 80eee207-09e6-45fb-bdac-ce473a02cc36
dms_data.aligned_sequence[253:end] ⊆ shape_data_500.aligned_sequences

# ╔═╡ b60a7dd7-9bdf-461a-b73b-7ba2ed3dc797
shape_data_500.aptamer_origin |> unique

# ╔═╡ 8a604d45-5890-4578-8a82-01065faf1a1c
# indices of DMS probed sequences in 500 (or nothing if it is not in 500)
_dms_500_indices = indexin(dms_data.aligned_sequence, map(string, shape_data_500.aligned_sequences[1:450]))

# ╔═╡ 1e02902f-da1e-4ff2-a9dd-4561ccc2bb53
count(shape_data_500.aptamer_origin[filter(!isnothing, _dms_500_indices)] .== "infernal")

# ╔═╡ f62fd44f-b0b8-42a7-b77f-5e5ab259293b
string.(shape_data_500.aligned_sequences[filter(!isnothing, _dms_500_indices)]) == dms_data.aligned_sequence[findall(!isnothing, _dms_500_indices)]

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
		count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (shape_data_500.aptamer_origin .== "rbm")[filter(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ 1fda3f8c-4477-4039-8636-14b2ae0a5df3
# per primer
[
	count(dms[findall(!isnothing, _dms_500_indices)] .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (dms_data_primers .== SamApp2024.primers20250607()[1])[findall(!isnothing, _dms_500_indices)])
	for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms, _responds_sam_nop_dms, _responds_sam_inconcl_dms)
]

# ╔═╡ afbf76df-ba31-44bc-945c-92efd2f7b9f9
md"## Comparison DMS+SHAPE on 500"

# ╔═╡ c479a6ff-8570-4e00-8e07-e04f87f7a953
x_mg_dms_plus_shape_500 = x_mg_dms[findall(!isnothing, _dms_500_indices)] + x_mg_500[filter(!isnothing, _dms_500_indices)];

# ╔═╡ e40dba2b-bb9a-4436-8554-14184873c223
x_sam_dms_plus_shape_500 = x_sam_dms[findall(!isnothing, _dms_500_indices)] + x_sam_500[filter(!isnothing, _dms_500_indices)];

# ╔═╡ df2043f0-38f4-4a8d-8bc7-454c981ef6ea
_responds_sam_yes_dms_plus_shape_500 = (x_mg_dms_plus_shape_500 .< -_thresh_rep0) .& (x_sam_dms_plus_shape_500 .> +_thresh_rep0);

# ╔═╡ eea8187c-bdd9-405f-bf0f-870300695133
_responds_sam_nop_dms_plus_shape_500 = (x_mg_dms_plus_shape_500 .> +_thresh_rep0) .| (x_sam_dms_plus_shape_500 .< -_thresh_rep0);

# ╔═╡ 2cc2ec0b-df1f-4d30-806b-f7c2792f6767
_responds_sam_inconcl_dms_plus_shape_500 = ((!).(_responds_sam_yes_dms_plus_shape_500)) .& ((!).(_responds_sam_nop_dms_plus_shape_500));

# ╔═╡ 2764d67d-4e27-48be-8163-c2e98f942cb5
[
	count(dms .&& shape500[filter(!isnothing, _dms_500_indices)])
	for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms_plus_shape_500, _responds_sam_nop_dms_plus_shape_500, _responds_sam_inconcl_dms_plus_shape_500)
]

# ╔═╡ bfcb3f30-83ec-4ab2-9ffa-7dd91333f33e
[ # RBM
	begin
		count(dms .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (shape_data_500.aptamer_origin .== "infernal")[filter(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms_plus_shape_500, _responds_sam_nop_dms_plus_shape_500, _responds_sam_inconcl_dms_plus_shape_500)
]

# ╔═╡ 292c6c81-bebd-4645-bcad-8ffd26338e21
[ # CM
	begin
		count(dms .&& shape500[filter(!isnothing, _dms_500_indices)] .&& (shape_data_500.aptamer_origin .== "infernal")[filter(!isnothing, _dms_500_indices)])
	end for shape500 = (_responds_sam_yes_500, _responds_sam_nop_500, _responds_sam_inconcl_500), dms = (_responds_sam_yes_dms_plus_shape_500, _responds_sam_nop_dms_plus_shape_500, _responds_sam_inconcl_dms_plus_shape_500)
]

# ╔═╡ 354c9b92-327c-4080-8f4a-29f38df0aefd
shape_data_500.aptamer_origin |> unique

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
# ╠═3eba092c-11a7-44c4-b712-17b64688f0ad
# ╠═c064d796-a9eb-4951-84a0-856e4c875e22
# ╠═e0367f59-2f02-4a09-8018-e71317695f3b
# ╠═af8c9beb-f844-4238-8967-2dbff72ac27c
# ╠═aae7a9e7-cf14-4c09-a6fb-93d6e1e19b3d
# ╠═e00ec7b5-4f06-4200-9bd0-9e486e4322bc
# ╠═2b18233e-d2ae-485e-9ad7-17e20a26fee7
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
# ╠═1844dcba-9000-4d9a-abb0-866c39976417
# ╠═b78e2d5c-6b1a-4aea-8b72-e926aa85457d
# ╠═4fa39043-b67e-4a7d-8a6e-e97a9acfcde5
# ╠═1a09c421-da49-48f8-ba08-3e9bb65e47cb
# ╠═3bea40c7-1c63-4bbd-a1fa-e28409a0c0bf
# ╠═c4fe44f9-fd26-4450-9479-afc596abe4ee
# ╠═f1f0126d-e212-47bc-8acc-736bb5905c0a
# ╠═5d9f93d5-007c-4d1f-9b9b-a913a7f4c00a
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
# ╠═4076ed6d-6f5f-4207-ac76-e045aa0b1c75
# ╠═306072f8-c317-40f3-aa43-5d23d24fbcc7
# ╠═3affcd4d-e17d-456f-b24d-d66236cf1c9d
# ╠═9b74ace3-7eac-461f-b6ad-fcd4513e42f4
# ╠═69de6e55-f5c6-4199-bb53-4469f1426d33
# ╠═9e889b57-de0e-4632-94d3-030a5ee58bd0
# ╠═2d7bd92e-dcb8-4987-a73d-e91ca9532a57
# ╠═ea8af6b3-b427-4930-b1d3-d11809c2d6c3
# ╠═a82a1fbc-a450-4a61-987d-1ae1a40a95e6
# ╠═72b9f9c4-64e8-4d32-836d-954843e5d860
# ╠═cf30e3f0-b478-43ae-88d4-433a2e5c54f5
# ╠═b1f6d289-3fb9-4d03-a4e8-f65c118f536e
# ╠═eafc4efa-6af5-4b02-8e2a-490cca3f9c7c
# ╠═7d04a758-95ab-46c1-ab13-9414afc5e560
# ╠═d3c68b7d-5f7c-4066-85cd-ba27b23d18bc
# ╠═469b2593-9b77-478c-ab8f-5c6a04ffca19
# ╠═4d03ee78-d457-4b49-a35e-4a49e351a254
# ╠═2794dc19-fb40-4e53-ac86-2b3494859810
# ╠═e19791c3-34e9-47bb-a818-f59bc961f005
# ╠═6f06ce95-c96e-4566-96b6-2c02db6fd12d
# ╠═449ae230-8736-4590-99a9-15be03527fb6
# ╠═bb1bbf19-5bce-4ce9-9bd5-b2440db7e216
# ╠═d2330e77-03ee-4d49-b717-1c0b005c34d1
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
# ╠═58c59835-fcce-47dd-aa1c-1d6dbb6c4616
# ╠═cf94281c-f97d-4966-bd3f-ddad7b8c21d0
# ╠═10e1166d-5e07-45aa-8adc-1b92f370d26d
# ╠═ef3554ab-8bed-40ba-9f3e-800159d899a4
# ╠═edc476ef-30aa-4fc7-bf26-041b768512e7
# ╠═1d37e22b-57b1-4ffb-a31b-229eb15eb4c1
# ╠═eb60b8a8-905e-43b3-a3a7-b940bd563d3c
# ╠═154f87e7-3199-4232-980c-f6f29b1db1fe
# ╠═80306132-cc00-4b7c-a012-3890aa00793f
# ╠═bba2554c-f415-473a-b442-a3ee19adb431
# ╠═f1f4f486-8867-4c24-bee1-08ac1c927c96
# ╠═0d79ddcc-f6cc-405a-8a1c-a6d3336c1476
# ╠═1d20ceef-c8ef-420d-9c31-e04db5d90a9b
# ╠═40c13e35-4c52-4051-8f85-fdb1ca99d3b1
# ╠═e727a7a7-5c40-43d0-9686-058483705c0e
# ╠═eca43047-bfeb-4e0a-b272-aa87167dfbbf
# ╠═882e33ac-2df0-4c02-9379-7896ca64aeea
# ╠═a64dbc0f-b9c8-4b59-a33a-92c6d732b792
# ╠═b0998da2-0b6f-469e-a2dc-7e73e809f1c2
# ╠═55708acd-bbfd-4b2b-b533-d5b1d86ecfbb
# ╠═07902e15-3f55-4355-a154-bf712aad41fe
# ╠═fcba9ea4-e623-4ebe-ba13-39bdf48984db
# ╠═a7e32428-c6e2-4bf1-bb99-73da2831e202
# ╠═c8149c54-d95a-4e2c-94ee-dc823219e381
# ╠═62e04c51-bfbd-4f3b-a057-ca0d75481ca8
# ╠═c86f4278-c152-4b5b-9660-3c3134049e46
# ╠═34159350-a262-4121-b600-98d3e0a2f1b2
# ╠═fe0ed922-c2b8-4f25-bd35-211c4f70e02e
# ╠═5bff2e5e-4aab-45df-a088-0d97cccc4def
# ╠═a8090cba-8d47-4371-95b8-3aacdac29594
# ╠═75902f05-30c0-4e0e-b696-eb1e9e96cd2e
# ╠═f696a001-9e73-4cef-8940-760718dbb4c7
# ╠═6aa691a1-6e62-4077-83f6-e168810f4ec2
# ╠═80eee207-09e6-45fb-bdac-ce473a02cc36
# ╠═b60a7dd7-9bdf-461a-b73b-7ba2ed3dc797
# ╠═1e02902f-da1e-4ff2-a9dd-4561ccc2bb53
# ╠═8a604d45-5890-4578-8a82-01065faf1a1c
# ╠═f62fd44f-b0b8-42a7-b77f-5e5ab259293b
# ╠═cc2bd6ed-6e09-453d-91ea-552aebe1ca1d
# ╠═deb8707e-3e53-4f95-ac51-e6ebf5bffe9c
# ╠═f366dabf-e8b7-4417-9a1e-99412ec6c582
# ╠═6ce57a7b-0e16-43f8-a840-d2ab4e9e03d3
# ╠═4aecd441-932b-44cb-a11b-3866f29f45d0
# ╠═1fda3f8c-4477-4039-8636-14b2ae0a5df3
# ╠═afbf76df-ba31-44bc-945c-92efd2f7b9f9
# ╠═c479a6ff-8570-4e00-8e07-e04f87f7a953
# ╠═e40dba2b-bb9a-4436-8554-14184873c223
# ╠═df2043f0-38f4-4a8d-8bc7-454c981ef6ea
# ╠═eea8187c-bdd9-405f-bf0f-870300695133
# ╠═2cc2ec0b-df1f-4d30-806b-f7c2792f6767
# ╠═2764d67d-4e27-48be-8163-c2e98f942cb5
# ╠═bfcb3f30-83ec-4ab2-9ffa-7dd91333f33e
# ╠═292c6c81-bebd-4645-bcad-8ffd26338e21
# ╠═354c9b92-327c-4080-8f4a-29f38df0aefd
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

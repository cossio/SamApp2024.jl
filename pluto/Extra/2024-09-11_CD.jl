### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ e7e0d865-a7d4-4aa8-a88e-a2e8dec2a6d4
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ fdf0ec94-8a2b-4bc1-a4dc-cfdfb2d4dc27
using BioSequences: LongRNA

# ╔═╡ e5cd7b7b-a53d-4e56-9522-fc495aeba660
using Distributions: Gamma

# ╔═╡ 4301a5ff-2d3a-42eb-9937-9d623612ec6c
using Makie: @L_str

# ╔═╡ fd1fd35c-aaa1-4988-ae03-cc512f60a7fa
using NaNStatistics: nanmean

# ╔═╡ d9ba83e0-ed7f-4d62-b436-a429de5c09fa
using NaNStatistics: nanstd

# ╔═╡ 2f00c339-fbad-4294-9e2e-4bbb1c51148e
using NaNStatistics: nansum

# ╔═╡ 9886be7b-68c2-41b3-964f-f729e39c8f83
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 44d319cc-2596-426f-9471-895803f4b66d
using Statistics: cor

# ╔═╡ 83bdb153-d516-4070-8e6a-9bd317ae2546
using Statistics: mean

# ╔═╡ a8209f85-fb44-4883-9dc5-1bb6411c599d
using DataFrames: DataFrame

# ╔═╡ e532bc4b-9876-4f4a-8007-f5de0207b4de
md"# Imports"

# ╔═╡ ab67353e-874a-4ce7-b0f7-945cbdf997a8
import PlutoUI

# ╔═╡ 5870f4e0-5683-4f06-bc9b-f91d2222006b
import CairoMakie

# ╔═╡ 64738c7e-24ab-4c77-8c29-8927e38dd52d
import Makie

# ╔═╡ 2d3d8134-0fa0-4d2f-a2da-21541b5e6819
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 5ba2b0d9-9c7f-4351-a0a1-53c1908b1d61
import SamApp2024

# ╔═╡ add115e1-e5ee-4e6b-87f0-5f91ba4678f6
import StatsBase

# ╔═╡ 0a3539fd-8557-44ca-b5ab-2aaa814225f3
import CSV

# ╔═╡ d0d45ec0-295a-482b-bd91-f822555a52bf
PlutoUI.TableOfContents()

# ╔═╡ f47c5268-b535-4f21-b63f-df0f2b4375be
md"# General data"

# ╔═╡ 60f6c3a7-f561-4646-931e-789b558d6851
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 5a85be27-eae3-4ace-834a-f38125d494fb
ss = SamApp2024.RF00162_sites_annotated_secondary_structure()

# ╔═╡ ef0da2da-847c-4b58-ade6-9eb0193be17c
# Hallmark sites
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ dee39b66-66e8-4929-9976-b28ca6fd06b4
# signifcance threshold for protection scores
_thresh = log(5)

# ╔═╡ 5fc7ed72-b26c-4bbc-8331-f527e538ef09
md"# Load data (Repl.0)"

# ╔═╡ 9b57bad1-99f9-42ed-92dc-a7850287c641
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 97572eb9-6064-4168-a0c9-896ce49ab87f
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 32449e0d-b43e-4d54-b41a-282c3778ea17
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ f0d5273b-29aa-4a03-88b3-2d314999596f
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ eb0270e1-65c8-4c5e-a3e4-5a52b600830e
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 06a82250-c089-4cad-a6ef-7f0009a52a3c
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ a533b2de-9fc3-4f95-8487-c3a3de831a64
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 4f2d2af7-0314-46f4-9500-f304d851cd79
inf_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 00065bc9-1fd3-4023-9d67-74ad37024d23
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 5f743baa-eb33-4d7c-87a4-7862c93ba150
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 8bd013f1-6e3b-428b-9c38-67c08b2e90b9
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ 65d19441-0cb9-4496-a883-6e5cc2dc1fe8
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ b985e5a1-5677-46af-a69c-c44e7ced82ad
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 48fe30e6-608e-49ac-88dc-b9adbcac0245
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ f9dfec41-e434-4c91-ade6-ea558013e486
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ e799ec52-8ea9-4904-920f-cb301aeb53e5
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3));

# ╔═╡ 2be6a398-e38c-4d6d-b09d-38726c172c01
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3));

# ╔═╡ 7c726486-616c-4c97-90e1-d34b827d5969
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 42b2e97e-547e-4cc4-9eb2-09d544c8270c
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ b0ec2e5e-7d50-4bb9-9893-4052301525c3
aptamer_rbm_energies_rep0 = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ 4f4c29ba-664b-4d2d-b7af-0b8fa751aeb7
md"# Load data (500 seqs)"

# ╔═╡ 63568092-6bfa-493a-9038-87210b09dddb
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ bd65f0ad-7d8e-4355-a942-9c69ef77150c
conds_sam_500, conds_mg_500, conds_30C_500 = [1,2], [4], [6];

# ╔═╡ b2c55a3c-0a19-4325-bd1a-b510af42e2f3
bps_reactivities_500 = shape_data_500.shape_reactivities[bps, :, conds_sam_500];

# ╔═╡ bbe19c79-93d1-4008-82d4-06bce3ca2592
nps_reactivities_500 = shape_data_500.shape_reactivities[nps, :, conds_sam_500];

# ╔═╡ d2cd622b-9289-4c6e-9b68-d2d82d8738a8
all_reactivities_500 = shape_data_500.shape_reactivities[:, :, conds_sam_500];

# ╔═╡ 3376b49e-4ca1-4548-921c-82eb681d5345
shape_stats_500 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = bps_reactivities_500,
    unpaired_reactivities = nps_reactivities_500,
    all_reactivities = all_reactivities_500,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ cf3d0b96-8d28-4ab1-905c-658be880e84f
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :,  conds_mg_500]; dim=(1,3));

# ╔═╡ 334237cc-11fa-48f0-8bc1-258128adc5f4
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ ad74cffd-0574-4e58-af8b-2cee2e4082d1
_responds_sam_yes_500 = (x_mg_500 .< -_thresh) .& (x_sam_500 .> +_thresh);

# ╔═╡ 8f34d835-271b-4195-8dcf-cd0c2cc3ae53
_responds_sam_nop_500 = (x_mg_500 .> +_thresh) .| (x_sam_500 .< -_thresh);

# ╔═╡ b118f1fb-1aac-4722-b5b7-a2bda4eaa7df
aptamer_rbm_energies_500 = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_500.aligned_sequences
];

# ╔═╡ 8cc027be-e750-4277-b3ac-38636ac5259b
md"# Protection scores of probed sequences"

# ╔═╡ 2bf7e337-900c-4195-b5cb-b1642551bec4
aptamers_rep0 = [only(findall(shape_data_rep0.aptamer_names .== ap)) for ap = ["APSAMS10", "APSAMS25", "APSAMN75", "APSAMN172"]]

# ╔═╡ c95ccad6-56cc-45ae-a571-15012cd09064
aptamers_500 = [143, 282]

# ╔═╡ 1f823392-02e0-4712-9cf5-0407b93556bb
md"## Protection scores of P1"

# ╔═╡ 1eb6ed79-a71d-409b-bbcb-6a080e60afd0
[
	nansum(shape_stats_rep0.shape_log_odds[ss.p1, aptamers_rep0,  conds_mg_rep0]; dim=(1,3)) ;;
	nansum(shape_stats_rep0.shape_log_odds[ss.p1, aptamers_rep0,  conds_sam_rep0]; dim=(1,3))
]

# ╔═╡ ae0a8987-a7ad-4579-9bab-438fd6a0e35b
[
	nansum(shape_stats_500.shape_log_odds[ss.p1, aptamers_500,  conds_mg_500]; dim=(1,3)) ;;
	nansum(shape_stats_500.shape_log_odds[ss.p1, aptamers_500,  conds_sam_500]; dim=(1,3))
]

# ╔═╡ 94b9eae6-f0db-48d9-a469-51abae74eda2
md"## Protection scores of Hallmark sites"

# ╔═╡ ebfc5dd0-f4e7-430d-a263-e417e46f1573
[
	nansum(shape_stats_rep0.shape_log_odds[_sites, aptamers_rep0,  conds_mg_rep0]; dim=(1,3)) ;;
	nansum(shape_stats_rep0.shape_log_odds[_sites, aptamers_rep0,  conds_sam_rep0]; dim=(1,3))
]

# ╔═╡ 5b089195-a7af-41c1-ab28-ec94d23e96a3
[
	nansum(shape_stats_500.shape_log_odds[_sites, aptamers_500,  conds_mg_500]; dim=(1,3)) ;;
	nansum(shape_stats_500.shape_log_odds[_sites, aptamers_500,  conds_sam_500]; dim=(1,3))
]

# ╔═╡ Cell order:
# ╠═e532bc4b-9876-4f4a-8007-f5de0207b4de
# ╠═e7e0d865-a7d4-4aa8-a88e-a2e8dec2a6d4
# ╠═ab67353e-874a-4ce7-b0f7-945cbdf997a8
# ╠═5870f4e0-5683-4f06-bc9b-f91d2222006b
# ╠═64738c7e-24ab-4c77-8c29-8927e38dd52d
# ╠═2d3d8134-0fa0-4d2f-a2da-21541b5e6819
# ╠═5ba2b0d9-9c7f-4351-a0a1-53c1908b1d61
# ╠═add115e1-e5ee-4e6b-87f0-5f91ba4678f6
# ╠═0a3539fd-8557-44ca-b5ab-2aaa814225f3
# ╠═fdf0ec94-8a2b-4bc1-a4dc-cfdfb2d4dc27
# ╠═e5cd7b7b-a53d-4e56-9522-fc495aeba660
# ╠═4301a5ff-2d3a-42eb-9937-9d623612ec6c
# ╠═fd1fd35c-aaa1-4988-ae03-cc512f60a7fa
# ╠═d9ba83e0-ed7f-4d62-b436-a429de5c09fa
# ╠═2f00c339-fbad-4294-9e2e-4bbb1c51148e
# ╠═9886be7b-68c2-41b3-964f-f729e39c8f83
# ╠═44d319cc-2596-426f-9471-895803f4b66d
# ╠═83bdb153-d516-4070-8e6a-9bd317ae2546
# ╠═a8209f85-fb44-4883-9dc5-1bb6411c599d
# ╠═d0d45ec0-295a-482b-bd91-f822555a52bf
# ╠═f47c5268-b535-4f21-b63f-df0f2b4375be
# ╠═60f6c3a7-f561-4646-931e-789b558d6851
# ╠═5a85be27-eae3-4ace-834a-f38125d494fb
# ╠═ef0da2da-847c-4b58-ade6-9eb0193be17c
# ╠═dee39b66-66e8-4929-9976-b28ca6fd06b4
# ╠═5fc7ed72-b26c-4bbc-8331-f527e538ef09
# ╠═9b57bad1-99f9-42ed-92dc-a7850287c641
# ╠═97572eb9-6064-4168-a0c9-896ce49ab87f
# ╠═32449e0d-b43e-4d54-b41a-282c3778ea17
# ╠═f0d5273b-29aa-4a03-88b3-2d314999596f
# ╠═eb0270e1-65c8-4c5e-a3e4-5a52b600830e
# ╠═06a82250-c089-4cad-a6ef-7f0009a52a3c
# ╠═a533b2de-9fc3-4f95-8487-c3a3de831a64
# ╠═4f2d2af7-0314-46f4-9500-f304d851cd79
# ╠═00065bc9-1fd3-4023-9d67-74ad37024d23
# ╠═5f743baa-eb33-4d7c-87a4-7862c93ba150
# ╠═8bd013f1-6e3b-428b-9c38-67c08b2e90b9
# ╠═65d19441-0cb9-4496-a883-6e5cc2dc1fe8
# ╠═b985e5a1-5677-46af-a69c-c44e7ced82ad
# ╠═48fe30e6-608e-49ac-88dc-b9adbcac0245
# ╠═f9dfec41-e434-4c91-ade6-ea558013e486
# ╠═e799ec52-8ea9-4904-920f-cb301aeb53e5
# ╠═2be6a398-e38c-4d6d-b09d-38726c172c01
# ╠═7c726486-616c-4c97-90e1-d34b827d5969
# ╠═42b2e97e-547e-4cc4-9eb2-09d544c8270c
# ╠═b0ec2e5e-7d50-4bb9-9893-4052301525c3
# ╠═4f4c29ba-664b-4d2d-b7af-0b8fa751aeb7
# ╠═63568092-6bfa-493a-9038-87210b09dddb
# ╠═bd65f0ad-7d8e-4355-a942-9c69ef77150c
# ╠═b2c55a3c-0a19-4325-bd1a-b510af42e2f3
# ╠═bbe19c79-93d1-4008-82d4-06bce3ca2592
# ╠═d2cd622b-9289-4c6e-9b68-d2d82d8738a8
# ╠═3376b49e-4ca1-4548-921c-82eb681d5345
# ╠═cf3d0b96-8d28-4ab1-905c-658be880e84f
# ╠═334237cc-11fa-48f0-8bc1-258128adc5f4
# ╠═ad74cffd-0574-4e58-af8b-2cee2e4082d1
# ╠═8f34d835-271b-4195-8dcf-cd0c2cc3ae53
# ╠═b118f1fb-1aac-4722-b5b7-a2bda4eaa7df
# ╠═8cc027be-e750-4277-b3ac-38636ac5259b
# ╠═2bf7e337-900c-4195-b5cb-b1642551bec4
# ╠═c95ccad6-56cc-45ae-a571-15012cd09064
# ╠═1f823392-02e0-4712-9cf5-0407b93556bb
# ╠═1eb6ed79-a71d-409b-bbcb-6a080e60afd0
# ╠═ae0a8987-a7ad-4579-9bab-438fd6a0e35b
# ╠═94b9eae6-f0db-48d9-a469-51abae74eda2
# ╠═ebfc5dd0-f4e7-430d-a263-e417e46f1573
# ╠═5b089195-a7af-41c1-ab28-ec94d23e96a3

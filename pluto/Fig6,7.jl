### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 77855f93-2b64-45bc-a307-df7f6e6187b3
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
using BioSequences: LongRNA

# ╔═╡ 55735555-c2f7-41f0-becd-cbad5717d7be
using DataFrames: DataFrame

# ╔═╡ eff20728-e3ee-470a-afc0-3e26a30c11be
using Distributions: Gamma

# ╔═╡ 746ffc2b-1be1-4f10-8d79-12364a8a8251
using Distributions: logpdf

# ╔═╡ b0997272-f623-43d8-a5ff-c1b1d6be99a4
using Distributions: pdf

# ╔═╡ a1807f68-6f16-4114-98ab-b075b45b7201
using Distributions: Poisson

# ╔═╡ 483d3154-5222-40b8-a380-aea037e4618c
using LinearAlgebra: Diagonal

# ╔═╡ 2132164b-51c2-4415-9d07-022a24011f3d
using LinearAlgebra: eigen

# ╔═╡ fb8e8cbd-050d-4d0e-81ac-32528f628a0e
using Makie: @L_str

# ╔═╡ 8973e48f-df81-4749-b71d-aea5ac4614b3
using NaNStatistics: nanmean, nansum

# ╔═╡ 96e03cbd-dc6f-4923-bb56-106ca66c1867
using Random: bitrand

# ╔═╡ e73b3886-162a-4a77-a28b-cdf269109b98
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2abcc581-8722-4cf7-bc09-8bf98a9b8648
using Statistics: cor

# ╔═╡ 6d6e3dc8-28eb-496d-9622-8128422ed2ed
using Statistics: mean

# ╔═╡ 7edc30fd-a2d6-4818-855c-c7f69c9f589b
using Statistics: middle

# ╔═╡ 4974c2e2-058d-41ca-924d-16709e4a58e6
using StatsBase: countmap

# ╔═╡ 94d99837-7415-4acb-b5d3-3b1dec5af05e
using Unitful: ustrip

# ╔═╡ 45907f4d-29ec-4be7-97e8-bfcb4695416b
import CairoMakie

# ╔═╡ f743167e-8552-4002-b9ca-c995cf3e9829
import CSV

# ╔═╡ 6be9532a-762f-4830-8d7f-175545fe75c9
import FASTX

# ╔═╡ 02b62a0f-159c-470d-91d0-673898ca277b
import HDF5

# ╔═╡ 651a9155-c7fc-4858-ae65-b433d1d116cc
import Infernal

# ╔═╡ c6371fcc-bb84-4ae4-a2d8-972b362c5350
import KernelDensity

# ╔═╡ 297c712c-e7e6-4d2f-a9c2-59b29665e6e3
import Makie

# ╔═╡ c53d715e-40d3-44cb-b85a-9c4c61d99819
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ f513852e-ce6a-41f0-98aa-155e06244aaf
import Rfam

# ╔═╡ ca499e53-296f-4cf3-9df1-5070e22fd6f0
import SamApp2024

# ╔═╡ 15974106-a5c6-4867-acfb-cdc69f5a57be
import StatsBase

# ╔═╡ b471828b-0822-4a08-a63c-5fa01e8d90b2
import ViennaRNA

# ╔═╡ 13f58dca-f8c6-40d9-bf2c-d20995181775
import ViennaRNA_jll

# ╔═╡ 2a3fff7d-a1e5-474e-8981-ac9c05494d4d
@show Rfam.get_rfam_directory();

# ╔═╡ 1679e271-e03c-4db7-ad7d-4aaf3a2c183a
@show Rfam.get_rfam_version();

# ╔═╡ c182ce25-a47f-4368-9ca7-9b5b93983dcb
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);



# ╔═╡ 3c96e229-e909-4de6-849e-753109b229d5
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 775e36c0-7c9b-47e7-adef-ac686ce0dc88
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 7b245a0d-f31b-41ae-9565-dd04f2a7ecfd
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 9d530b52-c9bd-4311-943f-61ec49786681
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 24905a44-cf3c-44a4-ae50-c067b1c2b763
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 59debd70-9eef-4a9d-b2f7-31337dac0664
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 47095bc6-336d-4097-bf61-8547fb9a4584
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ b3797c90-ff24-4c98-a05a-147dca6e66e7
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ 749e241a-769c-4f81-8e07-de89cfa222d2
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 012b1a9b-722e-4301-9043-ff1580c75236
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ ba2ea96a-aa4c-467d-82b5-688b46cc67bd
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 33b24793-a890-4b04-abe5-cbc2dc126633
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 6dcdfbb7-afaa-43e9-9f7d-a6691901e5a7
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ f74e788f-47a3-46d6-80da-8e3e8b50a029
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 802a5655-9b21-448b-933d-0471350f1058
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 9504ca2a-adad-4a04-a943-9b4cee91f834
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ e5230e7c-4c61-457f-b57c-07a5272074b5
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0];

# ╔═╡ fcd6da1f-2afb-4b47-8a39-4aacb2b48fb2
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0];

# ╔═╡ cb886986-01e2-494b-8d5d-0ceacaa6d886
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ 92975e4a-13a1-42b2-8219-a1266b7189a1
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ ac0f6b9c-abdc-4aa3-bf00-29f63330d219
_thresh = log(5)

# ╔═╡ 3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
_sites = SamApp.hallmark_sites_20230507;

# ╔═╡ 3bdc35bb-0b7d-48fe-87e4-3732fdb1feb8
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

sum(_responds_sam_yes_rep0)

# ╔═╡ 41caae29-9181-4cf8-a814-64868a9b8e90
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ cc12dd0a-c217-436f-870c-601f16ce735b
wuss = SamApp.rfam_ss("RF00162"; inserts=false)
ss = SamApp.clean_wuss(wuss)

# ╔═╡ 935b6930-8748-47c8-8db6-ee5557ef3423
p1_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p1;
p2_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p2;
p3_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p3;
p4_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p4;
pk_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().pk;

# ╔═╡ 387b740c-8475-4cdf-a08c-d8e5136bc7f2
ss_without_P1 = join([i ∈ p1_pos ? '.' : c for (i,c) in enumerate(ss)]);
ss_without_P2 = join([i ∈ p2_pos ? '.' : c for (i,c) in enumerate(ss)]);
ss_without_P3 = join([i ∈ p3_pos ? '.' : c for (i,c) in enumerate(ss)]);
ss_without_P4 = join([i ∈ p4_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 6212efdb-8f9a-4fc8-824c-0623ac71d2ce
ss_pk_only = replace(wuss, r"\(|\)|\[|\]|\{|\}|\<|\>|\-|\_|\," => '.', 'A' => '(', 'a' => ')')

# ╔═╡ a03f024e-d55d-4ffc-aedd-d96f5b539fdd
sampled_v = SamApp2024.rbm2022samples();

# ╔═╡ Cell order:
# ╠═77855f93-2b64-45bc-a307-df7f6e6187b3
# ╠═45907f4d-29ec-4be7-97e8-bfcb4695416b
# ╠═f743167e-8552-4002-b9ca-c995cf3e9829
# ╠═6be9532a-762f-4830-8d7f-175545fe75c9
# ╠═02b62a0f-159c-470d-91d0-673898ca277b
# ╠═651a9155-c7fc-4858-ae65-b433d1d116cc
# ╠═c6371fcc-bb84-4ae4-a2d8-972b362c5350
# ╠═297c712c-e7e6-4d2f-a9c2-59b29665e6e3
# ╠═c53d715e-40d3-44cb-b85a-9c4c61d99819
# ╠═f513852e-ce6a-41f0-98aa-155e06244aaf
# ╠═ca499e53-296f-4cf3-9df1-5070e22fd6f0
# ╠═15974106-a5c6-4867-acfb-cdc69f5a57be
# ╠═b471828b-0822-4a08-a63c-5fa01e8d90b2
# ╠═13f58dca-f8c6-40d9-bf2c-d20995181775
# ╠═c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
# ╠═55735555-c2f7-41f0-becd-cbad5717d7be
# ╠═eff20728-e3ee-470a-afc0-3e26a30c11be
# ╠═746ffc2b-1be1-4f10-8d79-12364a8a8251
# ╠═b0997272-f623-43d8-a5ff-c1b1d6be99a4
# ╠═a1807f68-6f16-4114-98ab-b075b45b7201
# ╠═483d3154-5222-40b8-a380-aea037e4618c
# ╠═2132164b-51c2-4415-9d07-022a24011f3d
# ╠═fb8e8cbd-050d-4d0e-81ac-32528f628a0e
# ╠═8973e48f-df81-4749-b71d-aea5ac4614b3
# ╠═96e03cbd-dc6f-4923-bb56-106ca66c1867
# ╠═e73b3886-162a-4a77-a28b-cdf269109b98
# ╠═2abcc581-8722-4cf7-bc09-8bf98a9b8648
# ╠═6d6e3dc8-28eb-496d-9622-8128422ed2ed
# ╠═7edc30fd-a2d6-4818-855c-c7f69c9f589b
# ╠═4974c2e2-058d-41ca-924d-16709e4a58e6
# ╠═94d99837-7415-4acb-b5d3-3b1dec5af05e
# ╠═2a3fff7d-a1e5-474e-8981-ac9c05494d4d
# ╠═1679e271-e03c-4db7-ad7d-4aaf3a2c183a
# ╠═c182ce25-a47f-4368-9ca7-9b5b93983dcb
# ╠═3c96e229-e909-4de6-849e-753109b229d5
# ╠═775e36c0-7c9b-47e7-adef-ac686ce0dc88
# ╠═7b245a0d-f31b-41ae-9565-dd04f2a7ecfd
# ╠═9d530b52-c9bd-4311-943f-61ec49786681
# ╠═24905a44-cf3c-44a4-ae50-c067b1c2b763
# ╠═59debd70-9eef-4a9d-b2f7-31337dac0664
# ╠═47095bc6-336d-4097-bf61-8547fb9a4584
# ╠═b3797c90-ff24-4c98-a05a-147dca6e66e7
# ╠═749e241a-769c-4f81-8e07-de89cfa222d2
# ╠═012b1a9b-722e-4301-9043-ff1580c75236
# ╠═ba2ea96a-aa4c-467d-82b5-688b46cc67bd
# ╠═33b24793-a890-4b04-abe5-cbc2dc126633
# ╠═6dcdfbb7-afaa-43e9-9f7d-a6691901e5a7
# ╠═f74e788f-47a3-46d6-80da-8e3e8b50a029
# ╠═802a5655-9b21-448b-933d-0471350f1058
# ╠═9504ca2a-adad-4a04-a943-9b4cee91f834
# ╠═e5230e7c-4c61-457f-b57c-07a5272074b5
# ╠═fcd6da1f-2afb-4b47-8a39-4aacb2b48fb2
# ╠═cb886986-01e2-494b-8d5d-0ceacaa6d886
# ╠═92975e4a-13a1-42b2-8219-a1266b7189a1
# ╠═ac0f6b9c-abdc-4aa3-bf00-29f63330d219
# ╠═3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
# ╠═3bdc35bb-0b7d-48fe-87e4-3732fdb1feb8
# ╠═41caae29-9181-4cf8-a814-64868a9b8e90
# ╠═cc12dd0a-c217-436f-870c-601f16ce735b
# ╠═935b6930-8748-47c8-8db6-ee5557ef3423
# ╠═387b740c-8475-4cdf-a08c-d8e5136bc7f2
# ╠═6212efdb-8f9a-4fc8-824c-0623ac71d2ce
# ╠═a03f024e-d55d-4ffc-aedd-d96f5b539fdd

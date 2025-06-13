### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ 77855f93-2b64-45bc-a307-df7f6e6187b3
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
using BioSequences: LongRNA

# ╔═╡ 55735555-c2f7-41f0-becd-cbad5717d7be
using DataFrames: DataFrame

# ╔═╡ fb8e8cbd-050d-4d0e-81ac-32528f628a0e
using Makie: @L_str

# ╔═╡ 8973e48f-df81-4749-b71d-aea5ac4614b3
using NaNStatistics: nanmean

# ╔═╡ 48ab42a5-3ea8-4b69-acfb-3e2471de7144
using NaNStatistics: nanstd

# ╔═╡ 61fcf845-c1b9-4798-8b11-1bee95c8d0e0
using NaNStatistics: nansum

# ╔═╡ e73b3886-162a-4a77-a28b-cdf269109b98
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2abcc581-8722-4cf7-bc09-8bf98a9b8648
using Statistics: cor

# ╔═╡ 6d6e3dc8-28eb-496d-9622-8128422ed2ed
using Statistics: mean

# ╔═╡ c7378daf-3efe-4317-a95e-cd538a64e1a6
using Statistics: std

# ╔═╡ 7edc30fd-a2d6-4818-855c-c7f69c9f589b
using Statistics: middle

# ╔═╡ 4974c2e2-058d-41ca-924d-16709e4a58e6
using StatsBase: countmap

# ╔═╡ 94d99837-7415-4acb-b5d3-3b1dec5af05e
using Unitful: ustrip

# ╔═╡ f6d726bd-4493-4aee-a824-a36c72b16e95
using SamApp2024: successes_tuple_str

# ╔═╡ a06dbeed-cd34-464f-95fc-f3659f95f760
md"# Imports"

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

# ╔═╡ 66de62df-32fb-4dc0-b3b5-05d0e5ca718c
import PlutoUI

# ╔═╡ 02bcf2aa-204e-4385-8790-3c76432deacf
PlutoUI.TableOfContents()

# ╔═╡ b02ae72e-1892-4fdc-8c71-4e2007c65895
md"# Load data"

# ╔═╡ 2a3fff7d-a1e5-474e-8981-ac9c05494d4d
@show Rfam.get_rfam_directory() Rfam.get_rfam_version();

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

# ╔═╡ 1978891c-c502-4a4d-bb0e-0b33635a3263
length(rbm_seqs)

# ╔═╡ ac0f6b9c-abdc-4aa3-bf00-29f63330d219
_thresh = log(5)

# ╔═╡ 3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ acf8bb76-bf87-4cf8-b280-e0e13c65fdff
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3))

# ╔═╡ 1a364c8f-49e3-4c33-b892-8f35106c4afa
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3))

# ╔═╡ cb29b45f-d7b1-4507-9531-5d439e4accbb
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 0e703742-3cb5-4701-97b9-6134d2ca30b6
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ e0223248-586e-48a1-9799-d3918f672e92
_inconclusive_rep0 = ((!).(_responds_sam_yes_rep0)) .& ((!).(_responds_sam_nop_rep0));

# ╔═╡ b8006eb5-b99f-4535-81ef-dfe669d0d332
_conclusive_rep0 = _responds_sam_yes_rep0 .| _responds_sam_nop_rep0;

# ╔═╡ db615f71-8858-47da-be2f-8196f7070e6d
sum(_responds_sam_yes_rep0), sum(_responds_sam_nop_rep0), sum(_inconclusive_rep0)

# ╔═╡ d896c790-02ab-4f57-801b-cf2eb25aefab
_responds_sam_yes_rep0[[299, 207]]

# ╔═╡ 41caae29-9181-4cf8-a814-64868a9b8e90
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ 3d0cd846-07a7-4910-a2a2-2ef660355ad5
wuss = SamApp2024.rfam_ss("RF00162"; inserts=false)

# ╔═╡ 74fbb354-2657-475e-8d1e-42d39e255396
ss = SamApp2024.clean_wuss(wuss)

# ╔═╡ 9f0baca9-46e2-4b20-9c0c-0dbbac36bbf0
p1_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p1;

# ╔═╡ ffcd720b-134f-4b83-adf7-ab39ab6dfba1
p2_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p2;

# ╔═╡ de174fc5-5bc0-4b5c-980d-88168f5355d1
p3_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p3;

# ╔═╡ bea93ba6-a1e5-444a-b419-f0ec290cb61f
p4_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().p4;

# ╔═╡ 2639031c-df64-47c4-9d75-4f2ad73291dd
pk_pos = SamApp2024.RF00162_sites_annotated_secondary_structure().pk;

# ╔═╡ 15ff286a-a8e4-4c8b-9c73-04450ad6c500
ss_without_P1 = join([i ∈ p1_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ b798af09-363d-4800-b37e-a8fe8ce36398
ss_without_P2 = join([i ∈ p2_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 1e382321-089d-407f-9d2b-e07dcdb00135
ss_without_P3 = join([i ∈ p3_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 445ab169-1024-4238-b3bb-35b9cbdc99e7
ss_without_P4 = join([i ∈ p4_pos ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ 6212efdb-8f9a-4fc8-824c-0623ac71d2ce
ss_pk_only = replace(wuss, r"\(|\)|\[|\]|\{|\}|\<|\>|\-|\_|\," => '.', 'A' => '(', 'a' => ')')

# ╔═╡ a03f024e-d55d-4ffc-aedd-d96f5b539fdd
sampled_v = SamApp2024.rbm2022samples();

# ╔═╡ 0f849d0c-0a3c-4a70-b437-98b6340ceb58
Vienna_energies_fold = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 7cc5d85d-27c7-4f4d-b62a-30925a58a98a
Vienna_energies_P1 = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss)) - ustrip(ViennaRNA.energy(string(seq), ss_without_P1)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 0d0bde49-85c7-4559-96c5-2602e59b505d
Vienna_energies_P2 = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss)) - ustrip(ViennaRNA.energy(string(seq), ss_without_P2)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 9bcf2ca1-b26b-453c-824f-43de4e31c967
Vienna_energies_P3 = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss)) - ustrip(ViennaRNA.energy(string(seq), ss_without_P3)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 94310e4e-644a-4b99-ad5f-d708a47ad7f8
Vienna_energies_P4 = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss)) - ustrip(ViennaRNA.energy(string(seq), ss_without_P4)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 7dcca23e-dfb0-407c-9e91-47f7e8390990
Vienna_energies_Pk = [ismissing(seq) ? missing : ustrip(ViennaRNA.energy(string(seq), ss_pk_only)) for seq = shape_data_rep0.aligned_sequences];

# ╔═╡ 461c7485-2551-4417-9ffc-7cad3f65d2bb
@time Vienna_energies_P1_RBM_samples = [
    ustrip(ViennaRNA.energy(string(seq), ss)) - ustrip(ViennaRNA.energy(string(seq), ss_without_P1))
    for seq = SamApp2024.rnaseq(sampled_v)
];

# ╔═╡ f547cf7b-943b-400e-aa4c-f1d60237c984
Vienna_energies_Pk_RBM_samples = [ustrip(ViennaRNA.energy(string(seq), ss_pk_only)) for seq = SamApp2024.rnaseq(sampled_v)];

# ╔═╡ 1f5d3f1b-5ac7-4cef-a746-eaadc19c94d0
Vienna_energies_Pk_RNAeval = [ismissing(seq) ? NaN : SamApp2024.vienna_pk_binding_energy_rnaeval(seq) for seq = shape_data_rep0.aligned_sequences]

# ╔═╡ 56444e37-03d6-4829-9776-f79dc1a839a2
Vienna_energies_Pk_RBM_samples_RNAeval = [SamApp2024.vienna_pk_binding_energy_rnaeval(string(seq)) for seq = SamApp2024.rnaseq(sampled_v)];

# ╔═╡ 83534334-93f5-422b-8280-a2ca8ccaba4f
# All merged data, for the reactivity profiles plots
shape_data_all_merged = SamApp2024.load_shapemapper_data_pierre_demux_20231027_repls_merged();

# ╔═╡ e346e39a-185b-4fae-9e16-3745bd82a906
conds_SAM_all_merged = map(identity, indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_allrep", "SAMAP_1M7_1SAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ ebede37d-ee91-4942-86bf-3924b0f8578f
conds_Mg_all_merged = map(identity, indexin(["SAMAP_1M7_noSAM_5Mg_T30C_allrep"], shape_data_all_merged.conditions));

# ╔═╡ ea6bb551-f3af-4647-9620-a3ecf90380b0
conds_SAM_all_merged, conds_Mg_all_merged

# ╔═╡ 8c4a1b15-13e8-43a9-9be2-cb614a8c1869
bps_reactivities_merged = shape_data_all_merged.shape_reactivities[bps, nat_seqs, conds_SAM_all_merged];

# ╔═╡ 44321d83-09e4-4928-ac36-b32a12297320
nps_reactivities_merged = shape_data_all_merged.shape_reactivities[nps, nat_seqs, conds_SAM_all_merged];

# ╔═╡ d5b6b49a-9721-4c4c-a735-743f3cd852ab
all_reactivities_merged = shape_data_all_merged.shape_reactivities[:, nat_seqs, conds_SAM_all_merged];

# ╔═╡ 309bf006-c0e1-46a6-9bc8-0834450340d9
shape_stats_merged = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_all_merged,
    paired_reactivities = bps_reactivities_merged,
    unpaired_reactivities = nps_reactivities_merged,
    all_reactivities = all_reactivities_merged,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 32d66706-b17d-4c58-8464-1128348ce06a
findall(_responds_sam_yes_rep0) ∩ rbm_seqs

# ╔═╡ c87fa471-dbb6-4802-8bd0-aae1eafd7fc0
md"# Plot"

# ╔═╡ 4c4144ca-9f36-44f0-b697-7411df2fac6a
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

# ╔═╡ 9fcca185-180f-4478-8fa2-d0cf37e1937b
let fig = Makie.Figure(; halign = :left)
	ax = Makie.Axis(fig[1,1][1,1]; halign=:left, width=200, height=200, xlabel="Protect. score(hallmark sites)", ylabel="RBM score", title="no SAM", xgridvisible=false, ygridvisible=false)
	#Makie.scatter!(ax, x_mg_rep0, -aptamer_rbm_energies, markersize=10, color=(:silver, 0.3), label="All probed")
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_nop_rep0) ∩ nat_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ nat_seqs], markersize=15, color=:green, marker='O')
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_yes_rep0) ∩ nat_seqs], -aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ nat_seqs], markersize=15, color=:green, marker='●')
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_nop_rep0) ∩ inf_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ inf_seqs], markersize=15, color=:red, marker='O')
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_yes_rep0) ∩ inf_seqs], -Float64.(aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ inf_seqs]), markersize=15, color=:red, marker='●') # empty
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_nop_rep0) ∩ rbm_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ rbm_seqs], markersize=15, color=:blue, marker='O')
	Makie.scatter!(ax, x_mg_rep0[findall(_responds_sam_yes_rep0) ∩ rbm_seqs], -aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ rbm_seqs], markersize=15, color=:blue, marker='●')
	Makie.vlines!(ax, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	
	ax = Makie.Axis(fig[1,1][1,2]; halign=:left, width=200, height=200, xlabel="Protect. score(hallmark sites)", ylabel="RBM score", title="with SAM", xgridvisible=false, ygridvisible=false)
	#Makie.scatter!(ax, x_sam_rep0, -aptamer_rbm_energies, markersize=10, color=(:gray, 0.5), label="All probed")
	plt1 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_nop_rep0) ∩ nat_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ nat_seqs], markersize=15, color=:green, marker='O', label="Nat. (❌)")
	plt2 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_yes_rep0) ∩ nat_seqs], -aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ nat_seqs], markersize=15, color=:green, marker='●', label="Nat. (✓)")
	plt3 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_nop_rep0) ∩ inf_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ inf_seqs], markersize=15, color=:red, marker='O', label="CM (❌)")
	plt4 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_yes_rep0) ∩ inf_seqs], -Float64.(aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ inf_seqs]), markersize=15, color=:red, marker='●', label="CM (✓)") # Empty
	plt5 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_nop_rep0) ∩ rbm_seqs], -aptamer_rbm_energies[findall(_responds_sam_nop_rep0) ∩ rbm_seqs], markersize=15, color=:blue, marker='O', label="RBM (❌)")
	plt6 = Makie.scatter!(ax, x_sam_rep0[findall(_responds_sam_yes_rep0) ∩ rbm_seqs], -aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ rbm_seqs], markersize=15, color=:blue, marker='●', label="RBM (✓)")
	Makie.vlines!(ax, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	#Makie.xlims!(ax, -78, 35)
	Makie.hideydecorations!(ax)
	#Makie.axislegend(ax, position=:rt, framevisible=false)
	Makie.Legend(fig[1,1][1,3],
	    [plt1, plt2, plt3, plt4, plt5, plt6],
	    ["Nat. (❌)", "Nat. (✓)", "rCM (❌)", "rCM (✓)", "RBM (❌)", "RBM (✓)"],
	    framevisible=false
	)
	
	# _dummy_ax = Makie.Axis(fig[1,2], width=400, height=300, xgridvisible=false, ygridvisible=false) # placeholder for the table
	# Makie.hidexdecorations!(_dummy_ax)
	# Makie.hideydecorations!(_dummy_ax)
	# Makie.hidespines!(_dummy_ax, :t, :b, :l, :r)
	
	n_ex_rbm = 299 # rbm example, switcher
	#n_ex_nat = 112 # natural example, responsive but not switcher
	#n_ex_nat = 101 # natural example, responsive but not switcher
	n_ex_rbm_2 = 207 # other RBM example, switcher
	
	_width = 550
	_height = 70
	
	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex_rbm, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex_rbm, only(conds_Mg_all_merged)]
	
	ax_react_1 = Makie.Axis(fig[2,1][1,1], width=_width, height=_height, xticks=5:10:108, ylabel="reactivity", xgridvisible=false, ygridvisible=false, yticks=0:2:5, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_1, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_1, 1:108, _R_mg, color=:gray, step=:center)
	Makie.stairs!(ax_react_1, 1:108, _R_sam, color=:purple, step=:center)
	Makie.xlims!(ax_react_1, 0.5, 108.5)
	Makie.ylims!(ax_react_1, -1, 4)
	Makie.hidespines!(ax_react_1, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_1)
	
	ax_diff_1 = Makie.Axis(fig[2,1][2,1]; width=_width, height=_height, xticks=5:10:108, xlabel="site", ylabel="Δreactivity", xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true, yticks=-2:2)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_1, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_1, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.xlims!(ax_diff_1, 0.5, 108.5)
	Makie.ylims!(ax_diff_1, -1.5, 1.5)
	Makie.hidespines!(ax_diff_1, :t, :b, :r)
	Makie.hidexdecorations!(ax_diff_1)
	
	_R_sam = shape_data_all_merged.shape_reactivities[:, n_ex_rbm_2, conds_SAM_all_merged[1]]
	_R_mg = shape_data_all_merged.shape_reactivities[:, n_ex_rbm_2, only(conds_Mg_all_merged)]
	
	ax_react_2 = Makie.Axis(fig[2,1][3,1], width=_width, height=_height, xticks=10:10:108, ylabel="reactivity", xgridvisible=false, ygridvisible=false, yticks=0:2:5, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_react_2, x0, xf; color=(color, alpha))
	end
	Makie.stairs!(ax_react_2, 1:108, _R_mg, color=:gray, step=:center)
	Makie.stairs!(ax_react_2, 1:108, _R_sam, color=:purple, step=:center)
	Makie.xlims!(ax_react_2, 0.5, 108.5)
	Makie.ylims!(ax_react_2, -1, 6)
	Makie.hidespines!(ax_react_2, :t, :r, :b)
	Makie.hidexdecorations!(ax_react_2)
	#Makie.hideydecorations!(ax_react_2)
	
	ax_diff_2 = Makie.Axis(fig[2,1][4,1]; width=_width, height=_height, xticks=5:10:108, xlabel="site", ylabel="Δreactivity", xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true, yticks=-2:2)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_diff_2, x0, xf; color=(color, alpha))
	end
	Makie.barplot!(ax_diff_2, 1:108, _R_sam - _R_mg, color=ifelse.(_R_sam - _R_mg .< 0, :green, :gray))
	Makie.scatter!(ax_diff_2, _sites, -1.4one.(_sites), markersize=7, color=:black, marker=:utriangle)
	Makie.xlims!(ax_diff_2, 0, 109)
	Makie.ylims!(ax_diff_2, -1.5, 1)
	Makie.hidespines!(ax_diff_2, :t, :r)
	#Makie.hideydecorations!(ax_diff_2)
	
	# Makie.scatter!(ax_diff_1, _sites, one.(_sites), color=:blue, markersize=5)
	# Makie.scatter!(ax_diff_2, _sites, one.(_sites), color=:blue, markersize=5)
	
	Makie.linkxaxes!(ax_react_1, ax_diff_1)
	Makie.linkxaxes!(ax_react_2, ax_diff_2)
	Makie.linkyaxes!(ax_react_1, ax_react_2)
	Makie.linkyaxes!(ax_diff_1, ax_diff_2)
	
	Makie.resize_to_layout!(fig)
	Makie.save("Figures/SAM response Repl0 v2.pdf", fig)
	fig
end

# ╔═╡ d3b705a7-03ac-4731-8c45-fa95c857d004
Float64.(aptamer_rbm_energies[findall(_responds_sam_yes_rep0) ∩ inf_seqs])

# ╔═╡ cb8c0927-cdda-4a9f-aad0-d6b8308ae933
md"# Pairing energies"

# ╔═╡ 935a3442-ed13-4239-962b-9c17dc86e792
let fig = Makie.Figure(; halign = :left)

	_sz = 200
	
	x_mg_p1 = nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_mg_rep0]; dim=(1,3))
	x_sam_p1 = nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_sam_rep0]; dim=(1,3))
	
	# responsive to SAM
	_resp_sam_P1 = (x_mg_p1 .< -_thresh) .& (x_sam_p1 .> _thresh)
	_stuck_open_P1 = (x_mg_p1 .< -_thresh) .& (x_sam_p1 .< -_thresh)
	_stuck_closed_P1 = (x_mg_p1 .> _thresh) .& (x_sam_p1 .> _thresh)
	
	_e_p1_bands = [-10, 0]
	
	ax_rbm_hist = Makie.Axis(fig[1,1]; halign=:left, width=75, height=_sz, yticks=-15:5:5, xticks=[-0.2, 0, 0.2], xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax_rbm_hist, Vienna_energies_P1_RBM_samples, color=:lightblue, normalization=:pdf, bins=-15:1:7, direction=:x, scale_to=-1)
	Makie.hlines!(ax_rbm_hist, _e_p1_bands, color=:blue, linestyle=:dash)
	Makie.hidespines!(ax_rbm_hist, :t, :r, :l, :b)
	Makie.hidexdecorations!(ax_rbm_hist)
	Makie.hideydecorations!(ax_rbm_hist)
	
	ax_mg = Makie.Axis(fig[1,2]; halign=:left, width=_sz, height=_sz, xlabel="Protect. score(P1)", ylabel="P1 pairing energy (kcal/mol)", xgridvisible=false, ygridvisible=false)
	#Makie.scatter!(ax, nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_mg_rep0]; dim=(1,3)), Vienna_energies_P1, color=(:gray, 0.5), markersize=15)
	plt1 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _stuck_closed_P1, conds_mg_rep0]; dim=(1,3)), Vienna_energies_P1[_stuck_closed_P1], color=:black, markersize=15, label="P1 closed", marker='O')
	plt2 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _stuck_open_P1, conds_mg_rep0]; dim=(1,3)), Vienna_energies_P1[_stuck_open_P1], color=:gray, markersize=15, label="P1 open", marker='O')
	plt3 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _resp_sam_P1, conds_mg_rep0]; dim=(1,3)), Vienna_energies_P1[_resp_sam_P1], color=:steelblue, markersize=15, label="P1 switch", marker='●')
	Makie.vlines!(ax_mg, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	Makie.hlines!(ax_mg, _e_p1_bands, color=:blue, linestyle=:dash)
	Makie.Legend(fig[1,4], [plt1, plt2, plt3], ["P1 closed", "P1 open", "P1 switch"], framevisible=false, valign=:top)
	
	ax_sam = Makie.Axis(fig[1,3]; halign=:left, width=_sz, height=_sz, xlabel="Protect. score(P1)", ylabel="P1 pairing energy (kcal/mol)", xgridvisible=false, ygridvisible=false)
	#Makie.scatter!(ax, nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_sam_rep0]; dim=(1,3)), Vienna_energies_P1, color=(:gray, 0.5),  markersize=15)
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _stuck_closed_P1, conds_sam_rep0]; dim=(1,3)), Vienna_energies_P1[_stuck_closed_P1], color=:black, markersize=15, marker='O')
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _stuck_open_P1, conds_sam_rep0]; dim=(1,3)), Vienna_energies_P1[_stuck_open_P1], color=:gray, markersize=15, marker='O')
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[p1_pos, _resp_sam_P1, conds_sam_rep0]; dim=(1,3)), Vienna_energies_P1[_resp_sam_P1], color=:steelblue, markersize=15, marker='●')
	# Makie.hlines!(ax, nanmean(Vienna_energies_P1[_stuck_open]), color=:red, linewidth=2)
	# Makie.hlines!(ax, nanmean(Vienna_energies_P1[_stuck_closed]), color=:green, linewidth=2)
	# Makie.hlines!(ax, nanmean(Vienna_energies_P1[_resp_sam]), color=:blue, linewidth=2)
	#nanmean(Vienna_energies_P1[_resp_sam]), nanmean(Vienna_energies_P1[_stuck_closed])
	Makie.vlines!(ax_sam, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	Makie.hlines!(ax_sam, _e_p1_bands, color=:blue, linestyle=:dash)
	Makie.hideydecorations!(ax_sam)
	
	Makie.linkyaxes!(ax_rbm_hist, ax_mg, ax_sam)
	
	
	x_mg_pk = nansum(shape_stats_rep0.shape_log_odds[pk_pos, :, conds_mg_rep0]; dim=(1,3))
	x_sam_pk = nansum(shape_stats_rep0.shape_log_odds[pk_pos, :, conds_sam_rep0]; dim=(1,3))
	
	# responsive to SAM
	_resp_sam_Pk = (x_mg_pk .< -_thresh) .& (x_sam_pk .> _thresh)
	_stuck_open_Pk = (x_mg_pk .< -_thresh) .& (x_sam_pk .< -_thresh)
	_stuck_closed_Pk = (x_mg_pk .> _thresh) .& (x_sam_pk .> _thresh)
	
	_e_pk_bands = [-8, -3]
	
	
	ax_rbm_hist = Makie.Axis(fig[2,1]; halign=:left, width=75, height=_sz, yticks=-15:5:5, xticks=[-0.2, 0, 0.2], xgridvisible=false, ygridvisible=false, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax_rbm_hist, Vienna_energies_P1_RBM_samples, color=:lightblue, normalization=:pdf, bins=-15:1:7, direction=:x, scale_to=-1)
	Makie.hlines!(ax_rbm_hist, _e_pk_bands, color=:blue, linestyle=:dash)
	Makie.hidespines!(ax_rbm_hist, :t, :b, :r, :l)
	Makie.hidexdecorations!(ax_rbm_hist)
	Makie.hideydecorations!(ax_rbm_hist)
	
	ax_mg = Makie.Axis(fig[2,2]; halign=:left, width=_sz, height=_sz, xlabel="Protect. score(Pk)", ylabel="Pk pairing energy (kcal/mol)", xgridvisible=false, ygridvisible=false)
	plt1 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _stuck_closed_Pk, conds_mg_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_stuck_closed_Pk], color=:black, markersize=15, label="Pk closed", marker='O')
	plt2 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _stuck_open_Pk, conds_mg_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_stuck_open_Pk], color=:gray, markersize=15, label="Pk open", marker='O')
	plt3 = Makie.scatter!(ax_mg, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _resp_sam_Pk, conds_mg_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_resp_sam_Pk], color=:steelblue, markersize=15, label="Pk switch", marker='●')
	Makie.vlines!(ax_mg, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	Makie.hlines!(ax_mg, _e_pk_bands, color=:blue, linestyle=:dash)
	Makie.Legend(fig[2,4], [plt1, plt2, plt3], ["Pk closed", "Pk open", "Pk switch"], framevisible=false, valign=:top)
	
	ax_sam = Makie.Axis(fig[2,3]; halign=:left, width=_sz, height=_sz, xlabel="Protect. score(Pk)", ylabel="Pk pairing energy (kcal/mol)", xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _stuck_closed_Pk, conds_sam_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_stuck_closed_Pk], color=:black, markersize=15, marker='O')
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _stuck_open_Pk, conds_sam_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_stuck_open_Pk], color=:gray, markersize=15, marker='O')
	Makie.scatter!(ax_sam, nansum(shape_stats_rep0.shape_log_odds[pk_pos, _resp_sam_Pk, conds_sam_rep0]; dim=(1,3)), Vienna_energies_Pk_RNAeval[_resp_sam_Pk], color=:steelblue, markersize=15, marker='●')
	Makie.vlines!(ax_sam, [-_thresh, _thresh], linestyle=:dash, color=:orange)
	Makie.hlines!(ax_sam, _e_pk_bands, color=:blue, linestyle=:dash)
	Makie.hideydecorations!(ax_sam)
	
	Makie.ylims!(ax_mg, -10, 5)
	Makie.ylims!(ax_sam, -10, 5)
	Makie.linkyaxes!(ax_rbm_hist, ax_mg, ax_sam)
	
	
	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/SAM response Repl0 -- Vienna.pdf", fig)
	fig
end

# ╔═╡ fbbd6109-ac4e-4676-97b1-e9fef5a3b729
length(Vienna_energies_P1)

# ╔═╡ 52bd54f9-0df6-4c26-849a-64928693fd9a
Vienna_energies_P1[rbm_seqs]

# ╔═╡ 8ff5340c-3ade-4c03-9a72-1c324aa0e882
Vienna_energies_P1[(!ismissing).(Vienna_energies_P1) .&& (Vienna_energies_P1 .> 0)]

# ╔═╡ 36cf3060-1bfd-4e2e-8009-a9d5bdaf19cb
shape_data_045.aptamer_origin[(!ismissing).(Vienna_energies_P1) .&& (Vienna_energies_P1 .> 0)]

# ╔═╡ 9255fa38-9004-43c9-8817-5657e6cfb2d5
md"# Read depths"

# ╔═╡ 6b38d182-a901-4b20-bfdf-9441e4586e7b
let fig = Makie.Figure()
	width = 550
	height = 100
	xticks = 5:10:108
	yticks = (0:10000:30000, ["0", "10000", "20000", "30000"])

	depth_M_avg = nanmean(shape_data_rep0.shape_M_depth; dim=(2,3))
	depth_U_avg = nanmean(shape_data_rep0.shape_U_depth; dim=(2,3))
	depth_D_avg = nanmean(shape_data_rep0.shape_D_depth; dim=(2,3))
	
	depth_M_err = nanstd(shape_data_rep0.shape_M_depth; dim=(2,3))
	depth_U_err = nanstd(shape_data_rep0.shape_U_depth; dim=(2,3))
	depth_D_err = nanstd(shape_data_rep0.shape_D_depth; dim=(2,3))
		
	ax_M = Makie.Axis(fig[1,1]; width, height, xticks, yticks, ylabel="Read depth (M)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_U = Makie.Axis(fig[2,1]; width, height, xticks, yticks, ylabel="Read depth (U)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_D = Makie.Axis(fig[3,1]; width, height, xticks, yticks, ylabel="Read depth (D)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_M, x0, xf; color=(color, alpha))
	    Makie.vspan!(ax_U, x0, xf; color=(color, alpha))
	    Makie.vspan!(ax_D, x0, xf; color=(color, alpha))
	end

	Makie.band!(ax_M, 1:108, depth_M_avg - depth_M_err, depth_M_avg + depth_M_err; color=(:gray, 0.3))
	Makie.band!(ax_U, 1:108, depth_U_avg - depth_U_err, depth_U_avg + depth_U_err; color=(:gray, 0.3))
	Makie.band!(ax_D, 1:108, depth_D_avg - depth_D_err, depth_D_avg + depth_D_err; color=(:gray, 0.3))
	
	Makie.stairs!(ax_M, 1:108, depth_M_avg, color=:black, step=:center)
	Makie.stairs!(ax_U, 1:108, depth_U_avg, color=:black, step=:center)
	Makie.stairs!(ax_D, 1:108, depth_D_avg, color=:black, step=:center)

	for ax = (ax_M, ax_U, ax_D)
		Makie.xlims!(ax, 0.5, 108.5)
	end
	#Makie.ylims!(ax_M, -1, 4)
	#Makie.hidespines!(ax_M, :t, :r, :b)
	#Makie.hidexdecorations!(ax_M)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 2cc26b9a-c18b-4321-ad79-a00f8085b525
let fig = Makie.Figure()
	width = 600
	height = 100
	xticks = 5:10:108
	yticks = (0:10000:30000, ["0", "10000", "20000", "30000"])

	depth_M_avg = dropdims(mean(replace(shape_data_rep0.shape_M_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	depth_U_avg = dropdims(mean(replace(shape_data_rep0.shape_U_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	depth_D_avg = dropdims(mean(replace(shape_data_rep0.shape_D_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	
	depth_M_err = dropdims(std(replace(shape_data_rep0.shape_M_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	depth_U_err = dropdims(std(replace(shape_data_rep0.shape_U_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	depth_D_err = dropdims(std(replace(shape_data_rep0.shape_D_depth, NaN => 0); dims=(2,3)); dims=(2,3))
	
	ax_M = Makie.Axis(fig[1,1]; width, height, xticks, yticks, ylabel="Read depth (M)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_U = Makie.Axis(fig[2,1]; width, height, xticks, yticks, ylabel="Read depth (U)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_D = Makie.Axis(fig[3,1]; width, height, xticks, yticks, ylabel="Read depth (D)", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	for (x0, xf, color, alpha) = struct_bands
	    Makie.vspan!(ax_M, x0, xf; color=(color, alpha))
	    Makie.vspan!(ax_U, x0, xf; color=(color, alpha))
	    Makie.vspan!(ax_D, x0, xf; color=(color, alpha))
	end

	Makie.band!(ax_M, 1:108, depth_M_avg - depth_M_err, depth_M_avg + depth_M_err; color=(:gray, 0.3))
	Makie.band!(ax_U, 1:108, depth_U_avg - depth_U_err, depth_U_avg + depth_U_err; color=(:gray, 0.3))
	Makie.band!(ax_D, 1:108, depth_D_avg - depth_D_err, depth_D_avg + depth_D_err; color=(:gray, 0.3))
	
	Makie.stairs!(ax_M, 1:108, depth_M_avg; color=:black, step=:center, linewidth=2)
	Makie.stairs!(ax_U, 1:108, depth_U_avg; color=:black, step=:center, linewidth=2)
	Makie.stairs!(ax_D, 1:108, depth_D_avg; color=:black, step=:center, linewidth=2)

	for ax = (ax_M, ax_U, ax_D)
		Makie.xlims!(ax, 0.5, 108.5)
	end
	#Makie.ylims!(ax_M, -1, 4)
	#Makie.hidespines!(ax_M, :t, :r, :b)
	#Makie.hidexdecorations!(ax_M)

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 14, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "B)", fontsize = 14, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[3,1][1, 1, Makie.TopLeft()], "C)", fontsize = 14, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA/cossio/SAM/2024/SamApp2024.jl/pluto/SI/Figures/read_depths.pdf", fig)
	fig
end

# ╔═╡ cd53cfd1-ec82-4c3b-93fe-62f0138c6f25
let fig = Makie.Figure()
	sz = 300
	cutoff = 20
	bins = 0:1e3:1e5
	
	ax_M = Makie.Axis(fig[1,1]; width=sz, height=sz, ylabel="Frequency", title="M", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_U = Makie.Axis(fig[1,2]; width=sz, height=sz, ylabel="Frequency", title="U", xgridvisible=false, ygridvisible=false, ytrimspine=true)
	ax_D = Makie.Axis(fig[1,3]; width=sz, height=sz, ylabel="Frequency", title="D", xgridvisible=false, ygridvisible=false, ytrimspine=true)

	Makie.stephist!(ax_M, vec(replace(shape_data_rep0.shape_M_depth[begin:cutoff, :, :], NaN => 0)); bins, normalization=:pdf, color=:red)
	Makie.stephist!(ax_M, vec(replace(shape_data_rep0.shape_M_depth[cutoff:end, :, :], NaN => 0)); bins, normalization=:pdf, color=:blue)

	Makie.stephist!(ax_U, filter(isfinite, shape_data_rep0.shape_U_depth[begin:cutoff, :, :]); bins, normalization=:pdf, color=:red)
	Makie.stephist!(ax_U, filter(isfinite, shape_data_rep0.shape_U_depth[cutoff:end, :, :]); bins, normalization=:pdf, color=:blue)

	Makie.stephist!(ax_D, filter(isfinite, shape_data_rep0.shape_D_depth[begin:cutoff, :, :]); bins, normalization=:pdf, color=:red)
	Makie.stephist!(ax_D, filter(isfinite, shape_data_rep0.shape_D_depth[cutoff:end, :, :]); bins, normalization=:pdf, color=:blue)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ d7cd6fac-95f6-4612-bd5a-bb6b839e3f2c
let fig = Makie.Figure()
	width = 550
	height = 100
	xticks = 0:5:108
	
	ax_M = Makie.Axis(fig[1,1]; width, height, xticks, ylabel="Null fraction", title="M", xgridvisible=false, ygridvisible=false, ytrimspine=true, yscale=log10)
	ax_U = Makie.Axis(fig[2,1]; width, height, xticks, ylabel="Null fraction", title="U", xgridvisible=false, ygridvisible=false, ytrimspine=true, yscale=log10)
	ax_D = Makie.Axis(fig[3,1]; width, height, xticks, ylabel="Null fraction", title="D", xgridvisible=false, ygridvisible=false, ytrimspine=true, yscale=log10)

	Makie.lines!(ax_M, 1:108, dropdims(1 .- mean((shape_data_rep0.shape_M_depth .> 0) .&& (shape_data_rep0.shape_M .> 0); dims=(2,3)); dims=(2,3)))
	Makie.lines!(ax_U, 1:108, dropdims(1 .- mean((shape_data_rep0.shape_U_depth .> 0) .&& (shape_data_rep0.shape_U .> 0); dims=(2,3)); dims=(2,3)))
	Makie.lines!(ax_D, 1:108, dropdims(1 .- mean((shape_data_rep0.shape_D_depth .> 0) .&& (shape_data_rep0.shape_D .> 0); dims=(2,3)); dims=(2,3)))

	for ax = (ax_M, ax_U, ax_D)
		Makie.xlims!(ax, 0.5, 108.5)
		Makie.ylims!(ax, 0.001, 1)
	end

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 34cc67fd-f567-425b-bad0-117b00c85a1d
md"# Counts"

# ╔═╡ dd25496f-e430-42df-a0b4-143bcb26f9f3
begin
	println("Replicate 0")

	println("Nat. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[nat_seqs]),
	        sum(_responds_sam_nop_rep0[nat_seqs])
	    )
	)
	println("Nat seed. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[seed_seqs]),
	        sum(_responds_sam_nop_rep0[seed_seqs])
	    )
	)
	println("Nat full. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[full_seqs]),
	        sum(_responds_sam_nop_rep0[full_seqs])
	    )
	)
	println("Nat. (E<-300): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ nat_seqs]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ nat_seqs])
	    )
	)
	println("Nat. (E<-310): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ nat_seqs]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ nat_seqs])
	    )
	)
	
	println("RBM (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[rbm_seqs]),
	        sum(_responds_sam_nop_rep0[rbm_seqs])
	    )
	)
	println("RBM (E<-300): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ rbm_seqs]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false)) ∩ rbm_seqs])
	    )
	)
	println("RBM (E<-310): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ rbm_seqs]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false)) ∩ rbm_seqs])
	    )
	)
	
	println("Infernal: ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[inf_seqs]),
	        sum(_responds_sam_nop_rep0[inf_seqs])
	    )
	)
	
	println("All. (all): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0),
	        sum(_responds_sam_nop_rep0)
	    )
	)
	println("All. (E<-300): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false))]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -300, missing => false))])
	    )
	)
	println("All. (E<-310): ", successes_tuple_str(
	        sum(_responds_sam_yes_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false))]),
	        sum(_responds_sam_nop_rep0[findall(replace(aptamer_rbm_energies .< -310, missing => false))])
	    )
	)
end

# ╔═╡ fddb0aa7-e722-4f8e-88f3-9255275f8da0
md"# Pairing energies, for inconclusives"

# ╔═╡ e84e8b68-bc2a-412b-b487-3b5163131ca5
let fig = Makie.Figure(; halign = :left)

	_sz = 300
	
	x_mg_p1 = nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_mg_rep0]; dim=(1,3))
	x_sam_p1 = nansum(shape_stats_rep0.shape_log_odds[p1_pos, :, conds_sam_rep0]; dim=(1,3))
	
	# responsive to SAM
	_resp_sam_P1_yes = (x_mg_p1 .< -_thresh) .& (x_sam_p1 .> _thresh)
	_resp_sam_P1_nop = (x_mg_p1 .> +_thresh) .| (x_sam_p1 .< -_thresh)
	_resp_sam_P1_incl = ((!).(_resp_sam_P1_yes)) .& ((!).(_resp_sam_P1_nop))

	ax = Makie.Axis(fig[1,1]; width=_sz, height=_sz, xlabel="P1 pairing energy", ylabel="RBM score")
	Makie.scatter!(Vienna_energies_P1, -aptamer_rbm_energies; color=:gray, label="all")
	Makie.scatter!(Vienna_energies_P1[_resp_sam_P1_incl], -aptamer_rbm_energies[_resp_sam_P1_incl]; color=:red, label="inconcl.")
	Makie.axislegend(ax, position=:rb)


	x_mg_pk = nansum(shape_stats_rep0.shape_log_odds[pk_pos, :, conds_mg_rep0]; dim=(1,3))
	x_sam_pk = nansum(shape_stats_rep0.shape_log_odds[pk_pos, :, conds_sam_rep0]; dim=(1,3))
	
	# responsive to SAM
	_resp_sam_Pk_yes = (x_mg_pk .< -_thresh) .& (x_sam_pk .> _thresh)
	_resp_sam_Pk_nop = (x_mg_pk .> +_thresh) .| (x_sam_pk .< -_thresh)
	_resp_sam_Pk_incl = ((!).(_resp_sam_Pk_yes)) .& ((!).(_resp_sam_Pk_nop))

	ax = Makie.Axis(fig[1,2]; width=_sz, height=_sz, xlabel="Pk pairing energy", ylabel="RBM score")
	Makie.scatter!(Vienna_energies_Pk_RNAeval, -aptamer_rbm_energies; color=:gray, label="all")
	Makie.scatter!(Vienna_energies_Pk_RNAeval[_resp_sam_Pk_incl], -aptamer_rbm_energies[_resp_sam_Pk_incl]; color=:red, label="inconcl.")
	Makie.xlims!(ax, -10, 5)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/SAM response Repl0 -- Vienna.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═a06dbeed-cd34-464f-95fc-f3659f95f760
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
# ╠═66de62df-32fb-4dc0-b3b5-05d0e5ca718c
# ╠═c50ef155-e21a-4fa6-9d59-4ff6aa870f1e
# ╠═55735555-c2f7-41f0-becd-cbad5717d7be
# ╠═fb8e8cbd-050d-4d0e-81ac-32528f628a0e
# ╠═8973e48f-df81-4749-b71d-aea5ac4614b3
# ╠═48ab42a5-3ea8-4b69-acfb-3e2471de7144
# ╠═61fcf845-c1b9-4798-8b11-1bee95c8d0e0
# ╠═e73b3886-162a-4a77-a28b-cdf269109b98
# ╠═2abcc581-8722-4cf7-bc09-8bf98a9b8648
# ╠═6d6e3dc8-28eb-496d-9622-8128422ed2ed
# ╠═c7378daf-3efe-4317-a95e-cd538a64e1a6
# ╠═7edc30fd-a2d6-4818-855c-c7f69c9f589b
# ╠═4974c2e2-058d-41ca-924d-16709e4a58e6
# ╠═94d99837-7415-4acb-b5d3-3b1dec5af05e
# ╠═f6d726bd-4493-4aee-a824-a36c72b16e95
# ╠═02bcf2aa-204e-4385-8790-3c76432deacf
# ╠═b02ae72e-1892-4fdc-8c71-4e2007c65895
# ╠═2a3fff7d-a1e5-474e-8981-ac9c05494d4d
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
# ╠═1978891c-c502-4a4d-bb0e-0b33635a3263
# ╠═ac0f6b9c-abdc-4aa3-bf00-29f63330d219
# ╠═3922c5e8-6d40-4cf7-ae2a-3bcb5ae48d73
# ╠═acf8bb76-bf87-4cf8-b280-e0e13c65fdff
# ╠═1a364c8f-49e3-4c33-b892-8f35106c4afa
# ╠═cb29b45f-d7b1-4507-9531-5d439e4accbb
# ╠═0e703742-3cb5-4701-97b9-6134d2ca30b6
# ╠═e0223248-586e-48a1-9799-d3918f672e92
# ╠═b8006eb5-b99f-4535-81ef-dfe669d0d332
# ╠═db615f71-8858-47da-be2f-8196f7070e6d
# ╠═d896c790-02ab-4f57-801b-cf2eb25aefab
# ╠═41caae29-9181-4cf8-a814-64868a9b8e90
# ╠═3d0cd846-07a7-4910-a2a2-2ef660355ad5
# ╠═74fbb354-2657-475e-8d1e-42d39e255396
# ╠═9f0baca9-46e2-4b20-9c0c-0dbbac36bbf0
# ╠═ffcd720b-134f-4b83-adf7-ab39ab6dfba1
# ╠═de174fc5-5bc0-4b5c-980d-88168f5355d1
# ╠═bea93ba6-a1e5-444a-b419-f0ec290cb61f
# ╠═2639031c-df64-47c4-9d75-4f2ad73291dd
# ╠═15ff286a-a8e4-4c8b-9c73-04450ad6c500
# ╠═b798af09-363d-4800-b37e-a8fe8ce36398
# ╠═1e382321-089d-407f-9d2b-e07dcdb00135
# ╠═445ab169-1024-4238-b3bb-35b9cbdc99e7
# ╠═6212efdb-8f9a-4fc8-824c-0623ac71d2ce
# ╠═a03f024e-d55d-4ffc-aedd-d96f5b539fdd
# ╠═0f849d0c-0a3c-4a70-b437-98b6340ceb58
# ╠═7cc5d85d-27c7-4f4d-b62a-30925a58a98a
# ╠═0d0bde49-85c7-4559-96c5-2602e59b505d
# ╠═9bcf2ca1-b26b-453c-824f-43de4e31c967
# ╠═94310e4e-644a-4b99-ad5f-d708a47ad7f8
# ╠═7dcca23e-dfb0-407c-9e91-47f7e8390990
# ╠═461c7485-2551-4417-9ffc-7cad3f65d2bb
# ╠═f547cf7b-943b-400e-aa4c-f1d60237c984
# ╠═1f5d3f1b-5ac7-4cef-a746-eaadc19c94d0
# ╠═56444e37-03d6-4829-9776-f79dc1a839a2
# ╠═83534334-93f5-422b-8280-a2ca8ccaba4f
# ╠═e346e39a-185b-4fae-9e16-3745bd82a906
# ╠═ebede37d-ee91-4942-86bf-3924b0f8578f
# ╠═ea6bb551-f3af-4647-9620-a3ecf90380b0
# ╠═8c4a1b15-13e8-43a9-9be2-cb614a8c1869
# ╠═44321d83-09e4-4928-ac36-b32a12297320
# ╠═d5b6b49a-9721-4c4c-a735-743f3cd852ab
# ╠═309bf006-c0e1-46a6-9bc8-0834450340d9
# ╠═32d66706-b17d-4c58-8464-1128348ce06a
# ╠═c87fa471-dbb6-4802-8bd0-aae1eafd7fc0
# ╠═4c4144ca-9f36-44f0-b697-7411df2fac6a
# ╠═9fcca185-180f-4478-8fa2-d0cf37e1937b
# ╠═d3b705a7-03ac-4731-8c45-fa95c857d004
# ╠═cb8c0927-cdda-4a9f-aad0-d6b8308ae933
# ╠═935a3442-ed13-4239-962b-9c17dc86e792
# ╠═fbbd6109-ac4e-4676-97b1-e9fef5a3b729
# ╠═52bd54f9-0df6-4c26-849a-64928693fd9a
# ╠═8ff5340c-3ade-4c03-9a72-1c324aa0e882
# ╠═36cf3060-1bfd-4e2e-8009-a9d5bdaf19cb
# ╠═9255fa38-9004-43c9-8817-5657e6cfb2d5
# ╠═6b38d182-a901-4b20-bfdf-9441e4586e7b
# ╠═2cc26b9a-c18b-4321-ad79-a00f8085b525
# ╠═cd53cfd1-ec82-4c3b-93fe-62f0138c6f25
# ╠═d7cd6fac-95f6-4612-bd5a-bb6b839e3f2c
# ╠═34cc67fd-f567-425b-bad0-117b00c85a1d
# ╠═dd25496f-e430-42df-a0b4-143bcb26f9f3
# ╠═fddb0aa7-e722-4f8e-88f3-9255275f8da0
# ╠═e84e8b68-bc2a-412b-b487-3b5163131ca5

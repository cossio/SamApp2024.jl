### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ 953d7355-78c1-4309-9f5f-d0c02abffc51
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ bb8a6eeb-32a3-46dc-982a-d826253291e3
using BioSequences: LongRNA

# ╔═╡ 54d6f1c4-bb92-457c-94ec-e67dd06fd5d0
using DataFrames: DataFrame

# ╔═╡ 80e407e4-f584-4253-b685-28e92b5c9119
using Distributions: Gamma

# ╔═╡ a06fb6da-12ce-40ee-888c-a214914c303f
using Distributions: Poisson

# ╔═╡ 29da1d84-e6bc-4890-9a02-9c49e6505681
using LinearAlgebra: Diagonal

# ╔═╡ d144cfc6-42c0-4e61-a80c-654715407b99
using LinearAlgebra: eigen

# ╔═╡ 9bee752e-9272-46af-8bc0-8c2a2a1807e2
using Makie: @L_str

# ╔═╡ 53a19766-1d66-43ee-b6df-3ef581746a78
using NaNStatistics: nansum, nanmean

# ╔═╡ 8c65263e-1e8e-40d8-b824-915b98a5ce57
using Random: bitrand

# ╔═╡ 75f00c6e-2859-4b9e-8e6e-cb102f2fd071
using Statistics: cor

# ╔═╡ 5d47b2ae-c30a-4d8e-8742-00fabad57f3c
using Statistics: mean

# ╔═╡ 2486a3e8-6e80-4d84-a1b3-e8952f3f7faa
using StatsBase: countmap

# ╔═╡ 5e9f596c-13f9-4e03-a08b-5903c4b28d15
using StatsBase: corspearman

# ╔═╡ 1c6856dd-20cb-449f-abcb-545316b28ab5
md"# Imports"

# ╔═╡ b8d4ed04-59c2-432f-b009-b648f27f89fc
import Makie

# ╔═╡ 949fc97f-e30d-4b6c-aa9a-1ec7f1275df6
import CairoMakie

# ╔═╡ b69fb966-3b81-48dd-bf9f-bc430295665b
import CSV

# ╔═╡ 4c50ae69-a07d-45c0-819d-d49f53e40a57
import HDF5

# ╔═╡ 07ca4371-d346-4e2c-9ca2-f6fcd53bf5a5
import FASTX

# ╔═╡ 444c7349-edab-4828-91c6-d2948fd6e7ee
import Infernal

# ╔═╡ df0195c0-6385-476a-962c-63585aafcc52
import KernelDensity

# ╔═╡ 2e51c390-d4ee-4362-b7cb-409d7cdc8df8
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ f706b593-446e-4554-ae6e-9a9d8bcc07f6
import Rfam

# ╔═╡ b087852c-f67d-4867-b7e2-5fcc2c8abfe2
import SamApp2024

# ╔═╡ d6ef5d92-5dac-4fdc-a5ff-46b2b1d8ed00
import StatsBase

# ╔═╡ 8945829b-7175-492a-8454-358874bb9eac
import Logomaker

# ╔═╡ 90b4f2c0-3d1c-453b-ba96-92e2ac9eb87f
import PlutoUI

# ╔═╡ 571ede6a-ad6b-4710-a17d-9c3c7951d719
PlutoUI.TableOfContents()

# ╔═╡ a768ec77-2f9f-4e18-b488-b0cb756a4b00
md"# General data"

# ╔═╡ 735efa70-2f53-45e0-a301-2ea716ea242a
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ aea6c6c9-faae-4c7e-bf35-b1faace52177
ENV["JULIA_RFAM_DIR"]

# ╔═╡ 84389a19-1f13-4162-a51b-19dcf16ea0c1
# Hallmark sites
_sites = SamApp2024.hallmark_sites_20230507;

# ╔═╡ ee44f513-ffa0-4183-bcf9-fdbcb73dba69
# signifcance threshold for protection scores
_thresh = log(5)

# ╔═╡ 781bdb61-4311-4e5c-95d4-4deb751c7e6e
md"# Load data"

# ╔═╡ aae30131-4654-4f98-83f1-9d0818189785
md"## Load data (Repl.0)"

# ╔═╡ 6781ca9d-85fc-46c9-8059-e62e13edabd2
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ ed515a86-9939-4123-9237-7747817b03d0
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ b222c039-9c93-4bcb-a4cb-bc4219b11f8f
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 392b2a08-597f-47e3-aeba-b45badcd43f5
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ a44fce13-f88e-456b-b870-983503331d16
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 341b5e17-056c-478d-8c3e-264bbf4bfeb0
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ 1294fa8e-d76b-406c-8d07-21e3aac0c59c
rbm_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 6e36ec71-a706-4dab-be55-f2082f4a21c3
inf_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 0afc5cf7-852d-4d05-b3ad-7477bf86c7b9
full_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ 7f93a86c-17e4-4472-a9db-c6b43fce0dca
seed_seqs_rep0 = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 4ff40840-74cb-4ac4-80b7-982c57478969
nat_seqs_rep0 = full_seqs_rep0 ∪ seed_seqs_rep0;

# ╔═╡ 600a0427-ad61-4317-aed9-41ba609c40c0
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 57444eac-7b7b-47c5-8e85-1ca01b4230ce
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 32abd57f-71b4-4188-9f28-0497c5ab7764
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs_rep0, conds_sam_rep0];

# ╔═╡ 6508d862-75c5-4076-9aec-29043a77493e
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ daa2c7e2-2e6a-4a37-b32d-8569f21800bf
x_mg_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :,  conds_mg_rep0]; dim=(1,3));

# ╔═╡ 5f1bb0e2-06e9-4476-8754-6175ed708b2e
x_sam_rep0 = nansum(shape_stats_rep0.shape_log_odds[_sites, :, conds_sam_rep0]; dim=(1,3));

# ╔═╡ a7484a8b-6feb-4ccc-987b-d87f800970f2
_responds_sam_yes_rep0 = (x_mg_rep0 .< -_thresh) .& (x_sam_rep0 .> +_thresh);

# ╔═╡ 430ded66-dc00-411a-9db7-130f2d6555bf
_responds_sam_nop_rep0 = (x_mg_rep0 .> +_thresh) .| (x_sam_rep0 .< -_thresh);

# ╔═╡ e1b3f8ef-9b44-4bec-8610-01145774be1a
aptamer_rbm_energies_rep0 = [
    ismissing(seq) ? missing : 
    RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_045.aligned_sequences
];

# ╔═╡ b0b44ab1-bd99-4239-b93f-3dea79d0c053
md"## Load data (500 seqs)"

# ╔═╡ a94baef6-c552-4b33-9139-0dc61dcf5aa9
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ 3cfc5665-936a-4621-8aa4-820c9a3e2c78
conds_sam_500, conds_mg_500, conds_30C_500 = [1,2], [4], [6];

# ╔═╡ 2e606abc-4ae4-4798-88ff-48153e41038b
bps_reactivities_500 = shape_data_500.shape_reactivities[bps, :, conds_sam_500];

# ╔═╡ d1b68ed8-f33e-4a02-b84f-03130ab96e58
nps_reactivities_500 = shape_data_500.shape_reactivities[nps, :, conds_sam_500];

# ╔═╡ 52a1f8c4-013a-4215-9b12-98f2dc23b2f3
all_reactivities_500 = shape_data_500.shape_reactivities[:, :, conds_sam_500];

# ╔═╡ f2cd1fd6-3b73-4615-9d3f-078dcdf50911
shape_stats_500 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_500,
    paired_reactivities = bps_reactivities_500,
    unpaired_reactivities = nps_reactivities_500,
    all_reactivities = all_reactivities_500,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ d51d3406-0e89-4a36-9591-199e3675d10f
x_mg_500 = nansum(shape_stats_500.shape_log_odds[_sites, :,  conds_mg_500]; dim=(1,3));

# ╔═╡ b9badf2d-b58e-4b12-919c-feb79f911926
x_sam_500 = nansum(shape_stats_500.shape_log_odds[_sites, :, conds_sam_500]; dim=(1,3));

# ╔═╡ d192c0af-8551-4cba-b3c3-e8ea98c1b2a3
_responds_sam_yes_500 = (x_mg_500 .< -_thresh) .& (x_sam_500 .> +_thresh);

# ╔═╡ cb4d2651-8471-4fe6-ae96-413306f9ecc8
_responds_sam_nop_500 = (x_mg_500 .> +_thresh) .| (x_sam_500 .< -_thresh);

# ╔═╡ 1639c8dc-e96d-4225-bcda-b838ad391d8b
aptamer_rbm_energies_500 = [
    ismissing(seq) ? missing : 
    RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in shape_data_500.aligned_sequences
];

# ╔═╡ 67ee804e-8877-4fc1-b2b5-a213bce545cd
md"# Make table"

# ╔═╡ 7e76da0f-91ca-4c5b-bf1d-ca10a2a6c5a8
_responsive_sam_rep0 = ifelse.(_responds_sam_yes_rep0, "Responsive", ifelse.(_responds_sam_nop_rep0, "Non-responsive", "Inconclusive"))

# ╔═╡ bfe80af1-a0ad-4702-9571-869aec104c28
_responsive_sam_500 = ifelse.(_responds_sam_yes_500, "Responsive", ifelse.(_responds_sam_nop_500, "Non-responsive", "Inconclusive"))

# ╔═╡ b5b6abab-40b0-430e-a012-59af58325a4d
ss = SamApp2024.RF00162_sites_annotated_secondary_structure()

# ╔═╡ 400df1dd-979a-45f4-9884-a6d9c012bd32
seq_groups_dfs = SamApp2024.artifact_load_sequencing_groups_2024_11_27()

# ╔═╡ f1d5d25d-aa49-4255-9b99-413356541622
md"# Add columns interesting for DMS"

# ╔═╡ 1cfc835b-679a-4966-9088-275ef43a0159
Extended_Hallmark_sites = [
	1:8; 101:108; # P1
	10:11; # SAM contact
    25:28; 77:80; # Pseudoknot
    34:37; # Kink-turn
    73:75; # A-minor
    46; 47; # SAM contact (bulge)
    24; 76; 100; # Base-triple
];

# ╔═╡ 51bae40e-9128-4675-b021-cda716620d05
sort([seq_groups_dfs["GP6-Synthetic-Set2-primer1"].name; seq_groups_dfs["GP7-Synthetic-Set2-primer2"].name; seq_groups_dfs["GP8-Synthetic-Set2-primer3"].name; seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name]) == ["APSAM-S2-" * lpad(n - 1, 3, "0") for n = 1:450]

# ╔═╡ 0840f95e-092f-4d89-a697-0cfc38e810cf
allunique([seq_groups_dfs["GP6-Synthetic-Set2-primer1"].sequence; seq_groups_dfs["GP7-Synthetic-Set2-primer2"].sequence; seq_groups_dfs["GP8-Synthetic-Set2-primer3"].sequence; seq_groups_dfs["GP9-Synthetic-Set2-primer5"].sequence])

# ╔═╡ dd49243c-3bbb-40fa-ac72-308c4ddebb57
group_names = sort(collect(keys(seq_groups_dfs)))

# ╔═╡ f6bb72e2-267f-45e6-b2f9-905bda307b06
begin
	df = DataFrame(;
		aptamer_names = [shape_data_rep0.aptamer_names; ["APSAM-S2-" * lpad(n - 1, 3, "0") for n = 1:500][shape_data_500.aptamer_origin .!= "Infrared"]],
	    aligned_sequences = [[ismissing(seq) ? missing : string(seq) for seq = shape_data_rep0.aligned_sequences]; string.(shape_data_500.aligned_sequences[shape_data_500.aptamer_origin .!= "Infrared"])],
	    aptamer_origin = [shape_data_rep0.aptamer_origin; shape_data_500.aptamer_origin[shape_data_500.aptamer_origin .!= "Infrared"]],
		experiment = [fill("Experiment_1", length(shape_data_rep0.aligned_sequences)); fill("Experiment_2", length(shape_data_500.aligned_sequences))[shape_data_500.aptamer_origin .!= "Infrared"]],
	    responsive = [_responsive_sam_rep0; _responsive_sam_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		RBM_score = [-aptamer_rbm_energies_rep0; -aptamer_rbm_energies_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		Protect_Score_Hallmark_Mg = [x_mg_rep0; x_mg_500[shape_data_500.aptamer_origin .!= "Infrared"]],
		Protect_Score_Hallmark_SAM = [x_sam_rep0; x_sam_500[shape_data_500.aptamer_origin .!= "Infrared"]]
	)

	RF00162_to_probed_distances = [ismissing(s2) ? missing : SamApp2024.hamming(s1, LongRNA{4}(s2)) for s1 = SamApp2024.rfam_RF00162_hits(), s2 = df.aligned_sequences]
	
	df.P4_length = [ismissing(seq) ? missing : length(replace(seq[ss.p4], '-' => "")) for seq = df.aligned_sequences]
	df.min_dist_to_natural = dropdims(minimum(RF00162_to_probed_distances; dims=1); dims=1)
	df.A_count = [ismissing(seq) ? missing : count(==('A'), seq) for seq = df.aligned_sequences]
	df.C_count = [ismissing(seq) ? missing : count(==('C'), seq) for seq = df.aligned_sequences]
	df.A_count_in_extended_hallamark = [ismissing(seq) ? missing : count(==('A'), seq[Extended_Hallmark_sites]) for seq = df.aligned_sequences]
	df.C_count_in_extended_hallamark = [ismissing(seq) ? missing : count(==('C'), seq[Extended_Hallmark_sites]) for seq = df.aligned_sequences]
	df.seq_lenth = [ismissing(seq) ? missing : length(replace(seq, '-' => "")) for seq = df.aligned_sequences]

	df.sequencing_group = [group_names[only(findall([n ∈ seq_groups_dfs[k].name for k = group_names]))] for n = df.aptamer_names]

	@assert count(df.experiment .== "Experiment_1") == length(shape_data_045.aligned_sequences) == 306
	for (seq1, seq2) = zip(df.aligned_sequences[df.experiment .== "Experiment_1"], shape_data_045.aligned_sequences)
		@assert ismissing(seq1) && ismissing(seq2) || seq1 == seq2
	end
	@assert df.experiment == [fill("Experiment_1", 306); fill("Experiment_2", 450)]
	
	for df = values(seq_groups_dfs)
		df.rna_seq = [replace(seq, 'T' => 'U') for seq = df.sequence] 
	end

	#df.read_depth_M = [nanmean(shape_data_rep0.shape_M_depth; dim=(1,3)); nanmean(shape_data_500.shape_M_depth[:, 1:450, :]; dim=(1,3))]
	df.read_depth_M = [
		[shape_data_rep0.shape_M_depth[:, n, :] for n = axes(shape_data_rep0.shape_M_depth, 2)];
		[shape_data_500.shape_M_depth[:, n, :] for n = 1:450]
	]
end

# ╔═╡ 9a8605cc-99a1-4fcc-9e8d-dd8a61fd26ca
[replace(string(seq), '-' => "") for seq = df.aligned_sequences[df.experiment .== "Experiment_2"]][
	indexin(
		[
			seq_groups_dfs["GP6-Synthetic-Set2-primer1"].name;
			seq_groups_dfs["GP7-Synthetic-Set2-primer2"].name;
			seq_groups_dfs["GP8-Synthetic-Set2-primer3"].name;
			seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name
		],
		df.aptamer_names[df.experiment .== "Experiment_2"]
)] == [seq_groups_dfs["GP6-Synthetic-Set2-primer1"].rna_seq; seq_groups_dfs["GP7-Synthetic-Set2-primer2"].rna_seq; seq_groups_dfs["GP8-Synthetic-Set2-primer3"].rna_seq; seq_groups_dfs["GP9-Synthetic-Set2-primer5"].rna_seq]

# ╔═╡ 0e568f56-d77b-4b67-9f45-5ccbe4b21841
indexin(
		[
			seq_groups_dfs["GP6-Synthetic-Set2-primer1"].name;
			seq_groups_dfs["GP7-Synthetic-Set2-primer2"].name;
			seq_groups_dfs["GP8-Synthetic-Set2-primer3"].name;
			seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name
		],
		df.aptamer_names[df.experiment .== "Experiment_2"]
	)

# ╔═╡ 6756bc38-4381-4d65-8af4-eddeaeb8cdc8
sort(df.aptamer_names[df.experiment .== "Experiment_2"]) == sort([
	seq_groups_dfs["GP6-Synthetic-Set2-primer1"].name;
	seq_groups_dfs["GP7-Synthetic-Set2-primer2"].name;
	seq_groups_dfs["GP8-Synthetic-Set2-primer3"].name;
	seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name
])

# ╔═╡ d80fdd97-e1c5-48aa-94c1-ed8dcf1b2b59
indexin(seq_groups_dfs["GP3-Natural-primer3"].name, shape_data_045.aptamer_names)

# ╔═╡ 79421a94-979a-48fc-939c-079202bf115c
seq_groups_dfs["GP3-Natural-primer3"].name

# ╔═╡ 528d061b-3147-47a0-be06-a8e673edc1d8
df.aptamer_names[df.sequencing_group .== "GP3-Natural-primer3"] == seq_groups_dfs["GP3-Natural-primer3"].name

# ╔═╡ 50f29fb4-00f6-4ead-970c-04f7aad653aa
seq_groups_dfs["GP1-Natural-primer1"].name

# ╔═╡ fec7d3da-3b86-43a5-886d-4e621155901d


# ╔═╡ 18582db7-f6a2-4892-9a9e-b31d31c05dbd
nanmean(shape_data_rep0.shape_M_depth[:, df.sequencing_group[1:306] .== "GP3-Natural-primer3", :])

# ╔═╡ 09a381b2-4564-4499-887b-14b7bdbaef00
nanmean(shape_data_rep0.shape_M_depth[:, indexin(seq_groups_dfs["GP3-Natural-primer3"].name, shape_data_045.aptamer_names), :])

# ╔═╡ f3b518da-4ffd-45cb-85ed-2d85cfdc2205
cor(
	[nanmean(mapreduce(vec, vcat, df.read_depth_M[df.sequencing_group .== gr])) for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))],
	[mean(df.responsive[df.sequencing_group .== gr] .== "Inconclusive") for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))]
)

# ╔═╡ afa683f5-4bc3-4bfd-8e73-ba3e4c911368
let fig = Makie.Figure()
	primer_colors = (:orange, :purple, :red, :blue, :teal)
	#markers = (:utriangle, :xcross, :dtriangle, :circle, :cross)
	markers = (:circle, :circle, :circle, :circle, :circle)
	
	ax = Makie.Axis(fig[1,1]; width=230, height=230, xlabel="Read depth (M)", ylabel="Inconclusives / Total")
	for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
		primer = parse(Int, last(gr))
		
		read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[df.sequencing_group .== gr]))
		inconclusive_rate = mean(df.responsive[df.sequencing_group .== gr] .== "Inconclusive")
		if gr_n ≤ 5
			Makie.scatter!(ax, read_depth, inconclusive_rate; color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
		elseif gr_n < 9 # skip last one which has a single sequence
			Makie.scatter!(ax, read_depth, inconclusive_rate; color=primer_colors[primer], marker=markers[primer], markersize=15)
		end
	end
	Makie.xlims!(ax, 0, 3.3e4)
	Makie.ylims!(ax, -0.05, 0.65)

	fig[1,2] = Makie.Legend(fig, ax, "Primers", framevisible = false)

	for (col, (origin, title)) = enumerate(zip([["RF00162_full30", "RF00162_seed70"], ["RF00162_syn_rbm", "rbm"], ["RF00162_syn_inf", "infernal"]], ["Naturals", "RBM", "CM"]))
		ax = Makie.Axis(fig[1,3][1,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. / Total", title, xticks=[1e4, 3e4])
		for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
			_flag = (df.sequencing_group .== gr) .&& (df.aptamer_origin .∈ Ref(origin))
			primer = parse(Int, last(gr))
			println("$gr Primer$primer origin $origin ", count(_flag))
			if any(_flag)
				read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[_flag]))
				resp_rate = mean(df.responsive[_flag] .== "Responsive")
				if gr_n ≤ 5
					Makie.scatter!(ax, read_depth, resp_rate; color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
				elseif gr_n < 9 # skip last one which has a single sequence
					Makie.scatter!(ax, read_depth, resp_rate; color=primer_colors[primer], marker=markers[primer], markersize=15)
				end
			end
		end
		Makie.xlims!(ax, 0, 3.3e4)
		Makie.ylims!(ax, -0.05, 0.75)

		ax = Makie.Axis(fig[1,3][2,col]; width=100, height=100, xlabel="Read depth (M)", ylabel="Resp. rate", xticks=[1e4, 3e4])
		for (gr_n, gr) = enumerate(sort(unique(df.sequencing_group)))
			_flag = (df.sequencing_group .== gr) .&& (df.aptamer_origin .∈ Ref(origin))
			if any(_flag)
				primer = parse(Int, last(gr))
				read_depth = nanmean(mapreduce(vec, vcat, df.read_depth_M[_flag]))
				resp_rate = mean(df.responsive[_flag] .== "Responsive")
				inconclusive_rate = mean(df.responsive[df.sequencing_group .== gr] .== "Inconclusive")
				if gr_n ≤ 5
					Makie.scatter!(ax, read_depth, resp_rate / (1 - inconclusive_rate); color=primer_colors[primer], marker=markers[primer], label="Primer $primer", markersize=15)
				elseif gr_n < 9 # skip last one which has a single sequence
					Makie.scatter!(ax, read_depth, resp_rate / (1 - inconclusive_rate); color=primer_colors[primer], marker=markers[primer], markersize=15)
				end
			end
		end
		Makie.xlims!(ax, 0, 3.3e4)
		Makie.ylims!(ax, -0.05, 0.75)

	end
	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ f35de085-a396-46b7-9461-de60d36b0068
df[(df.aptamer_origin .== "infernal" .|| df.aptamer_origin .== "RF00162_syn_inf"), :]

# ╔═╡ 321e2dcc-09c6-4f8e-aece-8a79488fd291
StatsBase.countmap(df[(df.aptamer_origin .== "infernal" .|| df.aptamer_origin .== "RF00162_syn_inf"), :].sequencing_group)

# ╔═╡ 8d9174fd-e8cb-479a-b567-084d6273f9da
df.responsive

# ╔═╡ f9314f91-6986-495f-8797-3a6b5588e4ad
unique(df.responsive)

# ╔═╡ c89ab7a9-a35c-474b-b9db-70d0a2c226c7
unique(df.sequencing_group)

# ╔═╡ 7452003d-b3f3-4e79-ad2d-638f53b0f91c
unique(df.aptamer_origin)

# ╔═╡ 61359389-76ca-4413-8fde-e7a5c8b2a23f
df.aptamer_origin[df.aptamer_origin .∈ Ref(["RF00162_syn_inf", "infernal"])]

# ╔═╡ afdba201-d750-42b5-9ba2-6bd1946d3a21
sort(unique(df.sequencing_group))

# ╔═╡ 8ea8987c-5b69-495e-af4f-e71f86a648da
names_500 = ["APSAM-S2-" * lpad(n - 1, 3, "0") for n = 1:500][shape_data_500.aptamer_origin .!= "Infrared"]

# ╔═╡ 3ff21fe5-49a0-4fc4-ad25-ca0fee168e88
read_depths = merge(
	Dict(gr => nanmean(shape_data_500.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :])
		for gr = ["GP6-Synthetic-Set2-primer1", "GP7-Synthetic-Set2-primer2", "GP8-Synthetic-Set2-primer3", "GP9-Synthetic-Set2-primer5"]),
	Dict(gr => nanmean(shape_data_rep0.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :])
		for gr = ["GP1-Natural-primer1", "GP2-Natural-primer2", "GP3-Natural-primer3", "GP4-Synthetic-Set1-primer4", "GP5-Synthetic-Set1-primer5"])
)

# ╔═╡ abfae16d-d438-40d1-aef1-10cd407cbc70
begin
	df_groups = DataFrame()
	for g = group_names
		group_name = g
		total_sequences = length(df.responsive[df.sequencing_group .== g, :])
		A_rate = sum(skipmissing(df.A_count[df.sequencing_group .== g, :])) / sum(skipmissing(df.seq_lenth[df.sequencing_group .== g]))
		C_rate = sum(skipmissing(df.C_count[df.sequencing_group .== g, :])) / sum(skipmissing(df.seq_lenth[df.sequencing_group .== g]))
		A_hallmark_ex_rate = sum(skipmissing(df.A_count_in_extended_hallamark[df.sequencing_group .== g, :])) / (total_sequences * length(Extended_Hallmark_sites))
		C_hallmark_ex_rate = sum(skipmissing(df.C_count_in_extended_hallamark[df.sequencing_group .== g, :])) / (total_sequences * length(Extended_Hallmark_sites))
		responsives = sum(df.responsive[df.sequencing_group .== g, :] .== "Responsive")
		inconclusives = sum(df.responsive[df.sequencing_group .== g, :] .== "Inconclusive")
		RBM_score_avg = mean(skipmissing(df.RBM_score[df.sequencing_group .== g]))
		RBM_score_min = minimum(skipmissing(df.RBM_score[df.sequencing_group .== g]))
		P4_len_avg = mean(skipmissing(df.P4_length[df.sequencing_group .== g]))
		P4_len_min = minimum(skipmissing(df.P4_length[df.sequencing_group .== g]))
		noP4_count = count(skipmissing(df.P4_length[df.sequencing_group .== g]) .≤ 1)
		min_dist_to_nat_avg = mean(skipmissing(df.min_dist_to_natural[df.sequencing_group .== g]))
		min_dist_to_nat_max = maximum(skipmissing(df.min_dist_to_natural[df.sequencing_group .== g]))
		response_rate = responsives ./ total_sequences
		inconclusive_rate = inconclusives ./ total_sequences
		push!(
			df_groups,
			(;
			 group_name, total_sequences, responsives, inconclusives, response_rate, inconclusive_rate, A_rate, C_rate, A_hallmark_ex_rate, C_hallmark_ex_rate, RBM_score_avg, RBM_score_min, P4_len_avg, 	P4_len_min, noP4_count, min_dist_to_nat_avg, min_dist_to_nat_max
			)
		)
	end

	df_groups.M_read_depths = [read_depths[gr] for gr = df_groups.group_name]
	#show(df_groups; allcols=true)
	df_groups
end

# ╔═╡ c0ee24d7-6560-496b-9c16-181d8577a422
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=230, height=230, xlabel="Read depth (M)", ylabel="Inconclusive rate")

	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], df_groups.inconclusive_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], df_groups.inconclusive_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], df_groups.inconclusive_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.inconclusive_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5]], df_groups.inconclusive_rate[[5]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.3e4)
	Makie.ylims!(ax, -0.05, 0.65)
		
	ax = Makie.Axis(fig[1,3]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate", title="Natural")
	Makie.scatter!(ax, df_groups.M_read_depths[[1]], df_groups.response_rate[[1]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2]], df_groups.response_rate[[2]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3]], df_groups.response_rate[[3]]; color=:red, label="Primer 3", markersize=15)
	Makie.xlims!(ax, 0, 3.3e4)
	Makie.ylims!(ax, -0.05, 0.65)

	ax = Makie.Axis(fig[1,4]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate", title="Artificial")
	Makie.scatter!(ax, df_groups.M_read_depths[[6]], df_groups.response_rate[[6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[7]], df_groups.response_rate[[7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[8]], df_groups.response_rate[[8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.response_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5]], df_groups.response_rate[[5]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.3e4)
	Makie.ylims!(ax, -0.05, 0.65)

	fig[1,2] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 8d003cc1-dbb8-4ae8-97d4-b74174c63b14
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], df_groups.response_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], df_groups.response_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], df_groups.response_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.response_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5]], df_groups.response_rate[[5]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	ax = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], df_groups.inconclusive_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], df_groups.inconclusive_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], df_groups.inconclusive_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.inconclusive_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5]], df_groups.inconclusive_rate[[5]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	ax = Makie.Axis(fig[1,3]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1]], df_groups.response_rate[[1]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2]], df_groups.response_rate[[2]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3]], df_groups.response_rate[[3]]; color=:red, label="Primer 3", markersize=15)
	Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	ax = Makie.Axis(fig[1,4]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1]], df_groups.inconclusive_rate[[1]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2]], df_groups.inconclusive_rate[[2]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3]], df_groups.inconclusive_rate[[3]]; color=:red, label="Primer 3", markersize=15)
	Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	fig[1,5] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 399072a0-2b27-4175-9cc8-479950d69215
let fig = Makie.Figure()
	
	renorm_response = df_groups.response_rate ./ (1 .- df_groups.inconclusive_rate)
	
	ax = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate / (1 - Incl. rate)")
	Makie.scatter!(ax, df_groups.M_read_depths[[1,6]], renorm_response[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2,7]], renorm_response[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3,8]], renorm_response[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], renorm_response[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5]], renorm_response[[5]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)

	ax = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="RBM score (avg.)", ylabel="Response rate / (1 - Incl. rate)")
	Makie.scatter!(ax, df_groups.RBM_score_avg[[1,6]], renorm_response[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[2,7]], renorm_response[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[3,8]], renorm_response[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[4]], renorm_response[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[5]], renorm_response[[5]]; color=:teal, label="Primer 5", markersize=15)
	#Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)

	ax = Makie.Axis(fig[1,3]; width=200, height=200, xlabel="RBM score (avg.)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.RBM_score_avg[[1,6]], df_groups.response_rate[[1,6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[2,7]], df_groups.response_rate[[2,7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[3,8]], df_groups.response_rate[[3,8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[4]], df_groups.response_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.RBM_score_avg[[5]], df_groups.response_rate[[5]]; color=:teal, label="Primer 5", markersize=15)
	#Makie.xlims!(ax, 0, 3.2e4)
	Makie.ylims!(ax, -0.05, 0.65)

	fig[1,5] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 824b1bf2-872e-40c8-9eb9-49adcc8cb9a9
df_groups.RBM_score_avg

# ╔═╡ bcb90e0e-4797-4160-9fb8-6d433bae5aaa
let fig = Makie.Figure() # natural only
	ax = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1]], df_groups.response_rate[[1]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2]], df_groups.response_rate[[2]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3]], df_groups.response_rate[[3]]; color=:red, label="Primer 3", markersize=15)
	Makie.xlims!(ax, 0, 3e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	ax = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[1]], df_groups.inconclusive_rate[[1]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[2]], df_groups.inconclusive_rate[[2]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[3]], df_groups.inconclusive_rate[[3]]; color=:red, label="Primer 3", markersize=15)
	Makie.xlims!(ax, 0, 3e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	fig[1,3] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ e142d322-2825-4dae-898b-76b0613be0d9
let fig = Makie.Figure() # synhtetic only
	ax = Makie.Axis(fig[1,1]; width=200, height=200, xlabel="Read depth (M)", ylabel="Response rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[6]], df_groups.response_rate[[6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[7]], df_groups.response_rate[[7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[8]], df_groups.response_rate[[8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.response_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5,9]], df_groups.response_rate[[5,9]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.1e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	ax = Makie.Axis(fig[1,2]; width=200, height=200, xlabel="Read depth (M)", ylabel="Inconclusive rate")
	Makie.scatter!(ax, df_groups.M_read_depths[[6]], df_groups.inconclusive_rate[[6]]; color=:orange, label="Primer 1", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[7]], df_groups.inconclusive_rate[[7]]; color=:purple, label="Primer 2", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[8]], df_groups.inconclusive_rate[[8]]; color=:red, label="Primer 3", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[4]], df_groups.inconclusive_rate[[4]]; color=:blue, label="Primer 4", markersize=15)
	Makie.scatter!(ax, df_groups.M_read_depths[[5,9]], df_groups.inconclusive_rate[[5,9]]; color=:teal, label="Primer 5", markersize=15)
	Makie.xlims!(ax, 0, 3.1e4)
	Makie.ylims!(ax, -0.05, 0.65)
	
	fig[1,3] = Makie.Legend(fig, ax, "Primers", framevisible = false)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 66cc828d-53c6-4fac-a564-25c1c19ed45d
df_groups.inconclusive_rate[[5,9]]

# ╔═╡ 900d9a2b-d9ab-4730-9163-474cddc25689
incl_rate = (df_groups.total_sequences - df_groups.responsives) ./ df_groups.total_sequences

# ╔═╡ 53f5d185-ade0-4a4c-a77f-3dd08aecb7ac
cor(log.(df_groups.M_read_depths), df_groups.response_rate)

# ╔═╡ 72d66520-8b99-4f49-b42d-24f6681a656e
CSV.write(tempname(), df_groups)

# ╔═╡ 31b5c866-fb43-42c6-a020-3fc027cf9e4f
df[df.sequencing_group .== "GP4-Synthetic-Set1-primer4", :]

# ╔═╡ 891a377f-a6bc-46f2-afdb-32d587e06cad
df[df.sequencing_group .== "GP5-Synthetic-Set1-primer5", :]

# ╔═╡ 6162ba9b-038c-4e7e-9c30-8045c32dede5
md"# Determine groups of the example aptamers in main figures"

# ╔═╡ 7cdd22de-8c5f-4f11-953a-fc5c12701b5a
# Second example from Fig.5 (Deltaproteobacteria bacterium)
df[df.aptamer_names .== "APSAMN7", :]

# ╔═╡ 67acffd5-6850-4d61-9745-6595c471aba4
# Examples from Fig.7
df[df.aptamer_names .== shape_data_045.aptamer_names[299], :]

# ╔═╡ 94b1b7f0-c636-46bd-ae67-fbdeebcf5c11
# Examples from Fig.7
df[df.aptamer_names .== shape_data_045.aptamer_names[207], :]

# ╔═╡ 737d32e6-4889-4704-bd0a-2a19d8e3780a
df[df.aptamer_names .== shape_data_045.aptamer_names[299], :].aligned_sequences

# ╔═╡ c9e418b6-00f2-4e76-a47b-67bc42b93467
# Example from Fig.9 (no P4)
df[df.aligned_sequences .=== string(shape_data_500.aligned_sequences[116]), :]

# ╔═╡ fc6019fd-13e3-453b-8ab5-161522b09a71
# Example from Fig.9 (distant)
df[df.aligned_sequences .=== string(shape_data_500.aligned_sequences[284]), :]

# ╔═╡ c62aa2db-f222-4e74-b1cb-9be741b60bc8
string(shape_data_500.aligned_sequences[284])

# ╔═╡ a2dee926-4434-4983-b70b-a449e512ad00
only(df.aligned_sequences[df.aligned_sequences .=== string(shape_data_500.aligned_sequences[284])])

# ╔═╡ Cell order:
# ╠═1c6856dd-20cb-449f-abcb-545316b28ab5
# ╠═953d7355-78c1-4309-9f5f-d0c02abffc51
# ╠═b8d4ed04-59c2-432f-b009-b648f27f89fc
# ╠═949fc97f-e30d-4b6c-aa9a-1ec7f1275df6
# ╠═b69fb966-3b81-48dd-bf9f-bc430295665b
# ╠═4c50ae69-a07d-45c0-819d-d49f53e40a57
# ╠═07ca4371-d346-4e2c-9ca2-f6fcd53bf5a5
# ╠═444c7349-edab-4828-91c6-d2948fd6e7ee
# ╠═df0195c0-6385-476a-962c-63585aafcc52
# ╠═2e51c390-d4ee-4362-b7cb-409d7cdc8df8
# ╠═f706b593-446e-4554-ae6e-9a9d8bcc07f6
# ╠═b087852c-f67d-4867-b7e2-5fcc2c8abfe2
# ╠═d6ef5d92-5dac-4fdc-a5ff-46b2b1d8ed00
# ╠═8945829b-7175-492a-8454-358874bb9eac
# ╠═90b4f2c0-3d1c-453b-ba96-92e2ac9eb87f
# ╠═bb8a6eeb-32a3-46dc-982a-d826253291e3
# ╠═54d6f1c4-bb92-457c-94ec-e67dd06fd5d0
# ╠═80e407e4-f584-4253-b685-28e92b5c9119
# ╠═a06fb6da-12ce-40ee-888c-a214914c303f
# ╠═29da1d84-e6bc-4890-9a02-9c49e6505681
# ╠═d144cfc6-42c0-4e61-a80c-654715407b99
# ╠═9bee752e-9272-46af-8bc0-8c2a2a1807e2
# ╠═53a19766-1d66-43ee-b6df-3ef581746a78
# ╠═8c65263e-1e8e-40d8-b824-915b98a5ce57
# ╠═75f00c6e-2859-4b9e-8e6e-cb102f2fd071
# ╠═5d47b2ae-c30a-4d8e-8742-00fabad57f3c
# ╠═2486a3e8-6e80-4d84-a1b3-e8952f3f7faa
# ╠═5e9f596c-13f9-4e03-a08b-5903c4b28d15
# ╠═571ede6a-ad6b-4710-a17d-9c3c7951d719
# ╠═a768ec77-2f9f-4e18-b488-b0cb756a4b00
# ╠═735efa70-2f53-45e0-a301-2ea716ea242a
# ╠═aea6c6c9-faae-4c7e-bf35-b1faace52177
# ╠═84389a19-1f13-4162-a51b-19dcf16ea0c1
# ╠═ee44f513-ffa0-4183-bcf9-fdbcb73dba69
# ╠═781bdb61-4311-4e5c-95d4-4deb751c7e6e
# ╠═aae30131-4654-4f98-83f1-9d0818189785
# ╠═6781ca9d-85fc-46c9-8059-e62e13edabd2
# ╠═ed515a86-9939-4123-9237-7747817b03d0
# ╠═b222c039-9c93-4bcb-a4cb-bc4219b11f8f
# ╠═392b2a08-597f-47e3-aeba-b45badcd43f5
# ╠═a44fce13-f88e-456b-b870-983503331d16
# ╠═341b5e17-056c-478d-8c3e-264bbf4bfeb0
# ╠═1294fa8e-d76b-406c-8d07-21e3aac0c59c
# ╠═6e36ec71-a706-4dab-be55-f2082f4a21c3
# ╠═0afc5cf7-852d-4d05-b3ad-7477bf86c7b9
# ╠═7f93a86c-17e4-4472-a9db-c6b43fce0dca
# ╠═4ff40840-74cb-4ac4-80b7-982c57478969
# ╠═600a0427-ad61-4317-aed9-41ba609c40c0
# ╠═57444eac-7b7b-47c5-8e85-1ca01b4230ce
# ╠═32abd57f-71b4-4188-9f28-0497c5ab7764
# ╠═6508d862-75c5-4076-9aec-29043a77493e
# ╠═daa2c7e2-2e6a-4a37-b32d-8569f21800bf
# ╠═5f1bb0e2-06e9-4476-8754-6175ed708b2e
# ╠═a7484a8b-6feb-4ccc-987b-d87f800970f2
# ╠═430ded66-dc00-411a-9db7-130f2d6555bf
# ╠═e1b3f8ef-9b44-4bec-8610-01145774be1a
# ╠═b0b44ab1-bd99-4239-b93f-3dea79d0c053
# ╠═a94baef6-c552-4b33-9139-0dc61dcf5aa9
# ╠═3cfc5665-936a-4621-8aa4-820c9a3e2c78
# ╠═2e606abc-4ae4-4798-88ff-48153e41038b
# ╠═d1b68ed8-f33e-4a02-b84f-03130ab96e58
# ╠═52a1f8c4-013a-4215-9b12-98f2dc23b2f3
# ╠═f2cd1fd6-3b73-4615-9d3f-078dcdf50911
# ╠═d51d3406-0e89-4a36-9591-199e3675d10f
# ╠═b9badf2d-b58e-4b12-919c-feb79f911926
# ╠═d192c0af-8551-4cba-b3c3-e8ea98c1b2a3
# ╠═cb4d2651-8471-4fe6-ae96-413306f9ecc8
# ╠═1639c8dc-e96d-4225-bcda-b838ad391d8b
# ╠═67ee804e-8877-4fc1-b2b5-a213bce545cd
# ╠═7e76da0f-91ca-4c5b-bf1d-ca10a2a6c5a8
# ╠═bfe80af1-a0ad-4702-9571-869aec104c28
# ╠═b5b6abab-40b0-430e-a012-59af58325a4d
# ╠═400df1dd-979a-45f4-9884-a6d9c012bd32
# ╠═f6bb72e2-267f-45e6-b2f9-905bda307b06
# ╠═f1d5d25d-aa49-4255-9b99-413356541622
# ╠═1cfc835b-679a-4966-9088-275ef43a0159
# ╠═51bae40e-9128-4675-b021-cda716620d05
# ╠═0840f95e-092f-4d89-a697-0cfc38e810cf
# ╠═9a8605cc-99a1-4fcc-9e8d-dd8a61fd26ca
# ╠═0e568f56-d77b-4b67-9f45-5ccbe4b21841
# ╠═6756bc38-4381-4d65-8af4-eddeaeb8cdc8
# ╠═dd49243c-3bbb-40fa-ac72-308c4ddebb57
# ╠═3ff21fe5-49a0-4fc4-ad25-ca0fee168e88
# ╠═d80fdd97-e1c5-48aa-94c1-ed8dcf1b2b59
# ╠═79421a94-979a-48fc-939c-079202bf115c
# ╠═528d061b-3147-47a0-be06-a8e673edc1d8
# ╠═50f29fb4-00f6-4ead-970c-04f7aad653aa
# ╠═fec7d3da-3b86-43a5-886d-4e621155901d
# ╠═18582db7-f6a2-4892-9a9e-b31d31c05dbd
# ╠═09a381b2-4564-4499-887b-14b7bdbaef00
# ╠═abfae16d-d438-40d1-aef1-10cd407cbc70
# ╠═c0ee24d7-6560-496b-9c16-181d8577a422
# ╠═f3b518da-4ffd-45cb-85ed-2d85cfdc2205
# ╠═afa683f5-4bc3-4bfd-8e73-ba3e4c911368
# ╠═f35de085-a396-46b7-9461-de60d36b0068
# ╠═321e2dcc-09c6-4f8e-aece-8a79488fd291
# ╠═8d9174fd-e8cb-479a-b567-084d6273f9da
# ╠═f9314f91-6986-495f-8797-3a6b5588e4ad
# ╠═c89ab7a9-a35c-474b-b9db-70d0a2c226c7
# ╠═7452003d-b3f3-4e79-ad2d-638f53b0f91c
# ╠═61359389-76ca-4413-8fde-e7a5c8b2a23f
# ╠═afdba201-d750-42b5-9ba2-6bd1946d3a21
# ╠═8d003cc1-dbb8-4ae8-97d4-b74174c63b14
# ╠═399072a0-2b27-4175-9cc8-479950d69215
# ╠═824b1bf2-872e-40c8-9eb9-49adcc8cb9a9
# ╠═bcb90e0e-4797-4160-9fb8-6d433bae5aaa
# ╠═e142d322-2825-4dae-898b-76b0613be0d9
# ╠═66cc828d-53c6-4fac-a564-25c1c19ed45d
# ╠═900d9a2b-d9ab-4730-9163-474cddc25689
# ╠═53f5d185-ade0-4a4c-a77f-3dd08aecb7ac
# ╠═8ea8987c-5b69-495e-af4f-e71f86a648da
# ╠═72d66520-8b99-4f49-b42d-24f6681a656e
# ╠═31b5c866-fb43-42c6-a020-3fc027cf9e4f
# ╠═891a377f-a6bc-46f2-afdb-32d587e06cad
# ╠═6162ba9b-038c-4e7e-9c30-8045c32dede5
# ╠═7cdd22de-8c5f-4f11-953a-fc5c12701b5a
# ╠═67acffd5-6850-4d61-9745-6595c471aba4
# ╠═94b1b7f0-c636-46bd-ae67-fbdeebcf5c11
# ╠═737d32e6-4889-4704-bd0a-2a19d8e3780a
# ╠═c9e418b6-00f2-4e76-a47b-67bc42b93467
# ╠═fc6019fd-13e3-453b-8ab5-161522b09a71
# ╠═c62aa2db-f222-4e74-b1cb-9be741b60bc8
# ╠═a2dee926-4434-4983-b70b-a449e512ad00

### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ 0fe3899a-1da4-4c5e-bec0-ba1fdeefbb24
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 845b6fe0-93c1-400e-97d3-ef78b27ed175
using DataFrames: DataFrame

# ╔═╡ 2386a288-a090-46ae-9b28-a1c90dbf5c25
using BioSequences: LongRNA

# ╔═╡ 9f34f37d-74b3-48c0-96eb-4f0d2c863704
using Makie: @L_str

# ╔═╡ 35675b3d-e73e-42e2-ac2d-7202483872cc
using Statistics: mean, cor

# ╔═╡ 2b3a2e14-9f10-4c55-b301-447ce730b15b
using NaNStatistics: nanmean, nanstd, nansum

# ╔═╡ 1ad1f988-4794-11f0-0ccb-e3e31a210177
md"# Imports"

# ╔═╡ e15eacbd-c78a-406c-9fbc-046d81ac0d03
import Makie, CairoMakie, PlutoUI

# ╔═╡ 14c7987f-8487-4e9c-9998-36714829fcad
import FASTX, Infernal, Rfam, ViennaRNA

# ╔═╡ 0e41f586-a85b-496f-9503-f456286a9d9c
import SamApp2024

# ╔═╡ 32433288-8243-4143-b301-69a30f60f3d2
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ a056c80e-f046-4fed-af48-49abec41474b
PlutoUI.TableOfContents()

# ╔═╡ f9c24c3e-5697-4ef0-ab9d-932b2d115406
md"# Data"

# ╔═╡ da043e8d-e05b-4861-baef-9ef666431cf7
_sites = SamApp2024.hallmark_sites_20230507

# ╔═╡ c15c31a3-3723-4daf-872e-b2f9806530ce
bps, nps, pks = SamApp2024.RF00162_sites_paired()

# ╔═╡ bce915b0-5133-4514-8199-896e8ffffd82
md"## DMS data"

# ╔═╡ f8acffcf-d97b-4439-92f3-5dea5d13378e
dms_df = SamApp2024.load_dms_data_sequences_table_20250303_with_aligned_sequences();

# ╔═╡ 1300bb9d-c260-47f7-806b-1ef6f42242bd
dms_data = SamApp2024.load_dms_data_20250303();

# ╔═╡ 0b8c5578-f139-4757-950c-3b81ea60ee3d
dms_data_primers = dms_df.primer[[only(findall(dms_df.name .== name)) for name = dms_data.aptamer_names]];

# ╔═╡ ab307bc8-00eb-4fc2-8d6b-1a840a9de24c
dms_num_sites, dms_num_seqs, dms_num_conds = size(dms_data.shape_reactivities);

# ╔═╡ 2dd203e7-b1db-4539-b304-ef045d3f7a47
all_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=1:dms_num_sites for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ aee76cf7-7464-427b-b892-ffbf69da555a
bps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=bps for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ e744dc3f-d1e5-4e1d-91f3-3dc5471ab4fc
nps_dms_reactivities = [dms_data.shape_reactivities[i,n,1] for i=nps for n=1:dms_num_seqs if (!ismissing(dms_data.aligned_sequence[n])) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')];

# ╔═╡ 137e9b1c-258a-415a-88af-82d7cc851782
dms_stats = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = dms_data,
    paired_reactivities = bps_dms_reactivities,
    unpaired_reactivities = nps_dms_reactivities,
    all_reactivities = all_dms_reactivities,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ e4ea0a52-3118-4fe9-8a29-73575f6bd608
dms_data.conditions

# ╔═╡ e3855c4c-4c7f-4918-a45e-bfcd48237382
x_mg_dms = [ismissing(dms_data.aligned_sequence[n]) ? NaN : nansum([dms_stats.shape_log_odds[i,n,2] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:400];

# ╔═╡ 5c1f1266-b359-4f08-808e-8f43cf93ce6b
x_sam_dms = [ismissing(dms_data.aligned_sequence[n]) ? NaN : nansum([dms_stats.shape_log_odds[i,n,1] for i = _sites if dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for n=1:400];

# ╔═╡ Cell order:
# ╠═1ad1f988-4794-11f0-0ccb-e3e31a210177
# ╠═0fe3899a-1da4-4c5e-bec0-ba1fdeefbb24
# ╠═e15eacbd-c78a-406c-9fbc-046d81ac0d03
# ╠═14c7987f-8487-4e9c-9998-36714829fcad
# ╠═0e41f586-a85b-496f-9503-f456286a9d9c
# ╠═32433288-8243-4143-b301-69a30f60f3d2
# ╠═845b6fe0-93c1-400e-97d3-ef78b27ed175
# ╠═2386a288-a090-46ae-9b28-a1c90dbf5c25
# ╠═9f34f37d-74b3-48c0-96eb-4f0d2c863704
# ╠═35675b3d-e73e-42e2-ac2d-7202483872cc
# ╠═2b3a2e14-9f10-4c55-b301-447ce730b15b
# ╠═a056c80e-f046-4fed-af48-49abec41474b
# ╠═f9c24c3e-5697-4ef0-ab9d-932b2d115406
# ╠═da043e8d-e05b-4861-baef-9ef666431cf7
# ╠═c15c31a3-3723-4daf-872e-b2f9806530ce
# ╠═bce915b0-5133-4514-8199-896e8ffffd82
# ╠═f8acffcf-d97b-4439-92f3-5dea5d13378e
# ╠═1300bb9d-c260-47f7-806b-1ef6f42242bd
# ╠═0b8c5578-f139-4757-950c-3b81ea60ee3d
# ╠═ab307bc8-00eb-4fc2-8d6b-1a840a9de24c
# ╠═2dd203e7-b1db-4539-b304-ef045d3f7a47
# ╠═aee76cf7-7464-427b-b892-ffbf69da555a
# ╠═e744dc3f-d1e5-4e1d-91f3-3dc5471ab4fc
# ╠═137e9b1c-258a-415a-88af-82d7cc851782
# ╠═e4ea0a52-3118-4fe9-8a29-73575f6bd608
# ╠═e3855c4c-4c7f-4918-a45e-bfcd48237382
# ╠═5c1f1266-b359-4f08-808e-8f43cf93ce6b

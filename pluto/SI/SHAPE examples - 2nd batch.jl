### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 277a57b4-a89f-4df9-9dd2-cd90886ed94b
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 25fcd520-29df-4dfa-aacf-ace666481dff
using BioSequences: LongRNA

# ╔═╡ f7e3be77-4136-4e49-8e4a-d2bc8ae61b7b
using DataFrames: DataFrame

# ╔═╡ 7fe7b225-7213-4062-8921-a5d09b885d83
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ ed87e4b5-2443-459d-86cb-e6ce38722181
using LinearAlgebra: Diagonal, eigen

# ╔═╡ eb89bfdf-9882-482f-bdb2-cc9175ff1913
using Makie: @L_str

# ╔═╡ 0b024795-9179-49d0-a881-356109db2524
using NaNStatistics: nansum

# ╔═╡ acf92e99-bca4-444b-bee7-b863bfd4e479
using Random: bitrand

# ╔═╡ 19e7563f-9e5c-4578-b89d-f61cf449dcb7
using Statistics: cor, mean

# ╔═╡ 13345063-8ebe-4c7e-954e-e609cff1bf92
using StatsBase: countmap

# ╔═╡ c17512b3-d42b-4cb5-a74d-1c6543a4ccc4
md"# Imports"

# ╔═╡ a9e64dab-d6ef-4768-a9e9-759c8aec4e5c
import Makie

# ╔═╡ fadf7faa-39f8-47fb-bd63-96edb35d0cc4
import CairoMakie

# ╔═╡ 7d94efd5-2571-4d74-9234-1685b3002de2
import CSV, HDF5

# ╔═╡ 5c6c89b5-f3eb-47b0-acd5-b7ac2a546d6d
import FASTX, Infernal

# ╔═╡ e9519585-697f-4935-a9c9-ea861a7850ca
import SamApp2024

# ╔═╡ a7665277-f88b-4438-88ae-3eabda660a49
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ cf036136-5f5c-4961-b50a-2f65cb1e7c21
import Rfam

# ╔═╡ c88b3fe0-20ad-4668-83ce-325dc6cc9e63
import PlutoUI

# ╔═╡ d9d0715c-9140-4125-9312-ccf1b64fa10c
import Unitful

# ╔═╡ c04f8c58-57a4-41d2-a47e-1e84d8b0ee10
import ViennaRNA

# ╔═╡ 656da8f3-b7ba-4f77-90c8-fd76baa192b5
import StatsBase

# ╔═╡ 92a651e3-976f-4989-9418-2306c89ea1eb
import KernelDensity

# ╔═╡ 9e79cd3d-9a7a-4484-9bda-dc52d1a1934a
PlutoUI.TableOfContents()

# ╔═╡ 020914fa-6247-406d-810f-c823de40cb01
md"# Load data"

# ╔═╡ b08754f8-3cd9-4dd9-b86e-8bc4b263ccaa
shape_data = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ f24fa72b-f44b-4c2a-9d6c-99dbe78cf83b
RF00162_hits = SamApp2024.rfam_RF00162_hits();

# ╔═╡ e41bdff5-4805-4fbf-bd2a-a960c06ec469
rbm_energies = free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(shape_data.aligned_sequences));

# ╔═╡ 66deb70b-af97-4c31-b20e-0e9c2cdbb99c
aptamer_natural_distances = SamApp2024.hamming(SamApp2024.onehot(shape_data.aligned_sequences), onehot(RF00162_hits));

# ╔═╡ cb48035d-af9b-4693-ada6-2cc9fb06923b
conds_sam = [1,2];

# ╔═╡ ddc6001d-367b-4e4c-84b8-09ead2a32bc4
conds_mg = [4];

# ╔═╡ eaea7032-c5d9-4852-aca0-12c9620296ec
conds_30C = [6];

# ╔═╡ 5a61fb4e-6047-4351-8c90-9abf7fa2d052
the_conds = vcat(conds_sam, conds_mg, conds_30C)

# ╔═╡ 9a8280ed-2327-4a6c-9b18-b491adb0f0af
shape_data.conditions[the_conds]

# ╔═╡ cb32e1e6-2be0-4370-94b6-abb65bd061f4
bps, nps, pks = SamApp2024.RF00162_sites_paired()

# ╔═╡ 942d6eb1-f07c-4675-a2fa-727b13db9636
ss_sites = SamApp2024.RF00162_sites_annotated_secondary_structure()

# ╔═╡ b15062d4-2b58-406e-8c2a-8a9fd1cc50c3
bps_reactivities = shape_data.shape_reactivities[bps, :, conds_sam];

# ╔═╡ 3a7cec1e-2444-4914-b22d-9a24f14aae04
nps_reactivities = shape_data.shape_reactivities[nps, :, conds_sam];

# ╔═╡ be0dd03e-0aae-4c94-b532-976d7e2e0385
all_reactivities = shape_data.shape_reactivities[:, :, conds_sam];

# ╔═╡ 619b8305-cd3d-47c1-aba9-bd4f363ed204
shape_stats = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data,
    paired_reactivities = bps_reactivities,
    unpaired_reactivities = nps_reactivities,
    all_reactivities = all_reactivities,
    only_hq_profile = true, p_thresh = 1e-2, nsamples=5000
);

# ╔═╡ 4c92b8bb-2885-4f63-b485-ea15e1bce1ba
_sites = SamApp2024.hallmark_sites_20230507

# ╔═╡ Cell order:
# ╠═c17512b3-d42b-4cb5-a74d-1c6543a4ccc4
# ╠═277a57b4-a89f-4df9-9dd2-cd90886ed94b
# ╠═a9e64dab-d6ef-4768-a9e9-759c8aec4e5c
# ╠═fadf7faa-39f8-47fb-bd63-96edb35d0cc4
# ╠═7d94efd5-2571-4d74-9234-1685b3002de2
# ╠═5c6c89b5-f3eb-47b0-acd5-b7ac2a546d6d
# ╠═e9519585-697f-4935-a9c9-ea861a7850ca
# ╠═a7665277-f88b-4438-88ae-3eabda660a49
# ╠═cf036136-5f5c-4961-b50a-2f65cb1e7c21
# ╠═c88b3fe0-20ad-4668-83ce-325dc6cc9e63
# ╠═d9d0715c-9140-4125-9312-ccf1b64fa10c
# ╠═c04f8c58-57a4-41d2-a47e-1e84d8b0ee10
# ╠═656da8f3-b7ba-4f77-90c8-fd76baa192b5
# ╠═92a651e3-976f-4989-9418-2306c89ea1eb
# ╠═25fcd520-29df-4dfa-aacf-ace666481dff
# ╠═f7e3be77-4136-4e49-8e4a-d2bc8ae61b7b
# ╠═7fe7b225-7213-4062-8921-a5d09b885d83
# ╠═ed87e4b5-2443-459d-86cb-e6ce38722181
# ╠═eb89bfdf-9882-482f-bdb2-cc9175ff1913
# ╠═0b024795-9179-49d0-a881-356109db2524
# ╠═acf92e99-bca4-444b-bee7-b863bfd4e479
# ╠═19e7563f-9e5c-4578-b89d-f61cf449dcb7
# ╠═13345063-8ebe-4c7e-954e-e609cff1bf92
# ╠═9e79cd3d-9a7a-4484-9bda-dc52d1a1934a
# ╠═020914fa-6247-406d-810f-c823de40cb01
# ╠═b08754f8-3cd9-4dd9-b86e-8bc4b263ccaa
# ╠═f24fa72b-f44b-4c2a-9d6c-99dbe78cf83b
# ╠═e41bdff5-4805-4fbf-bd2a-a960c06ec469
# ╠═66deb70b-af97-4c31-b20e-0e9c2cdbb99c
# ╠═cb48035d-af9b-4693-ada6-2cc9fb06923b
# ╠═ddc6001d-367b-4e4c-84b8-09ead2a32bc4
# ╠═eaea7032-c5d9-4852-aca0-12c9620296ec
# ╠═5a61fb4e-6047-4351-8c90-9abf7fa2d052
# ╠═9a8280ed-2327-4a6c-9b18-b491adb0f0af
# ╠═cb32e1e6-2be0-4370-94b6-abb65bd061f4
# ╠═942d6eb1-f07c-4675-a2fa-727b13db9636
# ╠═b15062d4-2b58-406e-8c2a-8a9fd1cc50c3
# ╠═3a7cec1e-2444-4914-b22d-9a24f14aae04
# ╠═be0dd03e-0aae-4c94-b532-976d7e2e0385
# ╠═619b8305-cd3d-47c1-aba9-bd4f363ed204
# ╠═4c92b8bb-2885-4f63-b485-ea15e1bce1ba

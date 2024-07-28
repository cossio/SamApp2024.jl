### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 0b2523a6-5530-4cf0-a22a-560a8f9549da
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 24f9e3a2-ae85-4d3d-a4cd-4281c07492c7
using BioSequences: LongRNA

# ╔═╡ 7ec13690-0410-4828-be59-49d698fc3448
using DataFrames: DataFrame

# ╔═╡ b70c3f5a-7f68-4ce9-891d-3eea9a2ba357
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ 9f724968-79b3-4f64-9471-33fcdd50b4c6
using LinearAlgebra: Diagonal, eigen

# ╔═╡ 55008aca-d6f9-4c50-9d84-84563463b235
using Makie: @L_str

# ╔═╡ b03079c8-a78a-4cca-9aa0-e77e19654f59
using NaNStatistics: nansum

# ╔═╡ 2e704b97-0baf-4cf2-b8e0-56e1a931c324
using Random: bitrand

# ╔═╡ c410f9a3-ff34-49ed-9a17-2d428f721d03
using Statistics: cor, mean

# ╔═╡ 7d343199-eb06-4377-903e-6bc45e60dfa3
using StatsBase: countmap

# ╔═╡ ac453f51-0ce2-4d9e-bc5b-1e4d6342170a
md"# Imports"

# ╔═╡ 60d9b880-9b75-431c-bb12-69d9db711b11
import Makie, CairoMakie

# ╔═╡ a7cf46e1-198d-4381-b5b3-cf781754f0e8
import CSV, HDF5

# ╔═╡ d7a11d40-768c-47c5-80b7-99fecd60f7fd
import FASTX, Infernal

# ╔═╡ af15dc75-cd8e-470f-a765-46173388a41a
import SamApp2024

# ╔═╡ bc474b44-5591-4893-8c69-fe66d55d2d9d
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 83974feb-e04a-4fe3-a4ce-e480578ba4ce
import Rfam

# ╔═╡ 0eede979-3e59-432d-9d02-7dd31d862c46
import PlutoUI

# ╔═╡ a6b2a41f-c516-4d47-b2e2-8454095d96e8
import StatsBase, KernelDensity

# ╔═╡ 88a9dd2e-5240-464f-8092-b1d2a16219da
PlutoUI.TableOfContents()

# ╔═╡ c460a3e0-9af4-40a8-8d14-b322d4047803
md"# Functions"

# ╔═╡ 40f244ae-92af-4bda-92e2-0ce55836578b
function hamming_nodiag(seqs)
    d = SamApp2024.hamming(seqs)
    return [d[i,j] for i in axes(d,1) for j in axes(d,2) if i ≠ j]
end

# ╔═╡ 4028ac22-a508-4cfb-8493-a9fc4e1421d3
md"# Load data"

# ╔═╡ ba074af4-5acd-48db-baa9-7e2fc7307e9d
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162")

# ╔═╡ 0a13b8b1-c704-43c8-94eb-1a1df31b88ff
RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")

# ╔═╡ f639c9db-374f-465e-ae28-2d25a0fde2e4
RF00162_seed_match_cols = findall(≠('.'), SamApp2024.stockholm_ss(RF00162_seed_stk.out));

# ╔═╡ 5f76b82f-7768-42f9-bb77-6cdbd5667ffb
RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM") # WARNING: this has inserts marked as '-'

# ╔═╡ b9f73363-7341-47f4-891d-957c2f224bb9
RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))

# ╔═╡ 51307a99-1ff4-45e4-8b25-1bddfa7f539f
RF00162_seed_seqs_noinserts = LongRNA{4}.([FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records]);

# ╔═╡ 59c76610-5bc1-4691-b91f-353429995bc2
# trimmed (no inserts) aligned fasta

# ╔═╡ f3423a49-e8d6-43f7-b0c6-acb800347040
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");
# these are already aligned and without inserts

# ╔═╡ c18fcd1c-efff-4b37-952d-4894870de951
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ 5e7c7128-93e8-48c6-adae-334690618757
rbm = SamApp2024.rbm2022();

# ╔═╡ c75b7167-6ed8-498a-9a39-748b37facda6
sampled_v = SamApp2024.rbm2022samples();

# ╔═╡ 410bb137-4e1d-4f63-a2e5-8e3c2d3206c5
md"# Compute distances"

# ╔═╡ fe3431c3-a6b6-4040-955d-9f255c51caab
RF00162_distances = hamming_nodiag(RF00162_hits_sequences);

# ╔═╡ Cell order:
# ╠═ac453f51-0ce2-4d9e-bc5b-1e4d6342170a
# ╠═0b2523a6-5530-4cf0-a22a-560a8f9549da
# ╠═60d9b880-9b75-431c-bb12-69d9db711b11
# ╠═a7cf46e1-198d-4381-b5b3-cf781754f0e8
# ╠═d7a11d40-768c-47c5-80b7-99fecd60f7fd
# ╠═af15dc75-cd8e-470f-a765-46173388a41a
# ╠═bc474b44-5591-4893-8c69-fe66d55d2d9d
# ╠═83974feb-e04a-4fe3-a4ce-e480578ba4ce
# ╠═0eede979-3e59-432d-9d02-7dd31d862c46
# ╠═a6b2a41f-c516-4d47-b2e2-8454095d96e8
# ╠═24f9e3a2-ae85-4d3d-a4cd-4281c07492c7
# ╠═7ec13690-0410-4828-be59-49d698fc3448
# ╠═b70c3f5a-7f68-4ce9-891d-3eea9a2ba357
# ╠═9f724968-79b3-4f64-9471-33fcdd50b4c6
# ╠═55008aca-d6f9-4c50-9d84-84563463b235
# ╠═b03079c8-a78a-4cca-9aa0-e77e19654f59
# ╠═2e704b97-0baf-4cf2-b8e0-56e1a931c324
# ╠═c410f9a3-ff34-49ed-9a17-2d428f721d03
# ╠═7d343199-eb06-4377-903e-6bc45e60dfa3
# ╠═88a9dd2e-5240-464f-8092-b1d2a16219da
# ╠═c460a3e0-9af4-40a8-8d14-b322d4047803
# ╠═40f244ae-92af-4bda-92e2-0ce55836578b
# ╠═4028ac22-a508-4cfb-8493-a9fc4e1421d3
# ╠═ba074af4-5acd-48db-baa9-7e2fc7307e9d
# ╠═0a13b8b1-c704-43c8-94eb-1a1df31b88ff
# ╠═f639c9db-374f-465e-ae28-2d25a0fde2e4
# ╠═5f76b82f-7768-42f9-bb77-6cdbd5667ffb
# ╠═b9f73363-7341-47f4-891d-957c2f224bb9
# ╠═51307a99-1ff4-45e4-8b25-1bddfa7f539f
# ╠═59c76610-5bc1-4691-b91f-353429995bc2
# ╠═f3423a49-e8d6-43f7-b0c6-acb800347040
# ╠═c18fcd1c-efff-4b37-952d-4894870de951
# ╠═5e7c7128-93e8-48c6-adae-334690618757
# ╠═c75b7167-6ed8-498a-9a39-748b37facda6
# ╠═410bb137-4e1d-4f63-a2e5-8e3c2d3206c5
# ╠═fe3431c3-a6b6-4040-955d-9f255c51caab

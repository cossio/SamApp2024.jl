### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 952a79aa-12fa-11ef-2326-d5fe00ff3dfe
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ adbed31f-5777-4b6f-acf9-f37a472ef5d8
using BioSequences: LongRNA

# ╔═╡ de469256-6471-4fa9-854b-7b7783295aef
using DataFrames: DataFrame

# ╔═╡ a9fb22f9-4c47-4430-a0aa-cc4075b86f1f
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ ebabc892-15f1-4b91-be5d-b85676d94fa7
using LinearAlgebra: Diagonal, eigen

# ╔═╡ 2a184e1f-bcba-48a3-a54b-793a25d9501e
using Makie: @L_str

# ╔═╡ 05b55f01-751c-4910-964f-aadb33f844c0
using NaNStatistics: nansum

# ╔═╡ a72c9391-c8f5-4bf5-94bd-3e9bbc326fa2
using Random: bitrand

# ╔═╡ bf22e1b2-f8bd-4990-938d-4db503660c3d
using Statistics: cor, mean

# ╔═╡ ac2f96ec-b31f-4f63-8fd3-4390d687755b
using StatsBase: countmap

# ╔═╡ 5d5eb95f-6801-47c3-a82e-afab35559c0e
md"# Imports"

# ╔═╡ f01ccb33-26f8-4fb9-94a2-88be47db3762
import Makie, CairoMakie

# ╔═╡ 897ffcd7-e59a-4121-975d-29bc9a27988b
import CSV, HDF5

# ╔═╡ 49b0a181-8c3b-41dc-b3e0-9e55efa87bd9
import FASTX, Infernal

# ╔═╡ 4a2a79bc-6fb1-45d3-ad4a-abfc889f2737
import SamApp2024

# ╔═╡ 5ca1d179-8dda-4e4b-9bab-b1cf4f6bafbe
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 8bed6401-4c1a-453d-ac16-ff4c81346936
import Rfam

# ╔═╡ 98c97d8b-9849-48dc-ae6c-c0924d14e560
import StatsBase, KernelDensity

# ╔═╡ aa8a002e-693a-402b-aaf1-576bf59e3a4e
md"# Plot"

# ╔═╡ 23712f17-3947-4ed3-85bd-00dceb854167
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162");

# ╔═╡ 660fa41b-c6c2-4b41-a698-82530c8dc4ef
RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")

# ╔═╡ 8399bdac-bd99-4a20-be2c-2055db32a5a6
RF00162_seed_match_cols = findall(≠('.'), SamApp2024.stockholm_ss(RF00162_seed_stk.out));

# ╔═╡ fd9539ff-e725-4b18-94fd-e6a75fc2ec73
RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM") # WARNING: this has inserts marked as '-'

# ╔═╡ 6538d3e8-54d8-484d-8baa-3d3cc44ed72c
RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))

# ╔═╡ 71b73684-6053-4efc-92b7-391267f26ad9
RF00162_seed_seqs_noinserts = LongRNA{4}.([FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records]);

# ╔═╡ e4309e31-4086-4b3f-baec-8798ab975c44
# trimmed (no inserts) aligned fasta

# ╔═╡ 3a5e2efa-a739-42c0-a76c-50c6c07487e9
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");
# these are already aligned and without inserts

# ╔═╡ 0d4f2776-4666-48dd-80df-990c479569df
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ 488c249d-9810-4697-a6ca-7ab326b53762
# emit sequences from Rfam CM model
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Rfam_cm.out; N=5000, aligned=true, outformat="AFA");

# ╔═╡ 0eeacfc3-2240-46e2-b851-9d574771d1d7
begin
	Rfam_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)));
	Rfam_cm_emitted_sequences = [filter(!=('.'), filter(!islowercase, seq)) for seq in Rfam_cm_emitted_sequences];
	Rfam_cm_emitted_sequences = LongRNA{4}.(Rfam_cm_emitted_sequences);
end

# ╔═╡ 896d95be-2a0c-4c4b-9c73-d6e2d1452306
begin
	# aligned hits, used to train a new noiseless CM model (in Stockholm format, without inserts!)
	RF00162_hits_stk = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true);
	# fit new CM model using full alignment (without inserts), and without entropic noise
	Refined_cm = Infernal.cmbuild(RF00162_hits_stk.out; enone=true);
	
	# emit sequences from Refined CM model
	Refined_cm_emitted_sequences_afa = Infernal.cmemit(Refined_cm.cmout; N=5000, aligned=true, outformat="AFA");
	Refined_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Refined_cm_emitted_sequences_afa.out)));
	# remove inserts
	Refined_cm_emitted_sequences = [filter(!=('.'), filter(!islowercase, seq)) for seq in Refined_cm_emitted_sequences];
	@assert only(unique(length.(Refined_cm_emitted_sequences))) == 108
	Refined_cm_emitted_sequences = LongRNA{4}.(Refined_cm_emitted_sequences);
end

# ╔═╡ 0bd3ffae-c009-4740-aa57-4dd5126eb848
# use saved RBM samples
sampled_v = SamApp2024.rbm2022samples(); # faster

# ╔═╡ Cell order:
# ╠═5d5eb95f-6801-47c3-a82e-afab35559c0e
# ╠═952a79aa-12fa-11ef-2326-d5fe00ff3dfe
# ╠═f01ccb33-26f8-4fb9-94a2-88be47db3762
# ╠═897ffcd7-e59a-4121-975d-29bc9a27988b
# ╠═49b0a181-8c3b-41dc-b3e0-9e55efa87bd9
# ╠═4a2a79bc-6fb1-45d3-ad4a-abfc889f2737
# ╠═5ca1d179-8dda-4e4b-9bab-b1cf4f6bafbe
# ╠═8bed6401-4c1a-453d-ac16-ff4c81346936
# ╠═98c97d8b-9849-48dc-ae6c-c0924d14e560
# ╠═adbed31f-5777-4b6f-acf9-f37a472ef5d8
# ╠═de469256-6471-4fa9-854b-7b7783295aef
# ╠═a9fb22f9-4c47-4430-a0aa-cc4075b86f1f
# ╠═ebabc892-15f1-4b91-be5d-b85676d94fa7
# ╠═2a184e1f-bcba-48a3-a54b-793a25d9501e
# ╠═05b55f01-751c-4910-964f-aadb33f844c0
# ╠═a72c9391-c8f5-4bf5-94bd-3e9bbc326fa2
# ╠═bf22e1b2-f8bd-4990-938d-4db503660c3d
# ╠═ac2f96ec-b31f-4f63-8fd3-4390d687755b
# ╠═aa8a002e-693a-402b-aaf1-576bf59e3a4e
# ╠═23712f17-3947-4ed3-85bd-00dceb854167
# ╠═660fa41b-c6c2-4b41-a698-82530c8dc4ef
# ╠═8399bdac-bd99-4a20-be2c-2055db32a5a6
# ╠═fd9539ff-e725-4b18-94fd-e6a75fc2ec73
# ╠═6538d3e8-54d8-484d-8baa-3d3cc44ed72c
# ╠═71b73684-6053-4efc-92b7-391267f26ad9
# ╠═e4309e31-4086-4b3f-baec-8798ab975c44
# ╠═3a5e2efa-a739-42c0-a76c-50c6c07487e9
# ╠═0d4f2776-4666-48dd-80df-990c479569df
# ╠═488c249d-9810-4697-a6ca-7ab326b53762
# ╠═0eeacfc3-2240-46e2-b851-9d574771d1d7
# ╠═896d95be-2a0c-4c4b-9c73-d6e2d1452306
# ╠═0bd3ffae-c009-4740-aa57-4dd5126eb848

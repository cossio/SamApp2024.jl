### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ d1c8de8f-c7f8-48d2-b25b-3037109a0995
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ eadccc83-ed19-427e-87b3-565714189ce8
using BioSequences: LongRNA

# ╔═╡ 774dd693-e5fe-41af-b317-8b37344452eb
using BioSequences: LongDNA

# ╔═╡ ff8f1c6c-576a-4f98-862f-1db75a308271
using BioSequences: translate

# ╔═╡ 10965851-2849-407d-a5a1-c0258f1b9750
using DataFrames: DataFrame

# ╔═╡ 7c6372e4-bef4-4a7b-80ea-87eadcb50878
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ d90e70d9-13f5-4822-a44b-94ff806d271f
using LinearAlgebra: Diagonal, eigen

# ╔═╡ 6bf569a2-edb3-45f0-8feb-ca441d54aedf
using Makie: @L_str

# ╔═╡ bb84dd4d-89d3-4a0a-8ac9-a65f0b1de8dd
using NaNStatistics: nansum

# ╔═╡ b67bca82-8fb2-421b-8eb5-a644622e6700
using Random: bitrand

# ╔═╡ 23e57160-8135-46c5-be49-5da2ba195909
using Statistics: cor, mean

# ╔═╡ 8d60afdc-642f-4cc7-aeb8-4180c97ed792
using StatsBase: countmap

# ╔═╡ eb11da0b-85e5-49fd-9aaf-2386cb654d4d
md"# Imports"

# ╔═╡ 2b14fd98-ab49-4592-be7d-031ec5240b1c
import Makie, CairoMakie

# ╔═╡ b09c0e3b-ed25-40f8-8183-f4bc5b3ebb3d
import CSV, HDF5

# ╔═╡ 089bdb81-cf76-467b-b183-8a080d8b4053
import FASTX, Infernal

# ╔═╡ 1770e39a-1b98-4e1f-8512-3f595860b809
import SamApp2024

# ╔═╡ b636df11-869c-47fd-ae41-c077185513f3
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ cfee5d6b-f94c-4776-afb3-57111f894865
import Rfam

# ╔═╡ 65ddb1dd-dc55-4115-bb28-ad2ae04f7548
import StatsBase, KernelDensity

# ╔═╡ 89cb99bd-6df1-43fc-bc7f-70f1b33cdf32
md"# Plots"

# ╔═╡ 2c06e1fc-da01-4461-8339-f4ffb20f4b58
Rfam_id = "RF00234"

# ╔═╡ 1e3b2229-e00a-4dc9-a639-95d97521e81a
Rfam_cm = Infernal.cmfetch(Rfam.cm(), Rfam_id).out

# ╔═╡ e623120c-55f0-4542-a5f4-e5eb29b61bec
rfamgen_seqs_fasta = "/DATA/cossio/SAM/2024/rfamgen/outputs/$Rfam_id/sampled_seq.fa"

# ╔═╡ dbdd126f-fa4e-4c8e-a431-82fb2be2bcad
rfamgen_seqs_raw = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(rfamgen_seqs_fasta))))

# ╔═╡ a4d55211-0a6c-4b42-a3f6-c1d100f4e6a6
rfamgen_seqs = SamApp2024.infernal_align_fasta_to_cm(rfamgen_seqs_fasta, Rfam_cm)

# ╔═╡ 8b352286-539f-4350-a0b2-a4a7e1bc575f
hits_afa = Infernal.cmalign(Rfam_cm, Rfam.fasta_file(Rfam_id); matchonly=true, outformat="AFA")

# ╔═╡ 0eff79c7-4f68-4458-96c5-42d3a3b1e3a0
hits_raw = LongRNA{4}.([replace(s, 'T' => 'U') for s = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam.fasta_file(Rfam_id))))])

# ╔═╡ 50ea4846-a140-46e6-9df5-1c74b014a0ef
hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(hits_afa.out))))

# ╔═╡ 92c8ef49-8302-487d-9f1b-5ff4d9d41ca2
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="Bit score", ylabel="frequency")
	Makie.hist!(ax, SamApp2024.infernal_score_sequences(Rfam_cm, rfamgen_seqs_raw).bit_sc; normalization=:pdf, label="RfamGen")
	Makie.hist!(ax, SamApp2024.infernal_score_sequences(Rfam_cm, hits_raw).bit_sc; normalization=:pdf, label="MSA")
	Makie.axislegend(ax)
	fig
end

# ╔═╡ b4407a24-ccab-49b8-a2dc-a48511e83733
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="Bit score", ylabel="frequency")
	Makie.hist!(ax, map(length, rfamgen_seqs_raw); normalization=:pdf, label="RfamGen")
	Makie.hist!(ax, map(length, hits_raw); normalization=:pdf, label="MSA")
	Makie.axislegend(ax)
	fig
end

# ╔═╡ Cell order:
# ╠═eb11da0b-85e5-49fd-9aaf-2386cb654d4d
# ╠═d1c8de8f-c7f8-48d2-b25b-3037109a0995
# ╠═2b14fd98-ab49-4592-be7d-031ec5240b1c
# ╠═b09c0e3b-ed25-40f8-8183-f4bc5b3ebb3d
# ╠═089bdb81-cf76-467b-b183-8a080d8b4053
# ╠═1770e39a-1b98-4e1f-8512-3f595860b809
# ╠═b636df11-869c-47fd-ae41-c077185513f3
# ╠═cfee5d6b-f94c-4776-afb3-57111f894865
# ╠═65ddb1dd-dc55-4115-bb28-ad2ae04f7548
# ╠═eadccc83-ed19-427e-87b3-565714189ce8
# ╠═774dd693-e5fe-41af-b317-8b37344452eb
# ╠═ff8f1c6c-576a-4f98-862f-1db75a308271
# ╠═10965851-2849-407d-a5a1-c0258f1b9750
# ╠═7c6372e4-bef4-4a7b-80ea-87eadcb50878
# ╠═d90e70d9-13f5-4822-a44b-94ff806d271f
# ╠═6bf569a2-edb3-45f0-8feb-ca441d54aedf
# ╠═bb84dd4d-89d3-4a0a-8ac9-a65f0b1de8dd
# ╠═b67bca82-8fb2-421b-8eb5-a644622e6700
# ╠═23e57160-8135-46c5-be49-5da2ba195909
# ╠═8d60afdc-642f-4cc7-aeb8-4180c97ed792
# ╠═89cb99bd-6df1-43fc-bc7f-70f1b33cdf32
# ╠═2c06e1fc-da01-4461-8339-f4ffb20f4b58
# ╠═1e3b2229-e00a-4dc9-a639-95d97521e81a
# ╠═e623120c-55f0-4542-a5f4-e5eb29b61bec
# ╠═dbdd126f-fa4e-4c8e-a431-82fb2be2bcad
# ╠═a4d55211-0a6c-4b42-a3f6-c1d100f4e6a6
# ╠═8b352286-539f-4350-a0b2-a4a7e1bc575f
# ╠═0eff79c7-4f68-4458-96c5-42d3a3b1e3a0
# ╠═50ea4846-a140-46e6-9df5-1c74b014a0ef
# ╠═92c8ef49-8302-487d-9f1b-5ff4d9d41ca2
# ╠═b4407a24-ccab-49b8-a2dc-a48511e83733

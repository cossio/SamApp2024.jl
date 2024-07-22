### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ c9ddb318-70db-4592-afde-80772aaa2244
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 794b3eb0-55b8-4f5c-ba03-675e23d2c471
using BioSequences: LongRNA

# ╔═╡ 62704882-f5f7-44d3-9c2b-a5a72212319f
using DataFrames: DataFrame

# ╔═╡ 96f37f86-ee7d-47b6-be86-8984823354a4
using Distributions: Gamma, logpdf, pdf, Poisson

# ╔═╡ 58e89990-976b-4308-86d8-cadf3986a98c
using LinearAlgebra: Diagonal, eigen

# ╔═╡ ef0e0da8-d64a-4496-94a5-74d683b8cf4d
using Makie: @L_str

# ╔═╡ f2ccacb2-dab7-4054-b25d-7075f8a20580
using NaNStatistics: nansum

# ╔═╡ 9058d796-56ba-4fed-bc10-2684bfc07add
using Random: bitrand

# ╔═╡ c7e6ceff-23f5-4ad5-be6f-a968f5bdf030
using Statistics: cor, mean

# ╔═╡ 2e28c501-fdca-46ee-b354-309d632c24bf
using StatsBase: countmap

# ╔═╡ e98903a9-5faf-4cc0-ab42-a7a7bd5ccfc9
md"# Imports" 

# ╔═╡ 7327676a-736e-498a-b80e-5d74dc163a35
import Makie, CairoMakie

# ╔═╡ f483c527-e02a-477e-a1a9-62aa75613ca0
import CSV, HDF5

# ╔═╡ f847fe47-b2b5-46c6-bf43-41ee3887c56b
import FASTX, Infernal

# ╔═╡ ccc08755-a3ca-43c3-bf3d-ebf6a32afa2c
import SamApp2024

# ╔═╡ f05273b2-c9b4-4b15-96e0-bdc3f45dfd63
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 5492b03c-c32e-4471-9748-ec46b7e1dbaa
import Rfam

# ╔═╡ 8616a37c-9f23-41c8-85b5-654b6ec7689b
import StatsBase, KernelDensity

# ╔═╡ 061e8cb0-6c3a-4e1a-93d3-eeb400fca037
md"# Load data"

# ╔═╡ b638a0dc-a03a-482e-8479-08da23b40494
Rfam_id = "RF00162"

# ╔═╡ 94c85a6b-8061-43cf-a92b-4e87cac1e147
Rfam_cm = Infernal.cmfetch(Rfam.cm(), Rfam_id).out

# ╔═╡ 9d2c22af-e34c-433e-8915-f7dd2f2214b6
rfamgen_seqs_fasta = "/DATA/cossio/SAM/2024/rfamgen/outputs/$Rfam_id/sampled_seq.fa"

# ╔═╡ ee4923e5-53c6-4410-a5c0-47e6c57ee4b6
rfamgen_seqs_aln = SamApp2024.infernal_align_fasta_to_cm(rfamgen_seqs_fasta, Rfam_cm)

# ╔═╡ 504ae038-7b2f-4aa3-b0dc-97372b506daa
rfamgen_seqs_raw = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(rfamgen_seqs_fasta))))

# ╔═╡ 1d590af4-a11f-4cc7-890e-2c54c5bebb60
rfamgen_seqs_fasta_sampl = "/DATA/cossio/SAM/2024/rfamgen/outputs/$Rfam_id/sampled_seq_sampl.fa"

# ╔═╡ c036fba4-baff-4e8c-973b-fb9f059251ae
rfamgen_seqs_raw_sampl = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(rfamgen_seqs_fasta_sampl))))

# ╔═╡ fc2d1581-84d8-4812-b596-250b781ff6c2
rfamgen_seqs_aln_sampl = SamApp2024.infernal_align_fasta_to_cm(rfamgen_seqs_fasta_sampl, Rfam_cm)

# ╔═╡ ff8c1333-50f5-40c0-bd4f-1cb2f7653699
hits_afa = Infernal.cmalign(Rfam_cm, Rfam.fasta_file(Rfam_id); matchonly=true, outformat="AFA")

# ╔═╡ f2421b0c-df06-4f73-9993-1ce54f96c3cf
hits_raw = LongRNA{4}.([replace(s, 'T' => 'U') for s = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam.fasta_file(Rfam_id))))])

# ╔═╡ 4ac19098-8342-4bef-844f-969462eb73e6
hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(hits_afa.out))))

# ╔═╡ 1899a3c0-4e1b-44b3-a7b5-efad4d6b8fb3
cm_score_rfamgen = SamApp2024.infernal_score_sequences(Rfam_cm, rfamgen_seqs_raw).bit_sc

# ╔═╡ 3a0e3cd7-d589-4e04-8de3-99e058c42c46
cm_score_rfamgen_sampl = SamApp2024.infernal_score_sequences(Rfam_cm, rfamgen_seqs_raw_sampl).bit_sc

# ╔═╡ b236dc9c-ad03-4392-9027-ca76513bc0d0
cm_score_hits = SamApp2024.infernal_score_sequences(Rfam_cm, hits_raw).bit_sc

# ╔═╡ 5166ec65-2c64-4f21-84e3-043280318006
rbm = SamApp2024.rbm2022()

# ╔═╡ 53f710ec-9193-43c7-ad0c-14b6af9ac892
md"# Intermediate plots"

# ╔═╡ 6610e810-6d94-4a8e-b3e8-7308ee6100de
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="RBM score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(hits_sequences)); normalization=:pdf, label="MSA", bins=200:2:370, color=(:green, 0.5))
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln_sampl)); normalization=:pdf, label="RfamGen (sampled)", bins=200:2:370, color=:brown, linewidth=2)
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln)); normalization=:pdf, label="RfamGen (argmax)", bins=200:2:370, color=:brown, linewidth=2, linestyle=:dash)
	Makie.axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 048284c6-3a38-4196-8da7-ba55dea5bee0
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="Infernal bit score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, cm_score_hits; normalization=:pdf, label="MSA", bins=30:120, color=(:green, 0.5))
	Makie.stephist!(ax, cm_score_rfamgen_sampl; normalization=:pdf, label="RfamGen (sampled)", bins=30:120, color=:brown, linewidth=3)
	Makie.stephist!(ax, cm_score_rfamgen; normalization=:pdf, label="RfamGen (argmax)", bins=30:120, color=:brown, linewidth=2, linestyle=:dash)
	Makie.axislegend(ax; position=:lt)
	fig
end

# ╔═╡ 3a2cdb94-1c0a-4798-a14d-80beb09e7716
md"# Final figure"

# ╔═╡ 2b90f592-5a28-485e-84bf-232d900f1abc
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="RBM score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(hits_sequences)); normalization=:pdf, label="MSA", bins=200:2:370, color=(:green, 0.5))
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln_sampl)); normalization=:pdf, label="RfamGen (sampled)", bins=200:2:370, color=:brown, linewidth=2)
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln)); normalization=:pdf, label="RfamGen (argmax)", bins=200:2:370, color=:brown, linewidth=2, linestyle=:dash)
	Makie.axislegend(ax, position=:lt)

	ax = Makie.Axis(fig[1,2]; width=400, height=400, xlabel="Infernal bit score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, cm_score_hits; normalization=:pdf, label="MSA", bins=30:120, color=(:green, 0.5))
	Makie.stephist!(ax, cm_score_rfamgen_sampl; normalization=:pdf, label="RfamGen (sampled)", bins=30:120, color=:brown, linewidth=2)
	Makie.stephist!(ax, cm_score_rfamgen; normalization=:pdf, label="RfamGen (argmax)", bins=30:120, color=:brown, linewidth=2, linestyle=:dash)
	#Makie.axislegend(ax; position=:lt)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═e98903a9-5faf-4cc0-ab42-a7a7bd5ccfc9
# ╠═c9ddb318-70db-4592-afde-80772aaa2244
# ╠═7327676a-736e-498a-b80e-5d74dc163a35
# ╠═f483c527-e02a-477e-a1a9-62aa75613ca0
# ╠═f847fe47-b2b5-46c6-bf43-41ee3887c56b
# ╠═ccc08755-a3ca-43c3-bf3d-ebf6a32afa2c
# ╠═f05273b2-c9b4-4b15-96e0-bdc3f45dfd63
# ╠═5492b03c-c32e-4471-9748-ec46b7e1dbaa
# ╠═8616a37c-9f23-41c8-85b5-654b6ec7689b
# ╠═794b3eb0-55b8-4f5c-ba03-675e23d2c471
# ╠═62704882-f5f7-44d3-9c2b-a5a72212319f
# ╠═96f37f86-ee7d-47b6-be86-8984823354a4
# ╠═58e89990-976b-4308-86d8-cadf3986a98c
# ╠═ef0e0da8-d64a-4496-94a5-74d683b8cf4d
# ╠═f2ccacb2-dab7-4054-b25d-7075f8a20580
# ╠═9058d796-56ba-4fed-bc10-2684bfc07add
# ╠═c7e6ceff-23f5-4ad5-be6f-a968f5bdf030
# ╠═2e28c501-fdca-46ee-b354-309d632c24bf
# ╠═061e8cb0-6c3a-4e1a-93d3-eeb400fca037
# ╠═b638a0dc-a03a-482e-8479-08da23b40494
# ╠═94c85a6b-8061-43cf-a92b-4e87cac1e147
# ╠═9d2c22af-e34c-433e-8915-f7dd2f2214b6
# ╠═ee4923e5-53c6-4410-a5c0-47e6c57ee4b6
# ╠═504ae038-7b2f-4aa3-b0dc-97372b506daa
# ╠═1d590af4-a11f-4cc7-890e-2c54c5bebb60
# ╠═c036fba4-baff-4e8c-973b-fb9f059251ae
# ╠═fc2d1581-84d8-4812-b596-250b781ff6c2
# ╠═ff8c1333-50f5-40c0-bd4f-1cb2f7653699
# ╠═f2421b0c-df06-4f73-9993-1ce54f96c3cf
# ╠═4ac19098-8342-4bef-844f-969462eb73e6
# ╠═1899a3c0-4e1b-44b3-a7b5-efad4d6b8fb3
# ╠═3a0e3cd7-d589-4e04-8de3-99e058c42c46
# ╠═b236dc9c-ad03-4392-9027-ca76513bc0d0
# ╠═5166ec65-2c64-4f21-84e3-043280318006
# ╠═53f710ec-9193-43c7-ad0c-14b6af9ac892
# ╠═6610e810-6d94-4a8e-b3e8-7308ee6100de
# ╠═048284c6-3a38-4196-8da7-ba55dea5bee0
# ╠═3a2cdb94-1c0a-4798-a14d-80beb09e7716
# ╠═2b90f592-5a28-485e-84bf-232d900f1abc

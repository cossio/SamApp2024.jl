### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ c9ddb318-70db-4592-afde-80772aaa2244
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 794b3eb0-55b8-4f5c-ba03-675e23d2c471
using BioSequences: LongRNA

# ╔═╡ ef0e0da8-d64a-4496-94a5-74d683b8cf4d
using Makie: @L_str

# ╔═╡ f2ccacb2-dab7-4054-b25d-7075f8a20580
using NaNStatistics: nansum

# ╔═╡ c7e6ceff-23f5-4ad5-be6f-a968f5bdf030
using Statistics: cor

# ╔═╡ cb1b7457-21be-4c1e-9363-6258f4b9377a
using Statistics: mean

# ╔═╡ 2e28c501-fdca-46ee-b354-309d632c24bf
using DelimitedFiles: readdlm

# ╔═╡ e98903a9-5faf-4cc0-ab42-a7a7bd5ccfc9
md"# Imports" 

# ╔═╡ 7327676a-736e-498a-b80e-5d74dc163a35
import Makie

# ╔═╡ 07d606f2-0dc0-470b-83ed-8d0eb1a1a2be
import CairoMakie

# ╔═╡ f847fe47-b2b5-46c6-bf43-41ee3887c56b
import FASTX

# ╔═╡ 8962064d-dfda-4a30-9194-af0a743d0273
import Infernal

# ╔═╡ ccc08755-a3ca-43c3-bf3d-ebf6a32afa2c
import SamApp2024

# ╔═╡ f05273b2-c9b4-4b15-96e0-bdc3f45dfd63
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 5492b03c-c32e-4471-9748-ec46b7e1dbaa
import Rfam

# ╔═╡ 8616a37c-9f23-41c8-85b5-654b6ec7689b
import StatsBase

# ╔═╡ 061e8cb0-6c3a-4e1a-93d3-eeb400fca037
md"# Load data"

# ╔═╡ b638a0dc-a03a-482e-8479-08da23b40494
Rfam_id = "RF00162"

# ╔═╡ 94c85a6b-8061-43cf-a92b-4e87cac1e147
Rfam_cm = Infernal.cmfetch(Rfam.cm(), Rfam_id).out

# ╔═╡ 17313d40-77c9-426a-a8d8-4946ec12e7a9
SamApp2024_jl_dirname = dirname(Base.current_project())

# ╔═╡ a53943fc-3019-49c2-8a6a-febc55afb0bd
rfamgen_basedir = "/DATA-SSD/cossio/projects/2022/SAM-I_Riboswitch/rfamgen"

# ╔═╡ 9d2c22af-e34c-433e-8915-f7dd2f2214b6
rfamgen_seqs_fasta = "$rfamgen_basedir/outputs/$Rfam_id/sampled_seq.fa"

# ╔═╡ ee4923e5-53c6-4410-a5c0-47e6c57ee4b6
rfamgen_seqs_aln = SamApp2024.infernal_align_fasta_to_cm(rfamgen_seqs_fasta, Rfam_cm)

# ╔═╡ 504ae038-7b2f-4aa3-b0dc-97372b506daa
rfamgen_seqs_raw = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(rfamgen_seqs_fasta))))

# ╔═╡ 1d590af4-a11f-4cc7-890e-2c54c5bebb60
rfamgen_seqs_fasta_sampl = "$rfamgen_basedir/outputs/$Rfam_id/sampled_seq_sampl.fa"

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

# ╔═╡ 4f6dff53-aa8d-4de8-84f9-450e0be3e45d
cm_score_rfamgen = SamApp2024.infernal_score_sequences(Rfam_cm, rfamgen_seqs_raw).bit_sc

# ╔═╡ d65b01bd-c65f-4321-9cc5-c4ba651b00c3
cm_score_rfamgen_sampl = SamApp2024.infernal_score_sequences(Rfam_cm, rfamgen_seqs_raw_sampl).bit_sc

# ╔═╡ c18c8a32-f3fb-4a9b-9ee7-a6d1fa6c2ffc
cm_score_hits = SamApp2024.infernal_score_sequences(Rfam_cm, hits_raw).bit_sc

# ╔═╡ 03413c4f-aab6-4f89-a4c4-d5361f3ce2f1
cm_score_rbm = SamApp2024.infernal_score_sequences(Rfam_cm, [replace(string(seq), '-' => "") for seq = SamApp2024.rnaseq(SamApp2024.rbm2022samples())]).bit_sc

# ╔═╡ 5166ec65-2c64-4f21-84e3-043280318006
rbm = SamApp2024.rbm2022()

# ╔═╡ c4e777fd-5e85-42d2-a72a-b46e84015a50
rfamgen_P_latent_MSA = vec(readdlm("$rfamgen_basedir/outputs/RF00162/log_p_latent_msa.txt"))

# ╔═╡ ad60eadc-a368-40e5-8ed6-6506dd75f48b
rfamgen_P_latent_RBM = vec(readdlm("$rfamgen_basedir/outputs/RF00162/log_p_latent_RBM.txt"))

# ╔═╡ be816ba7-5d00-4f92-91d7-1725e75e2110
rfamgen_P_latent_rCM = vec(readdlm("$rfamgen_basedir/outputs/RF00162/log_p_latent_rCM.txt"))

# ╔═╡ 034ced77-82be-41cc-a2eb-2e23d47fec9c
md"# CM emitted sequences"

# ╔═╡ 5e65934c-ed29-4bb0-a59d-5017dfcceb30
# emit sequences from Rfam CM model
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Rfam_cm; N=5000, aligned=true, outformat="AFA");

# ╔═╡ b2bbe86b-fbed-4eb2-9c2e-44b4b29b6f99
Rfam_cm_emitted_sequences_with_inserts = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)));

# ╔═╡ 2bfc01aa-d4fb-4cfa-be07-617652ae7127
# remove inserts
Rfam_cm_emitted_sequences = LongRNA{4}.([filter(!=('.'), filter(!islowercase, seq)) for seq = Rfam_cm_emitted_sequences_with_inserts]);

# ╔═╡ 6ccb1164-714c-4743-bfba-60201d288602
Rfam_cm_emitted_sequences_infernal_scores = SamApp2024.infernal_score_sequences(Rfam_cm, [replace(string(seq), '-' => "") for seq = Rfam_cm_emitted_sequences]; informat="FASTA", notrunc=false).bit_sc

# ╔═╡ 53f710ec-9193-43c7-ad0c-14b6af9ac892
md"# Intermediate plots"

# ╔═╡ 6610e810-6d94-4a8e-b3e8-7308ee6100de
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="RBM score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(hits_sequences)); normalization=:pdf, label="MSA", bins=200:2:370, color=(:green, 0.5))
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln_sampl)); normalization=:pdf, label="RfamGen", bins=200:2:370, color=:brown, linewidth=2)
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln)); normalization=:pdf, label="RfamGen (argmax)", bins=200:2:370, color=:brown, linewidth=2, linestyle=:dash)
	Makie.axislegend(ax, position=:lt)
	fig
end

# ╔═╡ 048284c6-3a38-4196-8da7-ba55dea5bee0
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400, xlabel="Infernal bit score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, cm_score_hits; normalization=:pdf, label="MSA", bins=30:120, color=(:green, 0.5))
	Makie.stephist!(ax, cm_score_rfamgen_sampl; normalization=:pdf, label="RfamGen", bins=30:120, color=:brown, linewidth=3)
	Makie.stephist!(ax, cm_score_rfamgen; normalization=:pdf, label="RfamGen (argmax)", bins=30:120, color=:brown, linewidth=2, linestyle=:dash)
	Makie.axislegend(ax; position=:lt)
	fig
end

# ╔═╡ 3a2cdb94-1c0a-4798-a14d-80beb09e7716
md"# Final figure"

# ╔═╡ 2b90f592-5a28-485e-84bf-232d900f1abc
let fig = Makie.Figure()
	_sz = 300
	
	ax = Makie.Axis(fig[1,1]; width=_sz, height=_sz, xlabel="RBM score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(hits_sequences)); normalization=:pdf, label="MSA", bins=200:2:370, color=(:gray, 0.5))
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln_sampl)); normalization=:pdf, label="RfamGen", bins=200:2:400, color=:brown, linewidth=2)
	#Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(rfamgen_seqs_aln)); normalization=:pdf, label="RfamGen (argmax)", bins=200:2:400, color=:brown, linewidth=2, linestyle=:dash)
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.rbm2022samples()); normalization=:pdf, label="RBM", bins=200:2:370, color=:blue, linewidth=2)
	Makie.stephist!(ax, -RBMs.free_energy(rbm, SamApp2024.onehot(Rfam_cm_emitted_sequences)); normalization=:pdf, label="rCM", bins=100:2:400, color=:red, linewidth=2)
	Makie.axislegend(ax; position=:lt, framevisible=false)
	Makie.xlims!(ax, 120, 370)
	Makie.ylims!(ax, 0, 0.04)

	ax = Makie.Axis(fig[1,2]; width=_sz, height=_sz, xlabel="CMVAE latent score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, rfamgen_P_latent_MSA; normalization=:pdf, label="MSA", color=(:gray, 0.5), bins=-70:1:0)
	Makie.stephist!(ax, rfamgen_P_latent_rCM; normalization=:pdf, label="rCM", color=:red, linewidth=2, bins=-70:1:0)
	Makie.stephist!(ax, rfamgen_P_latent_RBM; normalization=:pdf, label="RBM", color=:blue, linewidth=2, bins=-70:1:0)
	Makie.stephist!(ax, -dropdims(sum(abs2, randn(16, length(rfamgen_P_latent_RBM)); dims=1); dims=1); normalization=:pdf, label="RfamGen", color=:brown, linewidth=2, bins=-70:1:0)

	ax = Makie.Axis(fig[1,3]; width=_sz, height=_sz, xlabel="Rfam CM bit score", ylabel="frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, cm_score_hits; normalization=:pdf, label="MSA", bins=30:2:150, color=(:gray, 0.5))
	Makie.stephist!(ax, cm_score_rfamgen_sampl; normalization=:pdf, label="RfamGen", bins=30:2:150, color=:brown, linewidth=2)
	#Makie.stephist!(ax, cm_score_rfamgen; normalization=:pdf, label="RfamGen (argmax)", bins=30:2:150, color=:brown, linewidth=2, linestyle=:dash)
	Makie.stephist!(ax, cm_score_rbm; normalization=:pdf, label="RBM", bins=30:2:120, color=:blue, linewidth=2)
	Makie.stephist!(ax, Rfam_cm_emitted_sequences_infernal_scores; normalization=:pdf, label="rCM", bins=-50:2:150, color=:red, linewidth=2)
	Makie.xlims!(ax, -1, 120)
	Makie.ylims!(ax, 0, 0.077)

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,2][1, 1, Makie.TopLeft()], "B)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,3][1, 1, Makie.TopLeft()], "C)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	Makie.save("$SamApp2024_jl_dirname/pluto/SI/Figures/fig-SI_RfamGen.pdf", fig)
	fig
end

# ╔═╡ 3768ddcd-57fb-4fc1-85fd-898c9b6193e1
L"\log(P_{\mathrm{CMVAE}}(z)) + \mathrm{const.}"

# ╔═╡ ce7a137e-7fe9-4cdf-a827-e6db6d93d069
md"# Figure using Denoised CM"

# ╔═╡ Cell order:
# ╠═e98903a9-5faf-4cc0-ab42-a7a7bd5ccfc9
# ╠═c9ddb318-70db-4592-afde-80772aaa2244
# ╠═7327676a-736e-498a-b80e-5d74dc163a35
# ╠═07d606f2-0dc0-470b-83ed-8d0eb1a1a2be
# ╠═f847fe47-b2b5-46c6-bf43-41ee3887c56b
# ╠═8962064d-dfda-4a30-9194-af0a743d0273
# ╠═ccc08755-a3ca-43c3-bf3d-ebf6a32afa2c
# ╠═f05273b2-c9b4-4b15-96e0-bdc3f45dfd63
# ╠═5492b03c-c32e-4471-9748-ec46b7e1dbaa
# ╠═8616a37c-9f23-41c8-85b5-654b6ec7689b
# ╠═794b3eb0-55b8-4f5c-ba03-675e23d2c471
# ╠═ef0e0da8-d64a-4496-94a5-74d683b8cf4d
# ╠═f2ccacb2-dab7-4054-b25d-7075f8a20580
# ╠═c7e6ceff-23f5-4ad5-be6f-a968f5bdf030
# ╠═cb1b7457-21be-4c1e-9363-6258f4b9377a
# ╠═2e28c501-fdca-46ee-b354-309d632c24bf
# ╠═061e8cb0-6c3a-4e1a-93d3-eeb400fca037
# ╠═b638a0dc-a03a-482e-8479-08da23b40494
# ╠═94c85a6b-8061-43cf-a92b-4e87cac1e147
# ╠═17313d40-77c9-426a-a8d8-4946ec12e7a9
# ╠═a53943fc-3019-49c2-8a6a-febc55afb0bd
# ╠═9d2c22af-e34c-433e-8915-f7dd2f2214b6
# ╠═ee4923e5-53c6-4410-a5c0-47e6c57ee4b6
# ╠═504ae038-7b2f-4aa3-b0dc-97372b506daa
# ╠═1d590af4-a11f-4cc7-890e-2c54c5bebb60
# ╠═c036fba4-baff-4e8c-973b-fb9f059251ae
# ╠═fc2d1581-84d8-4812-b596-250b781ff6c2
# ╠═ff8c1333-50f5-40c0-bd4f-1cb2f7653699
# ╠═f2421b0c-df06-4f73-9993-1ce54f96c3cf
# ╠═4ac19098-8342-4bef-844f-969462eb73e6
# ╠═4f6dff53-aa8d-4de8-84f9-450e0be3e45d
# ╠═d65b01bd-c65f-4321-9cc5-c4ba651b00c3
# ╠═c18c8a32-f3fb-4a9b-9ee7-a6d1fa6c2ffc
# ╠═03413c4f-aab6-4f89-a4c4-d5361f3ce2f1
# ╠═5166ec65-2c64-4f21-84e3-043280318006
# ╠═c4e777fd-5e85-42d2-a72a-b46e84015a50
# ╠═ad60eadc-a368-40e5-8ed6-6506dd75f48b
# ╠═be816ba7-5d00-4f92-91d7-1725e75e2110
# ╠═034ced77-82be-41cc-a2eb-2e23d47fec9c
# ╠═5e65934c-ed29-4bb0-a59d-5017dfcceb30
# ╠═b2bbe86b-fbed-4eb2-9c2e-44b4b29b6f99
# ╠═2bfc01aa-d4fb-4cfa-be07-617652ae7127
# ╠═6ccb1164-714c-4743-bfba-60201d288602
# ╠═53f710ec-9193-43c7-ad0c-14b6af9ac892
# ╠═6610e810-6d94-4a8e-b3e8-7308ee6100de
# ╠═048284c6-3a38-4196-8da7-ba55dea5bee0
# ╠═3a2cdb94-1c0a-4798-a14d-80beb09e7716
# ╠═2b90f592-5a28-485e-84bf-232d900f1abc
# ╠═3768ddcd-57fb-4fc1-85fd-898c9b6193e1
# ╠═ce7a137e-7fe9-4cdf-a827-e6db6d93d069

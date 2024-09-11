### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 99970519-3d0b-4a0f-9882-3cf521b4d3b4
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 69a7f451-75c5-4945-a3b5-d6a3f7a14098
using BioSequences: LongRNA

# ╔═╡ d4f706c2-068e-4b94-8514-26a7ab451e3e
using Distributions: Gamma

# ╔═╡ 8161c255-8dca-475a-b946-a57da4a6056e
using Makie: @L_str

# ╔═╡ 29084884-9839-4e2a-83fe-ddaf36c851b9
using NaNStatistics: nanmean

# ╔═╡ 9a51ab6f-8a32-485b-9902-58cb5f62b6e2
using NaNStatistics: nanstd

# ╔═╡ 57780f8a-afae-4034-956b-3b27e028133e
using NaNStatistics: nansum

# ╔═╡ cbc7bac9-aa75-4df2-be8d-fc92e061f421
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 59c353e9-8112-4605-83ed-7bf75d11362f
using Statistics: cor

# ╔═╡ 1acf35b0-29be-48d1-824f-f18d223efc6a
using Statistics: mean

# ╔═╡ 58d6e0ee-b4df-4c1d-9230-96bd06981c0e
using DataFrames: DataFrame

# ╔═╡ a716a269-d78e-4dc6-9101-91d098254ed9
md"# Imports"

# ╔═╡ b88f1dbd-2822-4abe-a205-9a5dde5286be
import PlutoUI

# ╔═╡ e58a23b2-bc2c-474f-bde2-de0b7854ff2f
import CairoMakie

# ╔═╡ 122b7db3-64d3-44cf-9ebd-0204c7066b99
import Makie

# ╔═╡ da01b751-8496-4300-ae0b-1a5d29caec09
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 55f08712-3972-46a5-b77d-19d2e1c7ef4c
import SamApp2024

# ╔═╡ 3698f5ae-4e80-4d6a-9ed1-bd84acb2198d
import StatsBase

# ╔═╡ f88863fc-cb59-4736-a343-9acdb6da3069
import CSV

# ╔═╡ 20b6134e-0d99-485c-a48e-289a2b9c998a
import FASTX

# ╔═╡ d643f048-a536-4069-a7bb-139c84729930
import Infernal

# ╔═╡ 95c19cf9-9f07-4127-9117-34a87783f192
import Rfam

# ╔═╡ 9ce0de3c-8b3a-4222-9428-61ecd3846af8
PlutoUI.TableOfContents()

# ╔═╡ 5f053617-0a7c-4b97-a56d-7c0edc789f83
md"# Functions"

# ╔═╡ b2b836e9-a112-4315-9e0d-d9675c44e0cb
function save_fasta(sequences::BitArray{3})
	@assert size(sequences, 1) == 5
	@assert size(sequences, 2) == 108

	# write to FASTA
	_tmpfasta = tempname()
	FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
		for (n, seq) in enumerate(SamApp2024.rnaseq(sequences))
			@assert !ismissing(seq)
			write(writer, FASTX.FASTA.Record(string(n), filter(!=('-'), string(seq))))
		end
	end

	return _tmpfasta
end

# ╔═╡ 3d67fe2d-eabd-4813-8753-9b39e7edb77e
function save_aligned_fasta(sequences::BitArray{3})
	@assert size(sequences, 1) == 5
	@assert size(sequences, 2) == 108

	# write to FASTA
	_tmpfasta = tempname()
	FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
		for (n, seq) in enumerate(SamApp2024.rnaseq(sequences))
			@assert !ismissing(seq)
			write(writer, FASTX.FASTA.Record(string(n), string(seq)))
		end
	end

	return _tmpfasta
end

# ╔═╡ 714fd72f-6e28-4869-bd5e-060e4bc376e2
md"# RBM samples"

# ╔═╡ 3d5b1a98-fe54-4db3-85b9-54044faf373d
# use saved RBM samples
sampled_v = SamApp2024.rbm2022samples(); # faster

# ╔═╡ 651292a0-886b-4721-83e6-2e6f0710912c
rbm_samples_afa = save_aligned_fasta(sampled_v)

# ╔═╡ d354a48a-3b84-4bc1-8da4-67ee1a7a5688
Infernal.esl_reformat("STOCKHOLM", rbm_samples_afa; informat="AFA")

# ╔═╡ cacff80c-b3b9-4e66-9393-349ad9baa232
md"# Rfam CM samples"

# ╔═╡ 7e5ef2df-172c-42d5-bf15-ed63d3355784
# emit sequences from Rfam CM model
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Infernal.cmfetch(Rfam.cm(), "RF00162").out; N=5000, aligned=true, outformat="AFA")

# ╔═╡ 8b0aa021-7149-46d3-b236-697e7ab868e6
Rfam_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)));

# ╔═╡ 72054d0b-f342-4b2e-b761-ba987e1226e2
Rfam_cm_emitted_sequences_match_only = [join(c for c = seq if isuppercase(c) || c == '-') for seq = Rfam_cm_emitted_sequences]

# ╔═╡ 47888306-48eb-491a-ad5b-bf404d3ed5e2
@assert all(length.(Rfam_cm_emitted_sequences_match_only) .== 108)

# ╔═╡ 314877c5-dd2c-4574-85a0-196b39fb356a
rCM_samples_afa = save_aligned_fasta(SamApp2024.onehot(LongRNA{4}.(Rfam_cm_emitted_sequences_match_only)))

# ╔═╡ b01cbc3c-dfa9-481d-97b4-8ebf55abb777
Infernal.esl_reformat("STOCKHOLM", rCM_samples_afa; informat="AFA")

# ╔═╡ e617c169-ecde-4b48-965c-5a517bf513d0
md"# Denoised CM samples"

# ╔═╡ a4cbd365-0d42-4eba-bc0f-13e03a0d30b7
Denoised_cm = SamApp2024.rfam_RF00162_denoised_cm()

# ╔═╡ 673e12a4-2d42-416a-b1d6-f4d54997013c
Denoised_cm_emitted_sequences_afa = Infernal.cmemit(Denoised_cm.cmout; N=5000, aligned=true, outformat="AFA");

# ╔═╡ 7ef29ee6-0a8d-42f8-abd1-412f4420ff31
Denoised_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Denoised_cm_emitted_sequences_afa.out)));

# ╔═╡ 3e7b9557-6573-485a-a9ae-47073b32bdbe
Denoised_cm_emitted_sequences_match_only = [join(c for c = seq if isuppercase(c) || c == '-') for seq = Denoised_cm_emitted_sequences]

# ╔═╡ 7ad10b49-e12c-4c9f-b864-71fc9196b8e5
@assert all(length.(Denoised_cm_emitted_sequences_match_only) .== 108)

# ╔═╡ 30f86d0a-d875-45de-b75a-a4f8d29b0b17
dCM_samples_afa = save_aligned_fasta(SamApp2024.onehot(LongRNA{4}.(Denoised_cm_emitted_sequences_match_only)))

# ╔═╡ Cell order:
# ╠═a716a269-d78e-4dc6-9101-91d098254ed9
# ╠═99970519-3d0b-4a0f-9882-3cf521b4d3b4
# ╠═b88f1dbd-2822-4abe-a205-9a5dde5286be
# ╠═e58a23b2-bc2c-474f-bde2-de0b7854ff2f
# ╠═122b7db3-64d3-44cf-9ebd-0204c7066b99
# ╠═da01b751-8496-4300-ae0b-1a5d29caec09
# ╠═55f08712-3972-46a5-b77d-19d2e1c7ef4c
# ╠═3698f5ae-4e80-4d6a-9ed1-bd84acb2198d
# ╠═f88863fc-cb59-4736-a343-9acdb6da3069
# ╠═20b6134e-0d99-485c-a48e-289a2b9c998a
# ╠═d643f048-a536-4069-a7bb-139c84729930
# ╠═95c19cf9-9f07-4127-9117-34a87783f192
# ╠═69a7f451-75c5-4945-a3b5-d6a3f7a14098
# ╠═d4f706c2-068e-4b94-8514-26a7ab451e3e
# ╠═8161c255-8dca-475a-b946-a57da4a6056e
# ╠═29084884-9839-4e2a-83fe-ddaf36c851b9
# ╠═9a51ab6f-8a32-485b-9902-58cb5f62b6e2
# ╠═57780f8a-afae-4034-956b-3b27e028133e
# ╠═cbc7bac9-aa75-4df2-be8d-fc92e061f421
# ╠═59c353e9-8112-4605-83ed-7bf75d11362f
# ╠═1acf35b0-29be-48d1-824f-f18d223efc6a
# ╠═58d6e0ee-b4df-4c1d-9230-96bd06981c0e
# ╠═9ce0de3c-8b3a-4222-9428-61ecd3846af8
# ╠═5f053617-0a7c-4b97-a56d-7c0edc789f83
# ╠═b2b836e9-a112-4315-9e0d-d9675c44e0cb
# ╠═3d67fe2d-eabd-4813-8753-9b39e7edb77e
# ╠═714fd72f-6e28-4869-bd5e-060e4bc376e2
# ╠═3d5b1a98-fe54-4db3-85b9-54044faf373d
# ╠═651292a0-886b-4721-83e6-2e6f0710912c
# ╠═d354a48a-3b84-4bc1-8da4-67ee1a7a5688
# ╠═cacff80c-b3b9-4e66-9393-349ad9baa232
# ╠═7e5ef2df-172c-42d5-bf15-ed63d3355784
# ╠═8b0aa021-7149-46d3-b236-697e7ab868e6
# ╠═72054d0b-f342-4b2e-b761-ba987e1226e2
# ╠═47888306-48eb-491a-ad5b-bf404d3ed5e2
# ╠═314877c5-dd2c-4574-85a0-196b39fb356a
# ╠═b01cbc3c-dfa9-481d-97b4-8ebf55abb777
# ╠═e617c169-ecde-4b48-965c-5a517bf513d0
# ╠═a4cbd365-0d42-4eba-bc0f-13e03a0d30b7
# ╠═673e12a4-2d42-416a-b1d6-f4d54997013c
# ╠═7ef29ee6-0a8d-42f8-abd1-412f4420ff31
# ╠═3e7b9557-6573-485a-a9ae-47073b32bdbe
# ╠═7ad10b49-e12c-4c9f-b864-71fc9196b8e5
# ╠═30f86d0a-d875-45de-b75a-a4f8d29b0b17

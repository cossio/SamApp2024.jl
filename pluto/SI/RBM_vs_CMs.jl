### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ f4e57f90-2c1f-49af-a246-0c6a765e1a12
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ b8692df8-26f5-4de4-9cbe-0a248cd55f59
using BioSequences: LongRNA

# ╔═╡ 50f1ea9c-c14e-4d8f-9048-c6647d8abee9
using Distributions: Gamma

# ╔═╡ ce363346-8eca-4e7d-9d61-47446db3cd67
using Makie: @L_str

# ╔═╡ 77892981-16dc-4402-b15a-3ab205ae93a4
using NaNStatistics: nanmean

# ╔═╡ e06158a1-f844-4e43-bfab-b8d4e0beac77
using NaNStatistics: nanstd

# ╔═╡ 52c1d889-080f-42b2-8cbc-a94be40ea43e
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ 2323766a-00fb-4226-bb57-77884fa3c43a
using Statistics: cor

# ╔═╡ dec068c4-eab2-4d05-bc9b-c19e7e8ee626
using Statistics: mean

# ╔═╡ 67082dd1-db8e-40e6-ae83-8029eb1221b1
using NaNStatistics: nansum

# ╔═╡ a9e38c8d-19ed-44a4-8374-a6bf8e857fb2
md"""
# Imports
"""

# ╔═╡ 7dd53826-d7c0-4bb3-93d4-20a3a87da167
import PlutoUI

# ╔═╡ dd79b460-6448-45c8-acb7-8cf0b7b03d3b
import CairoMakie

# ╔═╡ 6a1c2140-2b8b-4f58-ae4e-8f8dc473c089
import Makie

# ╔═╡ 7bc08d72-e56c-4467-94a0-582dbded9ad4
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 86c7fcda-2d67-4b7f-a735-ce300c0227ff
import SamApp2024

# ╔═╡ b374aa45-88b5-413e-afb8-37cc97c9faaa
import StatsBase

# ╔═╡ 67d956b2-18c9-4e12-9d9e-e4a54b056f16
import Infernal

# ╔═╡ 72bb655b-e9db-4834-8d16-7aa51c6bba95
import Rfam

# ╔═╡ 7829c3b2-83d3-47bc-bbf7-e31751031bcd
import FASTX

# ╔═╡ bb893d2b-25aa-49ef-bd50-a2c6ada2448c
PlutoUI.TableOfContents()

# ╔═╡ 85c19ce6-4096-4331-aec2-082cc4c5f24a
md"""
# Load data
"""

# ╔═╡ 3e66fef9-612a-4a9b-89aa-eefb2b82cd8c
# RBM samples
sampled_v = SamApp2024.rbm2022samples(); # faster

# ╔═╡ 066e9cd6-afc9-49f9-8bcf-e17d94dae992
# CM model from Rfam (this has the noisy floor!)
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162");

# ╔═╡ fe359d48-38d6-40b8-bc17-4646f8d154a2
RF00162_seed_stk = Infernal.esl_afetch(Rfam.seed(), "RF00162")

# ╔═╡ 35ff35dd-ec34-4080-9831-45ac5d044356
RF00162_seed_match_cols = findall(≠('.'), SamApp2024.stockholm_ss(RF00162_seed_stk.out));

# ╔═╡ 7db1e9ea-d181-419c-9b02-c30b77c6e197
RF00162_seed_afa = Infernal.esl_reformat("AFA", RF00162_seed_stk.out; informat="STOCKHOLM") # WARNING: this has inserts marked as '-'

# ╔═╡ cdd43222-a6b0-426b-96d0-407d5b7106c8
RF00162_seed_records = collect(FASTX.FASTA.Reader(open(RF00162_seed_afa.out)))

# ╔═╡ 66617a9a-b345-4ad8-8522-80665c81c66d
RF00162_seed_seqs_noinserts = LongRNA{4}.([FASTX.sequence(record)[RF00162_seed_match_cols] for record in RF00162_seed_records]);

# ╔═╡ 9b352035-98c7-416f-8a9f-5c47c62ebda8
# trimmed (no inserts) aligned fasta
RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");

# ╔═╡ cad7df2d-d465-45d2-9ce6-4737adb18fcb
# these are already aligned and without inserts
RF00162_hits_sequences = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa.out))));

# ╔═╡ 7bd2c5f8-16bf-4d98-aa19-96b97fc7332f
# aligned hits, used to train a new noiseless CM model (in Stockholm format, without inserts!)
RF00162_hits_stk = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true);

# ╔═╡ fbce70a0-3b12-4ce6-87e3-34bfbd884382
# fit new CM model using full alignment (without inserts), and without entropic noise
Denoised_cm = Infernal.cmbuild(RF00162_hits_stk.out; enone=true);

# ╔═╡ cec2e5cc-2ebf-41a5-b678-49fff9939131
# emit sequences from Rfam CM model
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Rfam_cm.out; N=5000, aligned=true, outformat="AFA");

# ╔═╡ 8de1da12-cf86-4b5a-81b5-0f4c12724f20
Rfam_cm_emitted_sequences_with_inserts = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)));

# ╔═╡ c9afcea6-36b8-4574-ac6e-6434fbfcbf2e
# remove inserts
Rfam_cm_emitted_sequences = LongRNA{4}.([filter(!=('.'), filter(!islowercase, seq)) for seq = Rfam_cm_emitted_sequences_with_inserts]);

# ╔═╡ b68ea3ea-8b92-4b71-9b6b-f5172efc3f14
# emit sequences from Denoised CM model
Denoised_cm_emitted_sequences_afa = Infernal.cmemit(Denoised_cm.cmout; N=5000, aligned=true, outformat="AFA");

# ╔═╡ 12fe5491-b37a-4052-95d8-19f1d53e6799
Denoised_cm_emitted_sequences_with_inserts = FASTX.sequence.(FASTX.FASTA.Reader(open(Denoised_cm_emitted_sequences_afa.out)))

# ╔═╡ 4611b398-5030-4d0d-bf5c-9b41d57967d4
# remove inserts
Denoised_cm_emitted_sequences = LongRNA{4}.([filter(!=('.'), filter(!islowercase, seq)) for seq in Denoised_cm_emitted_sequences_with_inserts])

# ╔═╡ 9ccdf111-1d2b-41e1-b330-5eb096f5f901
md"""
# Untangled CM
"""

# ╔═╡ 457541e2-2638-4fe1-a9ba-1b58d96da6c0
# Consensus secondary structure of RF00162 in WUSS format
wuss = SamApp2024.rfam_ss("RF00162")

# ╔═╡ f5d1965d-21b8-4395-ab37-6d59a7d5658f
perm = [
    1:findlast(==('A'), wuss); # first segment up to first branch of pseudoknot
    findfirst(==('a'), wuss):findlast(==('a'), wuss); # second branch of pseudoknot
    (findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1); # what's in between the pseudoknot
    (findlast(==('a'), wuss) + 1):108 # what's after the pseudoknot
];

# ╔═╡ 8dc56958-fac1-4971-ad6b-8f8c8c79fc48
# build training alignment for the Untangled CM
RF00162_hits_stk_permuted = tempname()

# ╔═╡ 01ba4b08-5135-48dc-bfa6-9cd7e8f9ee8e
# RF00162_hits_stk has the hits in Stockholm format, with one line per sequence
open(RF00162_hits_stk_permuted, "w") do file
    for (line_index, line) in enumerate(eachline(RF00162_hits_stk.out))
        if startswith(line, "# STOCKHOLM") || startswith(line, "#=GF") || startswith(line, "#=GS") || isempty(line) || line == "//"
            write(file, line, '\n') # these comment lines we just copy
        elseif startswith(line, "#=GR")
            continue # skip GR annotations
        elseif startswith(line, "#=GC SS_cons") # consensus secondary structure
            _ss = line[41:end]
            @assert _ss  == "((((((((,,,,<<<<<---<<<_____>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<----<<<<<<_____>>>>>>-->))))))))"
            @assert wuss == "((((((((,,,,<<<<<---<<<_AAAA>>>------>>>>><<<<-<<<<<<_______>>>>-->>>>>>,,,<aaaa<<<<<<_____>>>>>>-->))))))))"
            @assert _ss == replace(wuss, 'A' => '_', 'a' => '-')
            # add parenthesis for pseudoknot
            _ss = _ss[1:(findfirst(==('A'), wuss) - 1)] * "((((" * _ss[(findlast(==('A'), wuss) + 1):(findfirst(==('a'), wuss) - 1)] * "))))" * _ss[(findlast(==('a'), wuss) + 1):end]
            @assert length(_ss) == 108
            write(file, line[1:40] * _ss[perm], '\n') # write untangled (permuted) secondary structure
        elseif startswith(line, "#=GC RF") # consensus sequence
            @assert line[41:end] == "cucUuAUcaAGAGgGGcgGAGGGAcuGGCCCuaUGAAgCCcCgGCAACCccccauaauaaggggAaGGUGCcAAuuCCugCcggccauuaaggccgGaaagAUaAgag"
            write(file, line[1:40] * line[41:end][perm], '\n')
        else # this is a sequence line
            @assert length(line[41:end]) == 108
            write(file, line[1:40] * line[41:end][perm], '\n')
        end
    end
end

# ╔═╡ ffc1c63c-6bc2-4ddf-a4e1-47afc8ac2b70
# fit new CM model using full alignment, with untangled pseudoknot (without inserts), and without entropic noise
Untangled_cm_permuted = Infernal.cmbuild(RF00162_hits_stk_permuted; enone=true)

# ╔═╡ 0f3a48b5-863c-4afe-8c5f-912e2b5e7c89
# emit sequences from CM model
Untangled_cm_permuted_emitted_sequences_afa = Infernal.cmemit(Untangled_cm_permuted.cmout; N=5000, aligned=true, outformat="AFA");

# ╔═╡ e8663b57-c9ea-421f-851b-e05e3d7c67e8
Untangled_cm_permuted_emitted_sequences_with_inserts = FASTX.sequence.(FASTX.FASTA.Reader(open(Untangled_cm_permuted_emitted_sequences_afa.out)));

# ╔═╡ 5ac9b927-fcbb-4ab7-903f-2d78ad6e6e3f
# remove inserts
Untangled_cm_permuted_emitted_sequences = LongRNA{4}.([filter(!=('.'), filter(!islowercase, seq)) for seq in Untangled_cm_permuted_emitted_sequences_with_inserts]);

# ╔═╡ 3f5b9040-6cb8-4445-a010-9215d8128456
# Permute back to correct column locations
Untangled_cm_emitted_sequences = [seq[invperm(perm)] for seq in Untangled_cm_permuted_emitted_sequences];

# ╔═╡ 9e2d19b1-581a-4529-8b25-6d537263224b
md"""
# Infernal scores
"""

# ╔═╡ be47ff53-67bc-4065-a675-7260f6bd2805
# Infernal scores of hits, using Rfam CM model
_tmp_rfam_fasta = Infernal.esl_reformat("FASTA", RF00162_hits_afa.out; informat="AFA")

# ╔═╡ 624388c1-f99c-4a57-a07b-9914875c300a
_tmp_rfam_cmalign = Infernal.cmalign(Rfam_cm.out, _tmp_rfam_fasta.out; glob=true, informat="FASTA");

# ╔═╡ 04954f12-aa09-4aea-8c19-889a4a14b2d1
_tmp_rfam_cmalign_df = Infernal.cmalign_parse_sfile(_tmp_rfam_cmalign.sfile);

# ╔═╡ 6abe643c-555f-43e1-b5ef-298e863b422c
RF00162_hits_Rfam_cm_scores = _tmp_rfam_cmalign_df.bit_sc;

# ╔═╡ 699d33dc-474b-4d15-b78c-5ab778c30a34
RF00162_hits_Rfam_cm_scores_2 = SamApp2024.infernal_score_sequences(Rfam_cm.out, [replace(string(seq), '-' => "") for seq = RF00162_hits_sequences]).bit_sc

# ╔═╡ b430b5b5-8a04-422c-970a-3e7715b8e2a0
Makie.scatter(RF00162_hits_Rfam_cm_scores, RF00162_hits_Rfam_cm_scores_2)

# ╔═╡ ce1e1bc1-cfbf-426e-9e03-4ef651fb030f
RF00162_hits_sequences

# ╔═╡ 64f26067-e3c6-4a0b-8fe1-e2b9f090d217
# Infernal scores of hits, using Denoised CM model
_tmp_denoised_fasta = Infernal.esl_reformat("FASTA", RF00162_hits_afa.out; informat="AFA")

# ╔═╡ a0ae0102-383d-4a16-89ff-294fc0e4cd1b
_tmp_denoised_cmalign = Infernal.cmalign(Denoised_cm.cmout, _tmp_denoised_fasta.out; glob=true, informat="FASTA");

# ╔═╡ 54ed110a-6600-4b28-85ec-ac45f588057d
_tmp_denoised_cmalign_df = Infernal.cmalign_parse_sfile(_tmp_denoised_cmalign.sfile);

# ╔═╡ 3f1c368a-d7e8-4eb7-aeae-0b3057e47386
RF00162_hits_Denoised_cm_scores = _tmp_denoised_cmalign_df.bit_sc;

# ╔═╡ 7f5553a2-1c9a-4907-96e1-bfc903b6ee7c
# Infernal scores of hits, using Untangled CM model
# first, must permute MSA columns
_tmpafa_perm = tempname()

# ╔═╡ 0319f490-0d0b-4811-9b27-988b723c9257
open(_tmpafa_perm, "w") do io
    for record in FASTX.FASTA.Reader(open(RF00162_hits_afa.out))
        write(io, ">" * FASTX.description(record) * '\n')
        write(io, FASTX.sequence(record)[perm] * '\n')
    end
end

# ╔═╡ c8ea00fb-c0ba-400d-a7e7-66c12aa02b33
_tmp_untangled_fasta = Infernal.esl_reformat("FASTA", _tmpafa_perm; informat="AFA")

# ╔═╡ d69de357-3d0e-439c-876c-11a6f6e99abf
_tmp_untangled_cmalign = Infernal.cmalign(Untangled_cm_permuted.cmout, _tmp_untangled_fasta.out; glob=true, informat="FASTA");

# ╔═╡ cb81f733-7c28-416a-86e2-95afe6c137af
_tmp_untangled_cmalign_df = Infernal.cmalign_parse_sfile(_tmp_untangled_cmalign.sfile);

# ╔═╡ d3a61998-e3fd-45d7-95f3-658d6eb7c554
RF00162_hits_Untangled_cm_scores = _tmp_untangled_cmalign_df.bit_sc;

# ╔═╡ 64c22582-6695-4057-bb93-ea331a24c8af
# Infernal scores of Rfam CM samples
_tmp_rfam_fasta = tempname()
FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
    for (n, seq) in enumerate(Rfam_cm_emitted_sequences)
        ismissing(seq) && continue
        write(writer, FASTX.FASTA.Record(string(n), filter(!=('-'), string(seq))))
    end
end
# unaligned fasta without inserts
_tmp_cmalign = Infernal.cmalign(Rfam_cm.out, _tmpfasta; glob=true, informat="FASTA");
# Infernal scores
_tmp_cmalign_df = Infernal.cmalign_parse_sfile(_tmp_cmalign.sfile);
Rfam_cm_emitted_sequences_infernal_scores = _tmp_cmalign_df.bit_sc;

# ╔═╡ Cell order:
# ╠═a9e38c8d-19ed-44a4-8374-a6bf8e857fb2
# ╠═f4e57f90-2c1f-49af-a246-0c6a765e1a12
# ╠═7dd53826-d7c0-4bb3-93d4-20a3a87da167
# ╠═dd79b460-6448-45c8-acb7-8cf0b7b03d3b
# ╠═6a1c2140-2b8b-4f58-ae4e-8f8dc473c089
# ╠═7bc08d72-e56c-4467-94a0-582dbded9ad4
# ╠═86c7fcda-2d67-4b7f-a735-ce300c0227ff
# ╠═b374aa45-88b5-413e-afb8-37cc97c9faaa
# ╠═67d956b2-18c9-4e12-9d9e-e4a54b056f16
# ╠═72bb655b-e9db-4834-8d16-7aa51c6bba95
# ╠═7829c3b2-83d3-47bc-bbf7-e31751031bcd
# ╠═b8692df8-26f5-4de4-9cbe-0a248cd55f59
# ╠═50f1ea9c-c14e-4d8f-9048-c6647d8abee9
# ╠═ce363346-8eca-4e7d-9d61-47446db3cd67
# ╠═77892981-16dc-4402-b15a-3ab205ae93a4
# ╠═e06158a1-f844-4e43-bfab-b8d4e0beac77
# ╠═52c1d889-080f-42b2-8cbc-a94be40ea43e
# ╠═2323766a-00fb-4226-bb57-77884fa3c43a
# ╠═dec068c4-eab2-4d05-bc9b-c19e7e8ee626
# ╠═67082dd1-db8e-40e6-ae83-8029eb1221b1
# ╠═bb893d2b-25aa-49ef-bd50-a2c6ada2448c
# ╠═85c19ce6-4096-4331-aec2-082cc4c5f24a
# ╠═3e66fef9-612a-4a9b-89aa-eefb2b82cd8c
# ╠═066e9cd6-afc9-49f9-8bcf-e17d94dae992
# ╠═fe359d48-38d6-40b8-bc17-4646f8d154a2
# ╠═35ff35dd-ec34-4080-9831-45ac5d044356
# ╠═7db1e9ea-d181-419c-9b02-c30b77c6e197
# ╠═cdd43222-a6b0-426b-96d0-407d5b7106c8
# ╠═66617a9a-b345-4ad8-8522-80665c81c66d
# ╠═9b352035-98c7-416f-8a9f-5c47c62ebda8
# ╠═cad7df2d-d465-45d2-9ce6-4737adb18fcb
# ╠═7bd2c5f8-16bf-4d98-aa19-96b97fc7332f
# ╠═fbce70a0-3b12-4ce6-87e3-34bfbd884382
# ╠═cec2e5cc-2ebf-41a5-b678-49fff9939131
# ╠═8de1da12-cf86-4b5a-81b5-0f4c12724f20
# ╠═c9afcea6-36b8-4574-ac6e-6434fbfcbf2e
# ╠═b68ea3ea-8b92-4b71-9b6b-f5172efc3f14
# ╠═12fe5491-b37a-4052-95d8-19f1d53e6799
# ╠═4611b398-5030-4d0d-bf5c-9b41d57967d4
# ╠═9ccdf111-1d2b-41e1-b330-5eb096f5f901
# ╠═457541e2-2638-4fe1-a9ba-1b58d96da6c0
# ╠═f5d1965d-21b8-4395-ab37-6d59a7d5658f
# ╠═8dc56958-fac1-4971-ad6b-8f8c8c79fc48
# ╠═01ba4b08-5135-48dc-bfa6-9cd7e8f9ee8e
# ╠═ffc1c63c-6bc2-4ddf-a4e1-47afc8ac2b70
# ╠═0f3a48b5-863c-4afe-8c5f-912e2b5e7c89
# ╠═e8663b57-c9ea-421f-851b-e05e3d7c67e8
# ╠═5ac9b927-fcbb-4ab7-903f-2d78ad6e6e3f
# ╠═3f5b9040-6cb8-4445-a010-9215d8128456
# ╠═9e2d19b1-581a-4529-8b25-6d537263224b
# ╠═be47ff53-67bc-4065-a675-7260f6bd2805
# ╠═624388c1-f99c-4a57-a07b-9914875c300a
# ╠═04954f12-aa09-4aea-8c19-889a4a14b2d1
# ╠═6abe643c-555f-43e1-b5ef-298e863b422c
# ╠═699d33dc-474b-4d15-b78c-5ab778c30a34
# ╠═b430b5b5-8a04-422c-970a-3e7715b8e2a0
# ╠═ce1e1bc1-cfbf-426e-9e03-4ef651fb030f
# ╠═64f26067-e3c6-4a0b-8fe1-e2b9f090d217
# ╠═a0ae0102-383d-4a16-89ff-294fc0e4cd1b
# ╠═54ed110a-6600-4b28-85ec-ac45f588057d
# ╠═3f1c368a-d7e8-4eb7-aeae-0b3057e47386
# ╠═7f5553a2-1c9a-4907-96e1-bfc903b6ee7c
# ╠═0319f490-0d0b-4811-9b27-988b723c9257
# ╠═c8ea00fb-c0ba-400d-a7e7-66c12aa02b33
# ╠═d69de357-3d0e-439c-876c-11a6f6e99abf
# ╠═cb81f733-7c28-416a-86e2-95afe6c137af
# ╠═d3a61998-e3fd-45d7-95f3-658d6eb7c554
# ╠═64c22582-6695-4057-bb93-ea331a24c8af

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

# ╔═╡ 261d7630-deec-4c57-9033-46509090c14f
sampled_distances = hamming_nodiag(sampled_v);

# ╔═╡ e68d577d-ae21-422c-a71e-8af1793ffe47
RF00162_to_sampled_distances = SamApp2024.hamming(SamApp2024.onehot(RF00162_hits_sequences), sampled_v);

# ╔═╡ 7648ec2f-5e08-4ffe-853e-8483d6fd1f35
md"# Site statistics"

# ╔═╡ 6f2074ff-be72-402f-9b76-7b5559c73675
RF00162_joint = reshape(SamApp2024.onehot(RF00162_hits_sequences), 5*108, :) * reshape(SamApp2024.onehot(RF00162_hits_sequences), 5*108, :)' / length(RF00162_hits_sequences)

# ╔═╡ 84cc7005-0aea-424a-8e89-2567add27e34
sampled_joint = reshape(sampled_v, 5*108, :) * reshape(sampled_v, 5*108, :)' / size(sampled_v, 3);

# ╔═╡ a7fe6d86-f213-4d95-b417-c92cf7789a55
RF00162_joint_flat = [reshape(RF00162_joint, 5, 108, 5, 108)[a,i,b,j] for a in 1:5 for b in 1:5 for i in 1:108 for j in 1:108 if i < j];

# ╔═╡ f29202a2-e26c-40bb-8a0c-8a85d2efed64
sampled_joint_flat = [reshape(sampled_joint, 5, 108, 5, 108)[a,i,b,j] for a in 1:5 for b in 1:5 for i in 1:108 for j in 1:108 if i < j];

# ╔═╡ ac611e17-8781-43ac-ad9e-f703ebd57fda
RF00162_joint_cov = RF00162_joint - mean(reshape(SamApp2024.onehot(RF00162_hits_sequences), 5*108, :); dims=2) * mean(reshape(SamApp2024.onehot(RF00162_hits_sequences), 5*108, :); dims=2)';

# ╔═╡ 12ecc8d3-8e3c-4089-8423-d564399673bc
sampled_joint_cov = sampled_joint - mean(reshape(sampled_v, 5*108, :); dims=2) * mean(reshape(sampled_v, 5*108, :); dims=2)';

# ╔═╡ 6e9078ae-1c2c-457e-bf8a-5671c9f16ac1
RF00162_joint_cov_flat = [reshape(RF00162_joint_cov, 5, 108, 5, 108)[a,i,b,j] for a in 1:5 for b in 1:5 for i in 1:108 for j in 1:108 if i < j];

# ╔═╡ 1846a3ff-6e1d-4bd2-9878-2b82fc587651
sampled_joint_cov_flat = [reshape(sampled_joint_cov, 5, 108, 5, 108)[a,i,b,j] for a in 1:5 for b in 1:5 for i in 1:108 for j in 1:108 if i < j];

# ╔═╡ 9dde2f3c-725e-4115-ba02-ac706dbdac9f
joint_kde = KernelDensity.InterpKDE(KernelDensity.kde(hcat(RF00162_joint_flat, sampled_joint_flat)));

# ╔═╡ c2ca4fd8-ece2-4997-832c-f088ca8f0937
joint_cov_kde = KernelDensity.InterpKDE(KernelDensity.kde(hcat(RF00162_joint_cov_flat, sampled_joint_cov_flat)));

# ╔═╡ 46a638e9-8aac-44f1-92f4-6df63b48d1f0
cor(vec(mean(SamApp2024.onehot(RF00162_hits_sequences); dims=3)), vec(mean(sampled_v; dims=3)))

# ╔═╡ a92f025a-dc7a-45ae-84c2-84b3958806b1
cor(RF00162_joint_flat, sampled_joint_flat)

# ╔═╡ ab28ea64-9538-44d6-a428-b4eca24169d4
cor(RF00162_joint_cov_flat, sampled_joint_cov_flat)

# ╔═╡ 7f4324b3-54d2-4242-bc55-4826520edd3f
md"# Figure"

# ╔═╡ 4c9acaa4-1b33-438a-9034-3cf4edc3e4df
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1], width=200, height=200, xlabel="MSA", ylabel="RBM", xticks=0:0.5:1, yticks=0:0.5:1, 
	    #xtrimspine=true, ytrimspine=true,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
		title="conservation (cor. $(round(cor(vec(mean(SamApp2024.onehot(RF00162_hits_sequences); dims=3)), vec(mean(sampled_v; dims=3))); digits=2)))"
	)
	Makie.scatter!(ax, vec(mean(SamApp2024.onehot(RF00162_hits_sequences); dims=3)), vec(mean(sampled_v; dims=3)), color=:red, markersize=7)
	Makie.xlims!(ax, 0, 1)
	Makie.ylims!(ax, 0, 1)
	
	_jnt_idx = rand(1:length(RF00162_joint_cov_flat), 10000)
	
	ax = Makie.Axis(fig[1,2], width=200, height=200, xlabel="sequence length", ylabel="frequency", xticks=60:20:108, yticks=0:0.1:0.15,
	    #xtrimspine=true, ytrimspine=true,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false)
	Makie.hist!(ax, vec(108 .- sum(SamApp2024.onehot(RF00162_hits_sequences)[5,:,:]; dims=1)), color=:gray, normalization=:pdf, label="MSA", bins=60:1:108,
	    strokewidth=0)
	Makie.stephist!(ax, vec(108 .- sum(sampled_v[5,:,:]; dims=1)), color=:red, normalization=:pdf, linewidth=3, label="RBM", bins=60:1:108)
	Makie.xlims!(60, 110)
	Makie.ylims!(0, 0.17)
	Makie.axislegend(ax, position=:lt)
	
	ax = Makie.Axis(fig[2,1], width=200, height=200, xlabel="RBM score", ylabel="frequency", xticks=250:50:350, yticks=0:0.02:0.04,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false)
	Makie.hist!(ax, -RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(RF00162_hits_sequences)), color=:gray, normalization=:pdf, label="MSA", bins=250:5:370)
	Makie.stephist!(ax, -RBMs.free_energy(SamApp2024.rbm2022(), sampled_v), color=:red, normalization=:pdf, linewidth=3, label="RBM", bins=250:5:370)
	Makie.xlims!(ax, 250, 370)
	Makie.ylims!(0, 0.04)
	
	ax = Makie.Axis(fig[1,3], width=200, height=200, xlabel="MSA", ylabel="RBM", xticks=-0.1:0.1:0.1, yticks=-0.1:0.1:0.1,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false,
		title="covariation (cor. $(round(cor(RF00162_joint_cov_flat, sampled_joint_cov_flat); digits=2)))"
	)
	_plt = Makie.scatter!(ax, RF00162_joint_cov_flat[_jnt_idx], sampled_joint_cov_flat[_jnt_idx], markersize=5,
	    color=[clamp(KernelDensity.pdf(joint_cov_kde, x,y), 0, 200) for (x,y) in zip(RF00162_joint_cov_flat[_jnt_idx], sampled_joint_cov_flat[_jnt_idx])],
	    colormap=Makie.Reverse(:redsblues)
	)
	# # Makie.scatter!(ax, RF00162_joint, sampled_joint, markersize=3, 
	# #     color=clamp.([KernelDensity.pdf(joint_kde, x,y) for (x,y) in zip(RF00162_joint, sampled_joint)], 0, 2),
	# #     colormap=:reds
	# # )
	# _ptran = 0:0.01:1
	
	# Makie.scatter!(ax, 
	#     vec([x for x in _ptran, y in _ptran]),
	#     vec([y for x in _ptran, y in _ptran]),
	#     markersize=3, 
	#     color=vec(clamp.([KernelDensity.pdf(joint_kde, x,y) for x in _ptran, y in _ptran], 0, 2)),
	#     colormap=Makie.cgrad([:white, :red], [0.1,1])
	# )
	# #Makie.hexbin!(ax, RF00162_joint, sampled_joint)
	Makie.xlims!(ax, -0.2, 0.2)
	Makie.ylims!(ax, -0.2, 0.2)
	Makie.Colorbar(fig[1,4], _plt)
	
	ax = Makie.Axis(fig[2,2], width=200, height=200, xlabel="distance", ylabel="frequency", xticks=0:50:100, yticks=0:0.02:0.04,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false)
	Makie.hist!(ax, RF00162_distances, color=:gray, normalization=:pdf, label="MSA", bins=0:5:100, strokewidth=0, strokecolor=:gray)
	Makie.stephist!(ax, sampled_distances, color=:red, normalization=:pdf, linewidth=3, label="RBM", bins=0:5:108)
	Makie.xlims!(0, 108)
	Makie.ylims!(0, 0.06)
	
	ax = Makie.Axis(fig[2,3], width=200, height=200, xlabel="distance to closest natural", ylabel="frequency", xticks=0:30:60, yticks=0:0.02:0.04,
	    xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false)
	Makie.stephist!(ax, vec(minimum(RF00162_to_sampled_distances, dims=1)), color=:red, normalization=:pdf, linewidth=3, bins=0:5:60)
	Makie.xlims!(0, 60)
	Makie.ylims!(0, 0.06)
	
	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,2][1, 1, Makie.TopLeft()], "B)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,3][1, 1, Makie.TopLeft()], "C)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "D)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,2][1, 1, Makie.TopLeft()], "E)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,3][1, 1, Makie.TopLeft()], "F)", fontsize = 16, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	#Makie.save(joinpath(paper_figs_dir, "rbm_generative_stats.pdf"), fig)
	fig
end

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
# ╠═261d7630-deec-4c57-9033-46509090c14f
# ╠═e68d577d-ae21-422c-a71e-8af1793ffe47
# ╠═7648ec2f-5e08-4ffe-853e-8483d6fd1f35
# ╠═6f2074ff-be72-402f-9b76-7b5559c73675
# ╠═84cc7005-0aea-424a-8e89-2567add27e34
# ╠═a7fe6d86-f213-4d95-b417-c92cf7789a55
# ╠═f29202a2-e26c-40bb-8a0c-8a85d2efed64
# ╠═ac611e17-8781-43ac-ad9e-f703ebd57fda
# ╠═12ecc8d3-8e3c-4089-8423-d564399673bc
# ╠═6e9078ae-1c2c-457e-bf8a-5671c9f16ac1
# ╠═1846a3ff-6e1d-4bd2-9878-2b82fc587651
# ╠═9dde2f3c-725e-4115-ba02-ac706dbdac9f
# ╠═c2ca4fd8-ece2-4997-832c-f088ca8f0937
# ╠═46a638e9-8aac-44f1-92f4-6df63b48d1f0
# ╠═a92f025a-dc7a-45ae-84c2-84b3958806b1
# ╠═ab28ea64-9538-44d6-a428-b4eca24169d4
# ╠═7f4324b3-54d2-4242-bc55-4826520edd3f
# ╠═4c9acaa4-1b33-438a-9034-3cf4edc3e4df

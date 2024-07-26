### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ c8af275e-a111-47f3-bd5f-0a9eb02eb2f4
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ d8752a7e-9514-4f92-a188-3d63de691744
using BioSequences: LongRNA

# ╔═╡ 770f0daa-c34f-4eab-ab0e-3062c33f3767
using Distributions: Gamma

# ╔═╡ 28e9bdbc-8569-4630-8f94-69b619a654e6
using Makie: @L_str

# ╔═╡ fdb15e3f-38af-4197-8bb5-200707b790d7
using NaNStatistics: nanmean

# ╔═╡ 752cac2b-5808-4995-8069-251921a00cff
using NaNStatistics: nanstd

# ╔═╡ 19820a3b-083a-4dd5-8fff-b7b4ad5f7e52
using RestrictedBoltzmannMachines: free_energy

# ╔═╡ baa679f8-51ff-4192-a458-d33f1eb60d0d
using Statistics: cor

# ╔═╡ 05876609-a565-4a67-9482-25c5e563e3b2
using Statistics: mean

# ╔═╡ bf297d64-72ab-4265-b13e-da4edb433583
using NaNStatistics: nansum

# ╔═╡ 9500714b-7e97-4a4c-8751-a57d51b6f2f7
md"""
# Imports
"""

# ╔═╡ f6c287bf-5bf4-4fec-8c8a-12773c86071d
import PlutoUI

# ╔═╡ d4dd6bb0-255e-4940-a3b7-c88ced90d837
import CairoMakie

# ╔═╡ 6e717083-0f4f-4e52-aa8c-e8b80857be95
import Makie

# ╔═╡ 6aa53ea2-cabc-4115-b564-e9bedd9165da
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ bf30ce70-4ebc-4210-98c1-2e53c9c0c807
import SamApp2024

# ╔═╡ 6a19ce1e-38d5-40c1-ae36-6b62f1667820
import StatsBase

# ╔═╡ d4946941-8ca8-423e-ad1a-45c95e8ece75
PlutoUI.TableOfContents()

# ╔═╡ 80600fef-ebe3-4e2d-8f72-0c2b00b26649
md"""
# Load data
"""

# ╔═╡ 781ad678-4bc3-4aec-960b-4e9fcdd11b41
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 2435abce-d06b-44b0-975c-0e879e7ba3ea
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 9e0fc8f3-5323-4016-a4f3-f4c0928c1187
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions));

# ╔═╡ 041fc0e7-faa7-4de5-a66b-77a75400f075
conds_sam_rep0 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep0", "SAMAP_1M7_0-5SAM_5Mg_T30C_rep0", "SAMAP_1M7_1SAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 3f72fc10-2488-4663-9139-f9d565f43e27
conds_mg_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ cfd80d8b-af4b-4dcb-8a06-51bb3a719199
conds_30C_rep0 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep0"], shape_data_rep0.conditions));

# ╔═╡ 7dbd535f-3585-4299-90e9-9eb1e3cc7674
conds_sam_rep45 = identity.(indexin(["SAMAP_1M7_0-1SAM_5Mg_T30C_rep45", "SAMAP_1M7_1SAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ adcf869e-3390-465d-a309-55d2e659a9b6
conds_mg_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_5Mg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ eda80d63-5cc5-4a5f-9b4e-42c48cbdc9a3
conds_30C_rep45 = identity.(indexin(["SAMAP_1M7_noSAM_noMg_T30C_rep45"], shape_data_rep45.conditions));

# ╔═╡ f89661ac-4eb9-4d93-b2ba-704f06db5227
@show conds_sam_rep0 conds_mg_rep0 conds_30C_rep0;

# ╔═╡ af17005f-8e2f-4bf5-8f92-caa82a61445b
@show conds_sam_rep45 conds_mg_rep45 conds_30C_rep45;

# ╔═╡ 76796d46-7818-4dca-ab7b-d932da8cdea2
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ e4270a30-4328-47d6-9c75-c46d6f272569
rbm_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_rbm")

# ╔═╡ 06218d42-b996-4f05-b59a-2a0a15022f4b
inf_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_syn_inf")

# ╔═╡ 1b5083b3-cd35-41d1-bb87-46a57d13e427
full_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_full30")

# ╔═╡ a36ff682-b1b7-434e-b5e4-041b509e3b83
seed_seqs = findall(shape_data_045.aptamer_origin .== "RF00162_seed70")

# ╔═╡ 57ac7759-51fb-4fb7-81b8-d4ca166fd544
nat_seqs = full_seqs ∪ seed_seqs;

# ╔═╡ 53727ce8-c6d4-48ff-ac81-c3bef1c7e93a
md"""
# Compute protection scores
"""

# ╔═╡ 862e0bab-7beb-494d-a2c6-b080b19dfddf
bps_reactivities_rep0 = shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0];

# ╔═╡ e69fd593-b0c3-4e9d-8962-87a158d43689
nps_reactivities_rep0 = shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0];

# ╔═╡ 161f4f92-84c8-4cc4-99a4-efb60708446b
all_reactivities_rep0 = shape_data_rep0.shape_reactivities[:, nat_seqs, conds_sam_rep0];

# ╔═╡ be50a498-689c-48c4-893b-f1649d5ede2b
shape_stats_rep0 = SamApp2024.shape_basepair_log_odds_v4(;
    shape_data = shape_data_rep0,
    paired_reactivities = bps_reactivities_rep0,
    unpaired_reactivities = nps_reactivities_rep0,
    all_reactivities = all_reactivities_rep0,
    only_hq_profile = true, p_thresh = 1e-3, nsamples = 1000
);

# ╔═╡ 3afc091b-a202-4290-847c-f699d54b2857
md"""
# Paired and unpaired histogram datasets
"""

# ╔═╡ 19eb0e8d-cfee-485b-9823-08e415ab4aba
shape_log_odds_bp_mg = filter(isfinite, shape_stats_rep0.shape_log_odds[bps, nat_seqs, conds_mg_rep0])

# ╔═╡ 721c8e96-22f3-4596-87af-bc338a79f49d
shape_log_odds_np_mg = filter(isfinite, shape_stats_rep0.shape_log_odds[nps, nat_seqs, conds_mg_rep0])

# ╔═╡ b7fea855-354f-444b-9f4a-cfc3d01d4915
shape_log_odds_bp_sam = filter(isfinite, shape_stats_rep0.shape_log_odds[bps, nat_seqs, conds_sam_rep0])

# ╔═╡ cf5261cf-4fac-4568-a76c-71896874237f
shape_log_odds_np_sam = filter(isfinite, shape_stats_rep0.shape_log_odds[nps, nat_seqs, conds_sam_rep0])

# ╔═╡ 10c7cb56-a987-4511-a30d-acbfda1f3f58
bps_reactivities_rep0_flat = filter(isfinite, shape_data_rep0.shape_reactivities[bps, nat_seqs, conds_sam_rep0]);

# ╔═╡ 902f9ef3-0aa4-4984-9ee7-d3f6bd714515
nps_reactivities_rep0_flat = filter(isfinite, shape_data_rep0.shape_reactivities[nps, nat_seqs, conds_sam_rep0]);

# ╔═╡ c287fa5f-dd02-4af6-b297-810bdbb96588
md"""
# Plots
"""

# ╔═╡ 52bf9835-1dff-47a1-9a85-e7605dced135
let fig = Makie.Figure()
	for (col, i) = enumerate((1,8))
		
		data_bp = dropdims(mean(rand(shape_log_odds_bp_sam, i, 10000); dims=1); dims=1)
		data_np = dropdims(mean(rand(shape_log_odds_np_sam, i, 10000); dims=1); dims=1)
		
		ax = Makie.Axis(fig[1, col], width=250, height=250, xlabel="Total protection score", ylabel="Frequency", xgridvisible=false, ygridvisible=false, title="$i sites")
		Makie.hist!(ax, data_bp, normalization=:pdf, bins=-5:0.01:10, color=(:blue, 0.3))
		Makie.hist!(ax, data_np, normalization=:pdf, bins=-5:0.01:10, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp, normalization=:pdf, bins=-5:0.05:10, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np, normalization=:pdf, bins=-5:0.05:10, color=:red, linewidth=2, label="unpaired")
		# Makie.stephist!(ax,
		# 	[dropdims(sum(rand(shape_log_odds_bp_sam, 8, 10000); dims=1); dims=1); dropdims(sum(rand(shape_log_odds_np_sam, 8, 10000); dims=1); dims=1)],
		# 	normalization=:pdf, bins=-11:0.05:5, color=:black, gap=-0.01
		# )
		if i == 1
			Makie.xlims!(ax, -3, 0.8)
		else
			Makie.xlims!(ax, -1.7, 0.7)
		end
		Makie.ylims!(ax, 0, 2.7)
		(col == 1) && Makie.axislegend(ax; position=:lt, framevisible=false)
	end
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 95ada45b-9cd3-4a3d-bea5-edb36697f048
let fig = Makie.Figure()
	for (row, i) = enumerate((1,8))

		data_bp = dropdims(mean(rand(bps_reactivities_rep0_flat, i, 10000); dims=1); dims=1)
		data_np = dropdims(mean(rand(nps_reactivities_rep0_flat, i, 10000); dims=1); dims=1)
		
		ax = Makie.Axis(fig[1, row], width=250, height=250, xlabel="Average reactivity", ylabel="Frequency", xgridvisible=false, ygridvisible=false, title="$i sites")
		Makie.hist!(ax, data_bp, normalization=:pdf, bins=-5:0.05:10, color=(:blue, 0.3))
		Makie.hist!(ax, data_np, normalization=:pdf, bins=-5:0.05:10, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp, normalization=:pdf, bins=-5:0.05:10, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np, normalization=:pdf, bins=-5:0.05:10, color=:red, linewidth=2, label="unpaired")
		# Makie.stephist!(ax,
		# 	[dropdims(sum(rand(shape_log_odds_bp_sam, 8, 10000); dims=1); dims=1); dropdims(sum(rand(shape_log_odds_np_sam, 8, 10000); dims=1); dims=1)],
		# 	normalization=:pdf, bins=-11:0.05:5, color=:black, gap=-0.01
		# )
		Makie.xlims!(ax, -0.5, 3)
		Makie.ylims!(ax, 0, 3.2)
		(row == 1) && Makie.axislegend(ax; position=:rt, framevisible=false)
	end
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 57e48547-f1fe-454f-bf9c-7ace090c238e
md"""
# Final figure
"""

# ╔═╡ 3fc412a9-0395-4375-8a15-d1add4333a8d
let fig = Makie.Figure()
	for (col, i) = enumerate((1,8))

		data_bp = dropdims(mean(rand(shape_log_odds_bp_sam, i, 10000); dims=1); dims=1)
		data_np = dropdims(mean(rand(shape_log_odds_np_sam, i, 10000); dims=1); dims=1)
		
		ax = Makie.Axis(fig[1, col], width=250, height=250, xlabel="Total protection score", ylabel="Frequency", xgridvisible=false, ygridvisible=false, title="$i sites")
		Makie.hist!(ax, data_bp, normalization=:pdf, bins=-5:0.01:10, color=(:blue, 0.3))
		Makie.hist!(ax, data_np, normalization=:pdf, bins=-5:0.01:10, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp, normalization=:pdf, bins=-5:0.05:10, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np, normalization=:pdf, bins=-5:0.05:10, color=:red, linewidth=2, label="unpaired")
		# Makie.stephist!(ax,
		# 	[dropdims(sum(rand(shape_log_odds_bp_sam, 8, 10000); dims=1); dims=1); dropdims(sum(rand(shape_log_odds_np_sam, 8, 10000); dims=1); dims=1)],
		# 	normalization=:pdf, bins=-11:0.05:5, color=:black, gap=-0.01
		# )
		if i == 1
			Makie.xlims!(ax, -3, 0.8)
		else
			Makie.xlims!(ax, -1.7, 0.7)
		end
		Makie.ylims!(ax, 0, 2.7)
		(col == 1) && Makie.axislegend(ax; position=:lt, framevisible=false)
		
		####
		# Second row
		####
		
		data_bp = dropdims(mean(rand(bps_reactivities_rep0_flat, i, 10000); dims=1); dims=1)
		data_np = dropdims(mean(rand(nps_reactivities_rep0_flat, i, 10000); dims=1); dims=1)
		
		ax = Makie.Axis(fig[2, col], width=250, height=250, xlabel="Average reactivity", ylabel="Frequency", xgridvisible=false, ygridvisible=false, title="$i sites")
		Makie.hist!(ax, data_bp, normalization=:pdf, bins=-5:0.01:10, color=(:blue, 0.3))
		Makie.hist!(ax, data_np, normalization=:pdf, bins=-5:0.01:10, color=(:red, 0.3))
		Makie.stephist!(ax, data_bp, normalization=:pdf, bins=-5:0.05:10, color=:blue, linewidth=2, label="paired")
		Makie.stephist!(ax, data_np, normalization=:pdf, bins=-5:0.05:10, color=:red, linewidth=2, label="unpaired")
		# Makie.stephist!(ax,
		# 	[dropdims(sum(rand(shape_log_odds_bp_sam, 8, 10000); dims=1); dims=1); dropdims(sum(rand(shape_log_odds_np_sam, 8, 10000); dims=1); dims=1)],
		# 	normalization=:pdf, bins=-11:0.05:5, color=:black, gap=-0.01
		# )
		Makie.xlims!(ax, -0.5, 3)
		Makie.ylims!(ax, 0, 3.2)
		(col == 1) && Makie.axislegend(ax; position=:rt, framevisible=false)
	end

	Makie.Label(fig[1,1][1, 1, Makie.TopLeft()], "A)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[1,2][1, 1, Makie.TopLeft()], "B)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,1][1, 1, Makie.TopLeft()], "C)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)
	Makie.Label(fig[2,2][1, 1, Makie.TopLeft()], "D)", fontsize = 18, font = :bold, padding = (0, 5, 5, 0), halign = :right)

	Makie.resize_to_layout!(fig)
	Makie.save("/DATA/cossio/SAM/2024/SamApp2024.jl/pluto/SI/Figures/SI_prot_score_additivity_separation.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═9500714b-7e97-4a4c-8751-a57d51b6f2f7
# ╠═c8af275e-a111-47f3-bd5f-0a9eb02eb2f4
# ╠═f6c287bf-5bf4-4fec-8c8a-12773c86071d
# ╠═d4dd6bb0-255e-4940-a3b7-c88ced90d837
# ╠═6e717083-0f4f-4e52-aa8c-e8b80857be95
# ╠═6aa53ea2-cabc-4115-b564-e9bedd9165da
# ╠═bf30ce70-4ebc-4210-98c1-2e53c9c0c807
# ╠═6a19ce1e-38d5-40c1-ae36-6b62f1667820
# ╠═d8752a7e-9514-4f92-a188-3d63de691744
# ╠═770f0daa-c34f-4eab-ab0e-3062c33f3767
# ╠═28e9bdbc-8569-4630-8f94-69b619a654e6
# ╠═fdb15e3f-38af-4197-8bb5-200707b790d7
# ╠═752cac2b-5808-4995-8069-251921a00cff
# ╠═19820a3b-083a-4dd5-8fff-b7b4ad5f7e52
# ╠═baa679f8-51ff-4192-a458-d33f1eb60d0d
# ╠═05876609-a565-4a67-9482-25c5e563e3b2
# ╠═bf297d64-72ab-4265-b13e-da4edb433583
# ╠═d4946941-8ca8-423e-ad1a-45c95e8ece75
# ╠═80600fef-ebe3-4e2d-8f72-0c2b00b26649
# ╠═781ad678-4bc3-4aec-960b-4e9fcdd11b41
# ╠═2435abce-d06b-44b0-975c-0e879e7ba3ea
# ╠═9e0fc8f3-5323-4016-a4f3-f4c0928c1187
# ╠═041fc0e7-faa7-4de5-a66b-77a75400f075
# ╠═3f72fc10-2488-4663-9139-f9d565f43e27
# ╠═cfd80d8b-af4b-4dcb-8a06-51bb3a719199
# ╠═7dbd535f-3585-4299-90e9-9eb1e3cc7674
# ╠═adcf869e-3390-465d-a309-55d2e659a9b6
# ╠═eda80d63-5cc5-4a5f-9b4e-42c48cbdc9a3
# ╠═f89661ac-4eb9-4d93-b2ba-704f06db5227
# ╠═af17005f-8e2f-4bf5-8f92-caa82a61445b
# ╠═76796d46-7818-4dca-ab7b-d932da8cdea2
# ╠═e4270a30-4328-47d6-9c75-c46d6f272569
# ╠═06218d42-b996-4f05-b59a-2a0a15022f4b
# ╠═1b5083b3-cd35-41d1-bb87-46a57d13e427
# ╠═a36ff682-b1b7-434e-b5e4-041b509e3b83
# ╠═57ac7759-51fb-4fb7-81b8-d4ca166fd544
# ╠═53727ce8-c6d4-48ff-ac81-c3bef1c7e93a
# ╠═862e0bab-7beb-494d-a2c6-b080b19dfddf
# ╠═e69fd593-b0c3-4e9d-8962-87a158d43689
# ╠═161f4f92-84c8-4cc4-99a4-efb60708446b
# ╠═be50a498-689c-48c4-893b-f1649d5ede2b
# ╠═3afc091b-a202-4290-847c-f699d54b2857
# ╠═19eb0e8d-cfee-485b-9823-08e415ab4aba
# ╠═721c8e96-22f3-4596-87af-bc338a79f49d
# ╠═b7fea855-354f-444b-9f4a-cfc3d01d4915
# ╠═cf5261cf-4fac-4568-a76c-71896874237f
# ╠═10c7cb56-a987-4511-a30d-acbfda1f3f58
# ╠═902f9ef3-0aa4-4984-9ee7-d3f6bd714515
# ╠═c287fa5f-dd02-4af6-b297-810bdbb96588
# ╠═52bf9835-1dff-47a1-9a85-e7605dced135
# ╠═95ada45b-9cd3-4a3d-bea5-edb36697f048
# ╠═57e48547-f1fe-454f-bf9c-7ace090c238e
# ╠═3fc412a9-0395-4375-8a15-d1add4333a8d

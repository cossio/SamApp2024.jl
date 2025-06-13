### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ e9b17c74-448d-11f0-0ba9-216d9ee5a763
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 90070e34-1811-4e1e-8103-44e5b26d3590
using BioSequences: LongRNA

# ╔═╡ bffcb030-630e-432f-9da7-2b86a49ab4cc
using DataFrames: DataFrame

# ╔═╡ f9b8b730-5a82-46f0-8286-4b2229a1b40a
using Makie: @L_str

# ╔═╡ b301fa5a-9946-4a2f-a010-ec6ccdfd4c99
using NaNStatistics: nansum

# ╔═╡ d1c984b3-ef44-41b9-87a1-da9525ae0085
using Random: bitrand

# ╔═╡ 01a93e83-c0ea-496c-be9d-983552394fc5
using Statistics: cor

# ╔═╡ dfc26c41-175b-42c0-8afa-7701ad7a8773
using Statistics: mean

# ╔═╡ 0d6625fb-5f9c-47f4-b3aa-43c4ee8471f5
using StatsBase: countmap

# ╔═╡ 9e903f2f-3347-4ac9-98a5-043b342086e8
using NaNStatistics: nanmean, nanstd

# ╔═╡ 8ead4fb8-9e99-4206-9d4d-4380efa8c844
md"# Imports"

# ╔═╡ cf3a75a2-d42a-4e23-8ea6-e4ac84084fbb
import Makie, CairoMakie, PlutoUI

# ╔═╡ 6e24bad5-e35e-40da-adc4-129758ef3d2b
import FASTX, Infernal, Rfam, ViennaRNA

# ╔═╡ 86bcc387-c449-4205-aaa2-23681d86e51b
import SamApp2024

# ╔═╡ bfc78de9-71fd-442b-8ad8-defefdc466ea
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 79c3322a-2b3f-4837-86fb-87be3cd786f5
import Unitful

# ╔═╡ 94872098-d780-4cae-9d6e-d59132df85fb
import StatsBase, KernelDensity

# ╔═╡ 9c9cf9c1-203c-40fe-b8ba-f20e4986913f
PlutoUI.TableOfContents()

# ╔═╡ e0f2aa85-e56c-47e6-8096-f6fdfff08087
md"# Load data"

# ╔═╡ f3adfea6-99fd-414c-8529-35e7f97b83a6
_sites = SamApp2024.hallmark_sites_20230507

# ╔═╡ 8c7baad1-d665-463d-8cc6-ca1b744d0257
bps, nps, pks = SamApp2024.RF00162_sites_paired()

# ╔═╡ 963dc356-fb03-4b5f-a7a3-c1ff48cbb5f6
dms_df = SamApp2024.load_dms_data_sequences_table_20250303_with_aligned_sequences();

# ╔═╡ 8f89748c-8da8-46d6-b39b-8d2480e8a3b2
dms_data = SamApp2024.load_dms_data_20250303();

# ╔═╡ ebef5a71-1759-41e6-9ba1-615534a6f033
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20240730_with_pdb();

# ╔═╡ 8f75e840-bca9-4088-abe6-cdedabebbbae
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ 76b10c95-17be-41f0-8e45-1304d5a1d46c
aptamer_names_500 = "APSAM-S2-" .* lpad.(0:499, 3, '0')

# ╔═╡ d50dffa9-837a-4045-a8c4-11e474cd88fc
dms_aptamer_origin = [
	begin
		if n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_rbm"]
			"rbm"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_seed70"]
			"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_full30"]
				"natural"
		elseif n ∈ shape_data_045.aptamer_names[shape_data_045.aptamer_origin .== "RF00162_syn_inf"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "rbm"]
			"rbm"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "infernal"]
			"infernal"
		elseif n ∈ aptamer_names_500[shape_data_500.aptamer_origin .== "Infrared"]
			"Infrared"
		end
	end for n = dms_data.aptamer_names
]

# ╔═╡ 5c25a6e5-674c-4463-b0cd-178fe6aca7db
aptamer_rbm_energies = [
    ismissing(seq) ? missing : 
    RBMs.free_energy(SamApp2024.rbm2022(), SamApp2024.onehot(LongRNA{4}(seq)))
    for seq in dms_data.aligned_sequence
];

# ╔═╡ 0c3e2755-0976-452e-a506-7b9ed175e0fa
nat_seqs = findall(dms_aptamer_origin .== "natural")

# ╔═╡ 1837c99a-5234-4827-99a3-d1dfe02d1732
_rbmlo = findall(dms_aptamer_origin .== "rbm") ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .< -300))

# ╔═╡ e7b9f53e-1cea-4520-8b00-18b693da5f05
_rbmhi = findall(dms_aptamer_origin .== "rbm") ∩ findall((!ismissing).(aptamer_rbm_energies) .&& (aptamer_rbm_energies .> -300))

# ╔═╡ 9d689c9f-bfa7-4fc6-9f74-9dc759f6e45a
conds_sam_dms = [1]

# ╔═╡ 824960b3-2398-4a48-b16e-e3644d91ece2
conds_mg_dms = [2]

# ╔═╡ 765c7d9b-a9ad-4788-b84b-921485e46dbb
ΔR_sam = dms_data.shape_reactivities[:, :, only(conds_sam_dms)] - dms_data.shape_reactivities[:, :, only(conds_mg_dms)];

# ╔═╡ d091f339-57a2-4fb7-b534-56c3cf2754b0
ΔR_sam_avg_natural = [nanmean([ΔR_sam[i,n] for n = nat_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 478ea11c-9d9e-4b89-a020-f37acfd762c5
ΔR_sam_std_natural = [nanstd([ΔR_sam[i,n] for n = nat_seqs if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ c4ee2016-a19c-4008-b6ce-8f96da501664
ΔR_sam_avg_rbmlo = [nanmean([ΔR_sam[i,n] for n = _rbmlo if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 30e40d75-b0c8-4848-8365-6745632046f6
ΔR_sam_std_rbmlo = [nanstd([ΔR_sam[i,n] for n = _rbmlo if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 720d6d0a-3670-4475-a222-08dbd4b6bdc7
ΔR_sam_avg_rbmhi = [nanmean([ΔR_sam[i,n] for n = _rbmhi if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ 21bf6c57-753d-4722-8d45-59787ce516e7
ΔR_sam_std_rbmhi = [nanstd([ΔR_sam[i,n] for n = _rbmhi if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]) for i = 1:108]

# ╔═╡ c968d393-85e1-4abd-9fa6-f8a00f0d5f99
AC_content = [mean(seq[i] ∈ ('A', 'C') for seq = skipmissing(dms_data.aligned_sequence)) for i = 1:108]

# ╔═╡ d129042c-b6e3-4bf8-bf00-7d04bbfb8279
md"# Histograms"

# ╔═╡ f299e745-9842-477f-b76c-28ed8b011a23
let fig = Makie.Figure()
	bp_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = nat_seqs for i = bps for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
	np_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = nat_seqs for i = nps for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]

	pk_nat_sam_reactivities = [dms_data.shape_reactivities[i,n,c] for n = nat_seqs for i = pks for c = conds_sam_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
	pk_nat_mg_reactivities = [dms_data.shape_reactivities[i,n,c] for n = nat_seqs for i = pks for c = conds_mg_dms if !ismissing(dms_data.aligned_sequence[n]) && dms_data.aligned_sequence[n][i] ∈ ('A', 'C')]
	
	ax = Makie.Axis(fig[1,1], width=300, height=300, xlabel="DMS reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:10, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), normalization=:pdf, bins=-1:0.025:2, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), normalization=:pdf, bins=-1:0.025:2, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), label="base paired", normalization=:pdf, bins=-1:0.025:2, linewidth=2, color=:teal)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), label="not paired", normalization=:pdf, bins=-1:0.025:2, linewidth=2, color=:orange)
	Makie.xlims!(-1, 2)
	Makie.ylims!(-0.07, 10)
	#Makie.axislegend(ax, framevisible=false, patchlabelgap=3, position=(-0.02, 1))
	Makie.axislegend(ax, position=(1, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r)
	
	_dummy_ax = Makie.Axis(fig[1,2], width=20, xgridvisible=false, ygridvisible=false)
	Makie.hidespines!(_dummy_ax, :t, :b, :r, :l)
	Makie.hidexdecorations!(_dummy_ax)
	Makie.hideydecorations!(_dummy_ax)
	
	ax = Makie.Axis(fig[1,3], width=300, height=300, xlabel="DMS reactivity", ylabel="density", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), label="b.p.", normalization=:pdf, bins=-1:0.025:2, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), label="n.p.", normalization=:pdf, bins=-1:0.025:2, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, pk_nat_mg_reactivities), label="pseudoknot", normalization=:pdf, bins=-1:0.025:2, linewidth=1, color=:black)
	Makie.xlims!(-1, 2)
	Makie.ylims!(-0.07, 10)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	ax = Makie.Axis(fig[1,4], width=300, height=300, xlabel="DMS reactivity", xgridvisible=false, ygridvisible=false, xticks=-2:1:6, yticks=0:2, xtrimspine=true, ytrimspine=true)
	Makie.hist!(ax, filter(x -> -2 < x < 6, bp_nat_sam_reactivities), normalization=:pdf, bins=-2:0.025:6, color=(:teal, 0.5), gap=-0.01)
	Makie.hist!(ax, filter(x -> -2 < x < 6, np_nat_sam_reactivities), normalization=:pdf, bins=-2:0.025:6, color=(:orange, 0.5), gap=-0.01)
	Makie.stephist!(ax, filter(x -> -2 < x < 6, pk_nat_sam_reactivities), label="pseudoknot", normalization=:pdf, bins=-1:0.025:2, linewidth=1, color=:black)
	Makie.xlims!(-1, 2)
	Makie.ylims!(-0.07, 10)
	Makie.axislegend(ax, position=(0.7, 0.2), framevisible=false)
	Makie.hidespines!(ax, :t, :r, :l)
	Makie.hideydecorations!(ax)
	
	_xs = 1:108

	_std_thresh = 0
	ΔR_sam_avg_natural_no0std = [ΔR_sam_std_natural[i] > _std_thresh ? ΔR_sam_avg_natural[i] : NaN for i = eachindex(ΔR_sam_avg_natural)]
	ΔR_sam_avg_rbmlo_no0std = [ΔR_sam_std_rbmlo[i] > _std_thresh ? ΔR_sam_avg_rbmlo[i] : NaN for i = eachindex(ΔR_sam_std_rbmlo)]
	ΔR_sam_avg_rbmhi_no0std = [ΔR_sam_std_rbmhi[i] > _std_thresh ? ΔR_sam_avg_rbmhi[i] : NaN for i = eachindex(ΔR_sam_avg_rbmhi)]
	
	ax = Makie.Axis(fig[2,:], width=900, height=150, xticks=5:10:108, yticks=-2:0.5:1, xgridvisible=false, ygridvisible=false, ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_natural_no0std - ΔR_sam_std_natural/2)[_xs], (ΔR_sam_avg_natural_no0std + ΔR_sam_std_natural/2)[_xs]; color=(:gray, 0.25))
	Makie.band!(ax, _xs, (ΔR_sam_avg_rbmlo_no0std - ΔR_sam_std_rbmlo/2)[_xs], (ΔR_sam_avg_rbmlo_no0std + ΔR_sam_std_rbmlo/2)[_xs], color=(:blue, 0.25))
	Makie.scatterlines!(ax, _xs, ΔR_sam_avg_natural_no0std; linewidth=0.5, markersize=5, color=:black, label="Natural")
	Makie.scatterlines!(ax, _xs, ΔR_sam_avg_rbmlo_no0std, linewidth=0.5, markersize=5, color=:blue, label="RBM (RBMscore>300)")
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax, :t, :r, :b)
	Makie.hidexdecorations!(ax)
	Makie.xlims!(1, 108)
	Makie.ylims!(-0.7, 0.7)

	ax = Makie.Axis(fig[3,:], width=900, height=150, xticks=5:10:108, yticks=-2:0.5:1, xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="Δreactivity", xtrimspine=true, ytrimspine=true)
	
	Makie.band!(ax, _xs, (ΔR_sam_avg_natural_no0std - ΔR_sam_std_natural/2)[_xs], (ΔR_sam_avg_natural_no0std + ΔR_sam_std_natural/2)[_xs], color=(:gray, 0.25))
	Makie.band!(ax, _xs, (ΔR_sam_avg_rbmhi - ΔR_sam_std_rbmhi/2)[_xs], (ΔR_sam_avg_rbmhi + ΔR_sam_std_rbmhi/2)[_xs], color=(:red, 0.25))
	Makie.scatterlines!(ax, _xs, ΔR_sam_avg_natural_no0std; markersize=5, linewidth=0.5, color=:black, label="Natural")
	Makie.scatterlines!(ax, _xs, ΔR_sam_avg_rbmhi_no0std; markersize=5, linewidth=0.5, color=:red, label="RBM (RBMscore<300)")
	
	Makie.scatter!(ax, SamApp2024.hallmark_sites_20230507, -0.65one.(SamApp2024.hallmark_sites_20230507), markersize=7, color=:black, marker=:utriangle)
	Makie.axislegend(ax, position=(0.5, 0), framevisible=false, patchlabelgap=-3)
	Makie.hidespines!(ax, :t, :r)
	Makie.xlims!(1, 108)
	Makie.ylims!(-0.7, 0.7)

	ax_AC = Makie.Axis(fig[4,:], width=900, height=40, xticks=5:10:108, xgridvisible=false, ygridvisible=false, xlabel="site", ylabel="A,C content", xtrimspine=true, ytrimspine=true)
	Makie.barplot!(ax_AC, AC_content; width=0.5, color=:black)
	Makie.xlims!(0, 109)
	Makie.hidespines!(ax_AC, :t, :r)
	Makie.linkxaxes!(ax, ax_AC)
	
	# Makie.Label(fig[1,1][1,1,Makie.TopLeft()], "A)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,2][1,1,Makie.TopLeft()], "B)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[1,3][1,1,Makie.TopLeft()], "C)", font=:bold, padding=(0,0,10,10))
	# Makie.Label(fig[2,:][1,1,Makie.TopLeft()], "D)", font=:bold, padding=(0,0,0,0))
	# Makie.Label(fig[3,:][1,1,Makie.TopLeft()], "E)", font=:bold, padding=(0,0,0,0))
	
	Makie.resize_to_layout!(fig)
	#Makie.save("/workspaces/SamApp.jl/notebooks/2024-03-14 New paper figures/Figures/SHAPE reactivities.pdf", fig)
	#Makie.save("Figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═8ead4fb8-9e99-4206-9d4d-4380efa8c844
# ╠═e9b17c74-448d-11f0-0ba9-216d9ee5a763
# ╠═cf3a75a2-d42a-4e23-8ea6-e4ac84084fbb
# ╠═6e24bad5-e35e-40da-adc4-129758ef3d2b
# ╠═86bcc387-c449-4205-aaa2-23681d86e51b
# ╠═bfc78de9-71fd-442b-8ad8-defefdc466ea
# ╠═79c3322a-2b3f-4837-86fb-87be3cd786f5
# ╠═94872098-d780-4cae-9d6e-d59132df85fb
# ╠═90070e34-1811-4e1e-8103-44e5b26d3590
# ╠═bffcb030-630e-432f-9da7-2b86a49ab4cc
# ╠═f9b8b730-5a82-46f0-8286-4b2229a1b40a
# ╠═b301fa5a-9946-4a2f-a010-ec6ccdfd4c99
# ╠═d1c984b3-ef44-41b9-87a1-da9525ae0085
# ╠═01a93e83-c0ea-496c-be9d-983552394fc5
# ╠═dfc26c41-175b-42c0-8afa-7701ad7a8773
# ╠═0d6625fb-5f9c-47f4-b3aa-43c4ee8471f5
# ╠═9e903f2f-3347-4ac9-98a5-043b342086e8
# ╠═9c9cf9c1-203c-40fe-b8ba-f20e4986913f
# ╠═e0f2aa85-e56c-47e6-8096-f6fdfff08087
# ╠═f3adfea6-99fd-414c-8529-35e7f97b83a6
# ╠═8c7baad1-d665-463d-8cc6-ca1b744d0257
# ╠═963dc356-fb03-4b5f-a7a3-c1ff48cbb5f6
# ╠═8f89748c-8da8-46d6-b39b-8d2480e8a3b2
# ╠═ebef5a71-1759-41e6-9ba1-615534a6f033
# ╠═8f75e840-bca9-4088-abe6-cdedabebbbae
# ╠═76b10c95-17be-41f0-8e45-1304d5a1d46c
# ╠═d50dffa9-837a-4045-a8c4-11e474cd88fc
# ╠═5c25a6e5-674c-4463-b0cd-178fe6aca7db
# ╠═0c3e2755-0976-452e-a506-7b9ed175e0fa
# ╠═1837c99a-5234-4827-99a3-d1dfe02d1732
# ╠═e7b9f53e-1cea-4520-8b00-18b693da5f05
# ╠═9d689c9f-bfa7-4fc6-9f74-9dc759f6e45a
# ╠═824960b3-2398-4a48-b16e-e3644d91ece2
# ╠═765c7d9b-a9ad-4788-b84b-921485e46dbb
# ╠═d091f339-57a2-4fb7-b534-56c3cf2754b0
# ╠═478ea11c-9d9e-4b89-a020-f37acfd762c5
# ╠═c4ee2016-a19c-4008-b6ce-8f96da501664
# ╠═30e40d75-b0c8-4848-8365-6745632046f6
# ╠═720d6d0a-3670-4475-a222-08dbd4b6bdc7
# ╠═21bf6c57-753d-4722-8d45-59787ce516e7
# ╠═c968d393-85e1-4abd-9fa6-f8a00f0d5f99
# ╠═d129042c-b6e3-4bf8-bf00-7d04bbfb8279
# ╠═f299e745-9842-477f-b76c-28ed8b011a23

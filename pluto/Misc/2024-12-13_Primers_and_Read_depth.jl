### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 107ae327-48d3-46de-9e1d-736736aeffc1
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ f8d39bff-fb97-4c50-8af4-e73b064e7b31
using BioSequences: LongRNA

# ╔═╡ 7edbd7f7-d082-426b-aa03-03c3c940fb43
using DataFrames: DataFrame

# ╔═╡ 84590e4c-575d-46d5-aeff-956bcdeb8cc6
using Distributions: Gamma

# ╔═╡ 0fc7fd33-6230-4498-b759-af5668244d7d
using Distributions: logpdf

# ╔═╡ e4270162-7ff1-4f6b-8d14-554e89053e18
using Distributions: pdf

# ╔═╡ 8f7e2cc8-676b-47d0-8155-b0cdee40e3fc
using Distributions: Poisson

# ╔═╡ a6a708ff-1188-4726-845b-d197b342b69a
using LinearAlgebra: Diagonal

# ╔═╡ 91acc0f5-77bf-4f8a-b934-d55c523f0aa0
using LinearAlgebra: eigen

# ╔═╡ aa76b309-7841-4fc4-88c9-783e82ea3e82
using Makie: @L_str

# ╔═╡ f2c0d959-a7a4-4c8e-921d-ecb49f1cbc6d
using Random: bitrand

# ╔═╡ 44d13119-ddc7-4fad-b66f-3a1ad28a1837
using Statistics: cor

# ╔═╡ bbd4a036-c657-47c7-90de-c170d35afc61
using Statistics: mean

# ╔═╡ 003bc9f5-5a36-4954-b8c5-3d8d672609eb
using StatsBase: countmap

# ╔═╡ 6841dbc9-2d6f-495b-9a20-73d7e6e71fe0
using NaNStatistics: nansum

# ╔═╡ 73fffa9f-730d-4a66-bb02-675bb3ca597d
using NaNStatistics: nanmean

# ╔═╡ 428714d8-b35b-47b7-aa27-b8beae1f48bf
using NaNStatistics: nanstd

# ╔═╡ 7ceb8b04-b944-11ef-020d-4b8a8c874f26
md"# Imports"

# ╔═╡ 031e448b-eb8d-4c02-9a6d-13907e70dd3b
import Makie

# ╔═╡ 710b8e9d-b88b-4d8f-bfd1-135f5cd32d41
import CairoMakie

# ╔═╡ 5db2bb41-706f-49ad-83cc-17f2aca939c0
import CSV

# ╔═╡ aaa594df-7b34-48dc-b286-b6d9f578aea3
import Infernal

# ╔═╡ 6f3da9d5-3916-4938-9730-1517cdd625e8
import KernelDensity

# ╔═╡ 11ce098d-7bc8-4ca3-95e8-69c906f21988
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ eab6becb-89f9-460d-8b15-b9c4aaa076fd
import Rfam

# ╔═╡ d29a7b49-327f-4787-8ca6-84bb0c334aa9
import HDF5

# ╔═╡ d82b2713-560d-467d-8a01-3483ee04da2a
import FASTX

# ╔═╡ 5435d700-b82a-4dd1-9b0d-df4c1d00cb4d
import SamApp2024

# ╔═╡ 5d919a3a-88ca-4940-aea3-936bda2b08b3
import StatsBase

# ╔═╡ 7ebaba22-65e6-4a86-be6b-ef9657d30fb7
import Logomaker

# ╔═╡ 0c041490-88e5-487c-ae41-bd49a81abcca
import PlutoUI

# ╔═╡ ba1af249-b24e-452d-b7d8-f69da56aba70
md"# Notebook setup"

# ╔═╡ 6895d85f-9fe1-48c3-864b-3ecee6af76b2
PlutoUI.TableOfContents()

# ╔═╡ 71f45265-b85d-4b6f-9cea-27329ad40db7
md"# Load data"

# ╔═╡ e879e0df-740c-47c3-a174-e1fa3891aa85
(; bps, nps, pks) = SamApp2024.RF00162_sites_paired()

# ╔═╡ 45167f0e-6b0f-4802-9e55-d865265b6cca
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ 54e7abed-838c-4c8c-825e-5f4daded922b
# split rep0 from rep4+5
shape_data_rep0 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ edbf16e1-338a-4263-9e3c-e8a89284dc7b
# split rep0 from rep4+5
shape_data_rep45 = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep45"), shape_data_045.conditions))

# ╔═╡ e3cc5d11-38bc-4f09-9088-ab70e08da111
shape_data_500 = SamApp2024.load_shapemapper_data_500v2_20240315();

# ╔═╡ a37cb0ef-a8bd-404d-b0cc-92b326cca485
seq_groups_dfs = SamApp2024.artifact_load_sequencing_groups_2024_11_27();

# ╔═╡ f7b6646d-45a9-4cc5-88a2-bc822c8c1c02
foreach(println, sort(collect(keys(seq_groups_dfs))))

# ╔═╡ c771e473-0140-48c2-9453-f3621361edb7
names_500 = ["APSAM-S2-" * lpad(n - 1, 3, "0") for n = 1:500][shape_data_500.aptamer_origin .!= "Infrared"]

# ╔═╡ e6be846b-82a8-46e8-833d-885671f10cef
md"# Plot"

# ╔═╡ 97385709-977d-49fe-8ab2-5d4579b26f28
let fig = Makie.Figure()
	groups = ["GP1-Natural-primer1", "GP2-Natural-primer2", "GP3-Natural-primer3", "GP4-Synthetic-Set1-primer4", "GP5-Synthetic-Set1-primer5"]
	colors = [:orange, :purple, :red, :blue, :teal]

	_width = 1000
	_height = 200
	
	ax = Makie.Axis(fig[1,1]; width=_width, height=_height, xticks=5:5:108, title="Modified condition (Repl.0)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep0.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep0.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)

	ax = Makie.Axis(fig[2,1]; width=_width, height=_height, xticks=5:5:108, title="Unmodified condition (Repl.0)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep0.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep0.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)
	
	ax = Makie.Axis(fig[3,1]; width=_width, height=_height, xticks=5:5:108, title="Denatured condition (Repl.0)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep0.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep0.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 40dcf092-28d9-4eab-a4bb-f4a424c3605c
let fig = Makie.Figure()
	groups = ["GP1-Natural-primer1", "GP2-Natural-primer2", "GP3-Natural-primer3", "GP4-Synthetic-Set1-primer4", "GP5-Synthetic-Set1-primer5"]
	colors = [:orange, :purple, :red, :blue, :teal]
	
	ax = Makie.Axis(fig[1,1]; width=1000, height=200, xticks=5:5:108, title="Modified condition (Repl.45)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep45.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep45.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)

	ax = Makie.Axis(fig[2,1]; width=1000, height=200, xticks=5:5:108, title="Unmodified condition (Repl.45)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep45.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep45.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)
	
	ax = Makie.Axis(fig[3,1]; width=1000, height=200, xticks=5:5:108, title="Denatured condition (Repl.45)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_rep45.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		σ = nanstd(shape_data_rep45.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, shape_data_045.aptamer_names), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ f52628d6-5a08-41a5-83f0-9b0a76cac4b2
let fig = Makie.Figure()
	groups = ["GP6-Synthetic-Set2-primer1", "GP7-Synthetic-Set2-primer2", "GP8-Synthetic-Set2-primer3", "GP9-Synthetic-Set2-primer5"]
	colors = [:orange, :purple, :red, :teal]
	
	ax = Makie.Axis(fig[1,1]; width=1000, height=200, xticks=5:5:108, title="Modified condition (Experiment 2)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_500.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		σ = nanstd(shape_data_500.shape_M_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		_idx = findall(!isnan, μ)
		_x = (1:108)[_idx]
		Makie.band!(ax, _x, (μ - σ/2)[_idx], (μ + σ/2)[_idx]; color=(color, 0.1))
		Makie.scatterlines!(ax, _x, μ[_idx]; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)

	ax = Makie.Axis(fig[2,1]; width=1000, height=200, xticks=5:5:108, title="Unmodified condition (Experiment 2)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_500.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		σ = nanstd(shape_data_500.shape_U_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)
	
	ax = Makie.Axis(fig[3,1]; width=1000, height=200, xticks=5:5:108, title="Denatured condition (Experiment 2)", xlabel="Position         ", ylabel="Read depth")
	for (gr, color) = zip(groups, colors)
		μ = nanmean(shape_data_500.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		σ = nanstd(shape_data_500.shape_D_depth[:, indexin(seq_groups_dfs[gr].name, names_500), :]; dim=(2,3))
		Makie.band!(ax, 1:108, μ - σ/2, μ + σ/2; color=(color, 0.1))
		Makie.scatterlines!(ax, 1:108, μ; label=gr, color)
	end
	Makie.xlims!(ax, 0.5, 140)
	Makie.axislegend(ax)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 98098fea-0939-46c6-b899-914e3402756c
nanmean(shape_data_500.shape_M_depth[:, indexin(seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name, names_500), :]; dim=(2,3))

# ╔═╡ 1f734bf4-98ab-4f8d-880e-721cedcb71c1
nanstd(shape_data_500.shape_M_depth[:, indexin(seq_groups_dfs["GP9-Synthetic-Set2-primer5"].name, names_500), :]; dim=(2,3))

# ╔═╡ Cell order:
# ╠═7ceb8b04-b944-11ef-020d-4b8a8c874f26
# ╠═107ae327-48d3-46de-9e1d-736736aeffc1
# ╠═031e448b-eb8d-4c02-9a6d-13907e70dd3b
# ╠═710b8e9d-b88b-4d8f-bfd1-135f5cd32d41
# ╠═5db2bb41-706f-49ad-83cc-17f2aca939c0
# ╠═aaa594df-7b34-48dc-b286-b6d9f578aea3
# ╠═6f3da9d5-3916-4938-9730-1517cdd625e8
# ╠═11ce098d-7bc8-4ca3-95e8-69c906f21988
# ╠═eab6becb-89f9-460d-8b15-b9c4aaa076fd
# ╠═d29a7b49-327f-4787-8ca6-84bb0c334aa9
# ╠═d82b2713-560d-467d-8a01-3483ee04da2a
# ╠═5435d700-b82a-4dd1-9b0d-df4c1d00cb4d
# ╠═5d919a3a-88ca-4940-aea3-936bda2b08b3
# ╠═7ebaba22-65e6-4a86-be6b-ef9657d30fb7
# ╠═0c041490-88e5-487c-ae41-bd49a81abcca
# ╠═f8d39bff-fb97-4c50-8af4-e73b064e7b31
# ╠═7edbd7f7-d082-426b-aa03-03c3c940fb43
# ╠═84590e4c-575d-46d5-aeff-956bcdeb8cc6
# ╠═0fc7fd33-6230-4498-b759-af5668244d7d
# ╠═e4270162-7ff1-4f6b-8d14-554e89053e18
# ╠═8f7e2cc8-676b-47d0-8155-b0cdee40e3fc
# ╠═a6a708ff-1188-4726-845b-d197b342b69a
# ╠═91acc0f5-77bf-4f8a-b934-d55c523f0aa0
# ╠═aa76b309-7841-4fc4-88c9-783e82ea3e82
# ╠═f2c0d959-a7a4-4c8e-921d-ecb49f1cbc6d
# ╠═44d13119-ddc7-4fad-b66f-3a1ad28a1837
# ╠═bbd4a036-c657-47c7-90de-c170d35afc61
# ╠═003bc9f5-5a36-4954-b8c5-3d8d672609eb
# ╠═6841dbc9-2d6f-495b-9a20-73d7e6e71fe0
# ╠═73fffa9f-730d-4a66-bb02-675bb3ca597d
# ╠═428714d8-b35b-47b7-aa27-b8beae1f48bf
# ╠═ba1af249-b24e-452d-b7d8-f69da56aba70
# ╠═6895d85f-9fe1-48c3-864b-3ecee6af76b2
# ╠═71f45265-b85d-4b6f-9cea-27329ad40db7
# ╠═e879e0df-740c-47c3-a174-e1fa3891aa85
# ╠═45167f0e-6b0f-4802-9e55-d865265b6cca
# ╠═54e7abed-838c-4c8c-825e-5f4daded922b
# ╠═edbf16e1-338a-4263-9e3c-e8a89284dc7b
# ╠═e3cc5d11-38bc-4f09-9088-ab70e08da111
# ╠═a37cb0ef-a8bd-404d-b0cc-92b326cca485
# ╠═f7b6646d-45a9-4cc5-88a2-bc822c8c1c02
# ╠═c771e473-0140-48c2-9453-f3621361edb7
# ╠═e6be846b-82a8-46e8-833d-885671f10cef
# ╠═97385709-977d-49fe-8ab2-5d4579b26f28
# ╠═40dcf092-28d9-4eab-a4bb-f4a424c3605c
# ╠═f52628d6-5a08-41a5-83f0-9b0a76cac4b2
# ╠═98098fea-0939-46c6-b899-914e3402756c
# ╠═1f734bf4-98ab-4f8d-880e-721cedcb71c1

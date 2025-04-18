### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ bc47d34d-3b75-4ddb-b729-892227278f7c
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 7720ea19-c7a9-4b59-afce-9922fd7d2524
using BioSequences: LongRNA

# ╔═╡ c3ec598b-8d33-496f-ab59-e788a39b5b9c
using DataFrames: DataFrame

# ╔═╡ 49196f78-1bc1-11f0-1926-e3f1b0558a8e
md"# Imports"

# ╔═╡ 61adb389-87e1-482c-b220-49dd4b8857b1
import Makie

# ╔═╡ f17e3dff-3817-402d-a750-64ef2eefea78
import CairoMakie

# ╔═╡ a2810e8a-e872-474d-8e2e-ba3d1084a4f0
import FASTX

# ╔═╡ e7b2a1f9-187e-49f4-bab4-3d96139da0ac
import Infernal

# ╔═╡ d7a12a22-4fdc-40b6-8c90-ef244bb3b118
import SamApp2024

# ╔═╡ 7814520b-c308-4588-9e7e-e6f735c5a0d5
import CSV

# ╔═╡ ace0c32f-a902-44db-adbd-fda782886af0
import Rfam

# ╔═╡ 22eed13e-8c38-47fb-ad2b-fadf3dd73638
import PlutoUI

# ╔═╡ 16086e7d-d5f2-4f61-913f-e82b6b2f13ee
PlutoUI.TableOfContents()

# ╔═╡ 952c645a-5458-4d77-9b5d-1911e5db8ca4
md"# Load data"

# ╔═╡ 86e28236-1af6-4a62-9e50-0d9338d174ce
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20240730_with_pdb();

# ╔═╡ 856cd0e1-8374-4a23-b8a2-124f8006a130
pdb_alignment = SamApp2024.shape_pdb_positions_mapping_20240731()[1,:]

# ╔═╡ cece3da9-9b99-433f-aea0-2e131de9e5a8
findall(pdb_alignment .=== 47)

# ╔═╡ 1dda4f09-32f8-4faa-9114-55f0db83d7e9
pdb_alignment

# ╔═╡ 7a1ffc95-9efc-4616-996d-c49f150d6e05
n_ex = only(findall(shape_data_045.aptamer_names .== "SAMAP-PDB0"))

# ╔═╡ f14e7003-ccc6-4117-b4f3-8bf283a04446
@show shape_data_045.aptamer_ids[n_ex]

# ╔═╡ 648b3d7b-1f48-4b88-9409-22a32ece39dd
size(shape_data_045.shape_reactivities[:, n_ex, :])

# ╔═╡ b606262e-8b5c-42c7-9b7e-dd6026620735
aptamer_name = "SAMAP-PDB0"

# ╔═╡ 415c0084-5f76-48c8-9d5d-47688af6a2d0
shape_dir = SamApp2024.shapemapper_data_pierre_demux_20230920_dir(; demux=true)

# ╔═╡ 2174ffe4-e40d-444b-8bfa-91db16694aa4
cond = first(filter(startswith("SAMAP"), readdir(shape_dir)))

# ╔═╡ 4956eb5e-c934-4097-95ec-568cf5372d54
profile_file = joinpath(shape_dir, cond, "$(cond)_$(aptamer_name)_profile.txt")

# ╔═╡ 6ead69f8-3486-4335-adec-55070d398eb4
prof_df = CSV.read(profile_file, DataFrame)

# ╔═╡ Cell order:
# ╠═49196f78-1bc1-11f0-1926-e3f1b0558a8e
# ╠═bc47d34d-3b75-4ddb-b729-892227278f7c
# ╠═61adb389-87e1-482c-b220-49dd4b8857b1
# ╠═f17e3dff-3817-402d-a750-64ef2eefea78
# ╠═a2810e8a-e872-474d-8e2e-ba3d1084a4f0
# ╠═e7b2a1f9-187e-49f4-bab4-3d96139da0ac
# ╠═d7a12a22-4fdc-40b6-8c90-ef244bb3b118
# ╠═7814520b-c308-4588-9e7e-e6f735c5a0d5
# ╠═ace0c32f-a902-44db-adbd-fda782886af0
# ╠═22eed13e-8c38-47fb-ad2b-fadf3dd73638
# ╠═7720ea19-c7a9-4b59-afce-9922fd7d2524
# ╠═c3ec598b-8d33-496f-ab59-e788a39b5b9c
# ╠═16086e7d-d5f2-4f61-913f-e82b6b2f13ee
# ╠═952c645a-5458-4d77-9b5d-1911e5db8ca4
# ╠═86e28236-1af6-4a62-9e50-0d9338d174ce
# ╠═856cd0e1-8374-4a23-b8a2-124f8006a130
# ╠═cece3da9-9b99-433f-aea0-2e131de9e5a8
# ╠═1dda4f09-32f8-4faa-9114-55f0db83d7e9
# ╠═7a1ffc95-9efc-4616-996d-c49f150d6e05
# ╠═f14e7003-ccc6-4117-b4f3-8bf283a04446
# ╠═648b3d7b-1f48-4b88-9409-22a32ece39dd
# ╠═b606262e-8b5c-42c7-9b7e-dd6026620735
# ╠═415c0084-5f76-48c8-9d5d-47688af6a2d0
# ╠═2174ffe4-e40d-444b-8bfa-91db16694aa4
# ╠═4956eb5e-c934-4097-95ec-568cf5372d54
# ╠═6ead69f8-3486-4335-adec-55070d398eb4

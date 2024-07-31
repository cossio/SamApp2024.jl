### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ e1633647-c7e2-46af-af5f-016f247e5d82
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 7bb51d21-270c-4725-a59d-d81bc676070b
md"""
# Imports
"""

# ╔═╡ 14a89917-47de-4a78-8d68-c6ee65205a60
import PlutoUI

# ╔═╡ 9262d7dc-c0f2-4711-89ad-570dfbf6e322
import Rfam

# ╔═╡ bd8ed5ad-2038-4a85-bced-37921b57c830
import Infernal

# ╔═╡ e498cf88-adf8-4a40-923a-cfe34602c033
import FASTX

# ╔═╡ 84bd0ac3-ea4a-4700-8dfb-0231405cf89a
import SamApp2024

# ╔═╡ fb4ba584-2243-4b3c-b21c-8c1afdf767c7
PlutoUI.TableOfContents()

# ╔═╡ 26b888bf-7627-40f1-bfa9-1a8f9054380a
md"""
# Load data
"""

# ╔═╡ a285a0bb-a8a8-44d4-9838-730e7393fb1e
natural_seed_stk_147 = Infernal.esl_afetch(Rfam.seed(; rfam_version="14.7"), "RF00162")

# ╔═╡ 78c131e2-73da-4d0d-b037-555f66f50057
natural_seed_afa_147 = Infernal.esl_reformat("afa", natural_seed_stk_147.out)

# ╔═╡ e1e46adb-7cf4-4801-b770-c8ca53e31bdc
natural_seed_ids_147 = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_seed_afa_147.out)))

# ╔═╡ 0f465125-385c-4061-95b9-7c3b8d6fbb28
natural_hits_afa_147 = Infernal.cmalign(
   Infernal.cmfetch(Rfam.cm(; rfam_version="14.7"), "RF00162").out,
   Rfam.fasta_file("RF00162"; rfam_version="14.7");
   outformat="afa"
)

# ╔═╡ 91419950-67dc-4a26-815b-94a778f48148
natural_hits_ids_147 = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_hits_afa_147.out)))

# ╔═╡ b1ff6d74-b5d1-4a8f-b959-b6ea9d8afdbb
natural_hits_dsc_147 = FASTX.description.(FASTX.FASTA.Reader(open(natural_hits_afa_147.out)))

# ╔═╡ 8806bda8-823c-43a1-8436-513f06d9d2b8
natural_hits_sequences_147 = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_hits_afa_147.out)))

# ╔═╡ 980de123-67b7-48de-9acd-7c4427a01395
natural_seed_stk_1410 = Infernal.esl_afetch(Rfam.seed(; rfam_version="14.10"), "RF00162")

# ╔═╡ 0640abdd-9d87-4b51-bcad-1573b7cb6d33
natural_seed_afa_1410 = Infernal.esl_reformat("afa", natural_seed_stk_1410.out)

# ╔═╡ ee666c37-ab51-49ec-b4e1-85b7cc957b78
natural_seed_ids_1410 = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_seed_afa_1410.out)))

# ╔═╡ e5420234-cd76-4d7b-b35e-374d07905881
natural_seed_ids_147 ⊆ natural_seed_ids_1410

# ╔═╡ 5ed6f4b9-6af3-4e07-9ae9-c7ab1985be04
Infernal.cmfetch(Rfam.cm(; rfam_version="14.7"), "RF00162").out

# ╔═╡ f5c34b82-9c1b-4ae6-a550-80d79580e5be
Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out

# ╔═╡ 87caec1b-89ac-4a5a-a6d9-8322faa34374
readlines(Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out)

# ╔═╡ ce6c6d44-875a-4d73-951e-4761578d56c4
natural_hits_afa_1410 = Infernal.cmalign(
   Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out,
   Rfam.fasta_file("RF00162"; rfam_version="14.10");
   outformat="afa"
)

# ╔═╡ fac749a6-245a-48fb-a72f-cccbdf80cc07
natural_hits_ids_1410 = FASTX.identifier.(FASTX.FASTA.Reader(open(natural_hits_afa_1410.out)))

# ╔═╡ 8e4c3995-a6e5-41fb-8c8a-fee3c8990075
natural_hits_dsc_1410 = FASTX.description.(FASTX.FASTA.Reader(open(natural_hits_afa_1410.out)))

# ╔═╡ 5a89855e-6218-461f-9231-c2adc49d2fea
natural_hits_sequences_1410 = FASTX.sequence.(FASTX.FASTA.Reader(open(natural_hits_afa_1410.out)))

# ╔═╡ f394c06c-73fc-448d-bc43-b29857fa30d4
[d for d = natural_hits_dsc_147 if occursin("2GIS", d)], [d for d = natural_hits_dsc_1410 if occursin("2GIS", d)]

# ╔═╡ 857da530-c8a4-4a2c-b855-f9bf0c636d2f
[d for d = natural_hits_dsc_147 if occursin("4KQY", d)], [d for d = natural_hits_dsc_1410 if occursin("4KQY", d)]

# ╔═╡ d4b342ff-9866-4db0-9696-badc5d041a62
[d for d = natural_hits_dsc_147 if occursin("6WLR", d)], [d for d = natural_hits_dsc_1410 if occursin("6WLR", d)]

# ╔═╡ 37eb44b3-7e1a-4e83-b1ea-698189536ee9
[d for d = natural_hits_dsc_147 if occursin("4OQU", d)], [d for d = natural_hits_dsc_1410 if occursin("4OQU", d)]

# ╔═╡ 19331e9c-5e63-4ad3-8d0a-d4e3ec4fb483
SamApp2024.probed_aptamers_table_20221027()

# ╔═╡ Cell order:
# ╠═7bb51d21-270c-4725-a59d-d81bc676070b
# ╠═e1633647-c7e2-46af-af5f-016f247e5d82
# ╠═14a89917-47de-4a78-8d68-c6ee65205a60
# ╠═9262d7dc-c0f2-4711-89ad-570dfbf6e322
# ╠═bd8ed5ad-2038-4a85-bced-37921b57c830
# ╠═e498cf88-adf8-4a40-923a-cfe34602c033
# ╠═84bd0ac3-ea4a-4700-8dfb-0231405cf89a
# ╠═fb4ba584-2243-4b3c-b21c-8c1afdf767c7
# ╠═26b888bf-7627-40f1-bfa9-1a8f9054380a
# ╠═a285a0bb-a8a8-44d4-9838-730e7393fb1e
# ╠═78c131e2-73da-4d0d-b037-555f66f50057
# ╠═e1e46adb-7cf4-4801-b770-c8ca53e31bdc
# ╠═0f465125-385c-4061-95b9-7c3b8d6fbb28
# ╠═91419950-67dc-4a26-815b-94a778f48148
# ╠═b1ff6d74-b5d1-4a8f-b959-b6ea9d8afdbb
# ╠═8806bda8-823c-43a1-8436-513f06d9d2b8
# ╠═980de123-67b7-48de-9acd-7c4427a01395
# ╠═0640abdd-9d87-4b51-bcad-1573b7cb6d33
# ╠═ee666c37-ab51-49ec-b4e1-85b7cc957b78
# ╠═e5420234-cd76-4d7b-b35e-374d07905881
# ╠═5ed6f4b9-6af3-4e07-9ae9-c7ab1985be04
# ╠═f5c34b82-9c1b-4ae6-a550-80d79580e5be
# ╠═87caec1b-89ac-4a5a-a6d9-8322faa34374
# ╠═ce6c6d44-875a-4d73-951e-4761578d56c4
# ╠═fac749a6-245a-48fb-a72f-cccbdf80cc07
# ╠═8e4c3995-a6e5-41fb-8c8a-fee3c8990075
# ╠═5a89855e-6218-461f-9231-c2adc49d2fea
# ╠═f394c06c-73fc-448d-bc43-b29857fa30d4
# ╠═857da530-c8a4-4a2c-b855-f9bf0c636d2f
# ╠═d4b342ff-9866-4db0-9696-badc5d041a62
# ╠═37eb44b3-7e1a-4e83-b1ea-698189536ee9
# ╠═19331e9c-5e63-4ad3-8d0a-d4e3ec4fb483

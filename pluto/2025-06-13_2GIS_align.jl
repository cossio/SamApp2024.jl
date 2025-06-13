### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ a10070de-4884-11f0-146d-8fab10b1f04e
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 38cd997c-1a61-426a-b5b4-545f65461217
md"# Imports"

# ╔═╡ b934cb60-4920-4596-9d7c-ecb468924232
import Makie, CairoMakie, PlutoUI, FASTX, Infernal, Rfam, ViennaRNA, BioSequences

# ╔═╡ cb6245f5-8c50-47ee-9fc0-ac5fa6622c56
# CM model from Rfam (this has the noisy floor!)
Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162");

# ╔═╡ 4bd65dc4-4fa2-418d-94ab-5aff80faf38e
Rfam.fasta_file("RF00162")

# ╔═╡ f04ed226-3ad4-47d5-82a9-00d50502c0ac


# ╔═╡ 73e4226f-95a5-4f65-943f-a2ba95436a7a
#RF00162_hits_afa = Infernal.cmalign(Rfam_cm.out, Rfam.fasta_file("RF00162"); matchonly=true, outformat="AFA");

# ╔═╡ Cell order:
# ╠═38cd997c-1a61-426a-b5b4-545f65461217
# ╠═a10070de-4884-11f0-146d-8fab10b1f04e
# ╠═b934cb60-4920-4596-9d7c-ecb468924232
# ╠═cb6245f5-8c50-47ee-9fc0-ac5fa6622c56
# ╠═4bd65dc4-4fa2-418d-94ab-5aff80faf38e
# ╠═f04ed226-3ad4-47d5-82a9-00d50502c0ac
# ╠═73e4226f-95a5-4f65-943f-a2ba95436a7a

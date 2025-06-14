### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# ╔═╡ a10070de-4884-11f0-146d-8fab10b1f04e
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 47bf2463-6930-4eb9-869a-ebbf9b73c3c3
using DataFrames: DataFrame

# ╔═╡ 38cd997c-1a61-426a-b5b4-545f65461217
md"# Imports"

# ╔═╡ b934cb60-4920-4596-9d7c-ecb468924232
import Makie, CairoMakie, PlutoUI, FASTX, Infernal, Rfam, ViennaRNA, BioSequences, BioAlignments, SamApp2024, Printf

# ╔═╡ 4eda5c48-9843-482b-8293-aefed4072805
pdb_2GIS_raw_seq = SamApp2024.artifact_load_pdb_2GIS_sequence()

# ╔═╡ 4eca60b5-3b85-4108-873f-da9158a44ea0
pdb_4KQY_raw_seq = SamApp2024.artifact_load_pdb_4KQY_sequence()

# ╔═╡ 058c7e1b-a56a-4d50-ab87-c4c6224b4d21
pdb_2GIS_aligned_fasta = Infernal.cmalign(
	Infernal.cmfetch(Rfam.cm(; rfam_version="14.0"), "RF00162").out,
	SamApp2024.artifact_path_pdb_2GIS_sequence();
	matchonly=false, outformat="afa"
 )

# ╔═╡ 99661603-7b29-4a33-9a66-3c21c432f324
pdb_4KQY_aligned_fasta = Infernal.cmalign(
	Infernal.cmfetch(Rfam.cm(; rfam_version="14.0"), "RF00162").out,
	SamApp2024.artifact_path_pdb_4KQY_sequence();
	matchonly=false, outformat="afa"
 )

# ╔═╡ 1f717749-9862-4ebd-a64c-fc222c07a3a4
pdb_2GIS_aln_seq_with_inserts = FASTX.FASTA.Reader(open(pdb_2GIS_aligned_fasta.out)) do reader
	only(FASTX.sequence.(reader))
end

# ╔═╡ 839ffaa0-1353-4104-958e-ff7b9a4a3e75
pdb_4KQY_aln_seq_with_inserts = FASTX.FASTA.Reader(open(pdb_4KQY_aligned_fasta.out)) do reader
	only(FASTX.sequence.(reader))
end

# ╔═╡ fb70097e-d7a7-4c51-947d-2b5f6d0f762d
pdb_2GIS_ref_seq = join(a == '-' ? 'N' : islowercase(a) ? '-' : a for a = pdb_2GIS_aln_seq_with_inserts)

# ╔═╡ 840749ac-6c91-4be5-9da1-a929c94aa824
pdb_4KQY_ref_seq = join(a == '-' ? 'N' : islowercase(a) ? '-' : a for a = pdb_4KQY_aln_seq_with_inserts)

# ╔═╡ bc4961dd-56b8-41e3-8da5-49d47194be93
bioaln_2GIS = BioAlignments.AlignedSequence(BioSequences.LongRNA{4}(pdb_2GIS_aln_seq_with_inserts), BioSequences.LongRNA{4}(pdb_2GIS_ref_seq))

# ╔═╡ f4bb0db9-af41-49e0-833f-6119f6af1aae
bioaln_4KQY = BioAlignments.AlignedSequence(BioSequences.LongRNA{4}(pdb_4KQY_aln_seq_with_inserts), BioSequences.LongRNA{4}(pdb_4KQY_ref_seq))

# ╔═╡ 0de74586-9147-4e0f-a298-80510621defe
function align_map(aln_seq_with_inserts::AbstractString)
	positions_in_aligned = Int[]
	anchor = 0
	for (i, letter) = enumerate(aln_seq_with_inserts)
		if isuppercase(letter)
			anchor += 1	
			push!(positions_in_aligned, anchor)
		elseif letter == '-'
			anchor += 1
		elseif islowercase(letter)
			push!(positions_in_aligned, anchor)
		else
			error("Unexpected")
		end
	end
	return positions_in_aligned
end

# ╔═╡ 175b595c-3cd9-46de-bbfb-d7454e0a30a1
align_map(pdb_2GIS_aln_seq_with_inserts)[32]

# ╔═╡ e3b02d6a-c1f0-4fde-aa46-01643257fb72
align_map(pdb_4KQY_aln_seq_with_inserts)[(7) .+ 1]

# ╔═╡ Cell order:
# ╠═38cd997c-1a61-426a-b5b4-545f65461217
# ╠═a10070de-4884-11f0-146d-8fab10b1f04e
# ╠═b934cb60-4920-4596-9d7c-ecb468924232
# ╠═47bf2463-6930-4eb9-869a-ebbf9b73c3c3
# ╠═4eda5c48-9843-482b-8293-aefed4072805
# ╠═4eca60b5-3b85-4108-873f-da9158a44ea0
# ╠═058c7e1b-a56a-4d50-ab87-c4c6224b4d21
# ╠═99661603-7b29-4a33-9a66-3c21c432f324
# ╠═1f717749-9862-4ebd-a64c-fc222c07a3a4
# ╠═839ffaa0-1353-4104-958e-ff7b9a4a3e75
# ╠═fb70097e-d7a7-4c51-947d-2b5f6d0f762d
# ╠═840749ac-6c91-4be5-9da1-a929c94aa824
# ╠═bc4961dd-56b8-41e3-8da5-49d47194be93
# ╠═f4bb0db9-af41-49e0-833f-6119f6af1aae
# ╠═0de74586-9147-4e0f-a298-80510621defe
# ╠═175b595c-3cd9-46de-bbfb-d7454e0a30a1
# ╠═e3b02d6a-c1f0-4fde-aa46-01643257fb72

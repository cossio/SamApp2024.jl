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

# ╔═╡ 73e4226f-95a5-4f65-943f-a2ba95436a7a
raw_seq = FASTX.FASTA.Reader(open(SamApp2024.artifact_path_pdb_2GIS_sequence())) do reader
	only(FASTX.sequence.(reader))
end

# ╔═╡ 058c7e1b-a56a-4d50-ab87-c4c6224b4d21
pdb_2gis_aligned = Infernal.cmalign(
	Infernal.cmfetch(Rfam.cm(; rfam_version="14.0"), "RF00162").out,
	SamApp2024.artifact_path_pdb_2GIS_sequence();
	matchonly=false, outformat="afa"
 )

# ╔═╡ cc8a161a-a25f-49e5-b81b-9c145b0b3a39
aln_seq_with_inserts = FASTX.FASTA.Reader(open(pdb_2gis_aligned.out)) do reader
	only(FASTX.sequence.(reader))
end

# ╔═╡ fb70097e-d7a7-4c51-947d-2b5f6d0f762d
ref_seq = join(a == '-' ? 'N' : islowercase(a) ? '-' : a for a = aln_seq_with_inserts)

# ╔═╡ c14e2c85-551d-4446-8f53-8f06b6fdbc72
length(ref_seq)

# ╔═╡ fa970632-fe01-4a7c-9d3e-ef94fe4ee570
length(filter(!=('-'), aln_seq_with_inserts)) == length(raw_seq)

# ╔═╡ d5728266-3ae7-424d-8478-d6b7b627bdbc
aln_seq_no_inserts = filter(!islowercase, aln_seq_with_inserts)

# ╔═╡ bc4961dd-56b8-41e3-8da5-49d47194be93
bioaln = BioAlignments.AlignedSequence(BioSequences.LongRNA{4}(aln_seq_with_inserts), BioSequences.LongRNA{4}(ref_seq))

# ╔═╡ 026307e6-a873-4755-a54d-bc9dba03999c
bioaln.aln[1]

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

# ╔═╡ 4c090621-2748-41f3-9086-c2108efccb47
align_map(aln_seq_with_inserts) |> length

# ╔═╡ a038a262-6418-4718-b402-150fa99e9b2e
length(raw_seq)

# ╔═╡ cf02edad-3739-480f-8c17-908c4887e61c
row1 = map(string, align_map(aln_seq_with_inserts))

# ╔═╡ 4f31fbc5-1d86-4775-bb60-9758f7809ddb
row2 = raw_seq

# ╔═╡ 1b4639cb-8a91-4602-90b6-4a92a0612720
widths = map((a,b) -> max(length(a), length(b)), row1, row2)

# ╔═╡ 6996e675-6a2f-45e4-832e-d49d4b3b1d47
fmt = join(["%-$(w)s" for w in widths], "  ") * "\n"  # e.g. "%-6s  %-6s  %-8s\n"

# ╔═╡ e8f2d822-b8fc-4529-9334-550d208d5d0e
for row = (row1, row2)
	for (width, word) = zip(widths, collect(row))
		print(rpad(word, width + 1))
	end
	println()
end

# ╔═╡ 175b595c-3cd9-46de-bbfb-d7454e0a30a1
align_map(aln_seq_with_inserts)[7]

# ╔═╡ 7fe859db-3853-4a93-82cd-ddb7e497fa52
cumsum(map(islowercase, collect(aln_seq_with_inserts)))

# ╔═╡ 1188b2b5-95d9-484f-9f83-3b612b679052
cumsum(map(==('-'), collect(aln_seq_with_inserts)))

# ╔═╡ Cell order:
# ╠═38cd997c-1a61-426a-b5b4-545f65461217
# ╠═a10070de-4884-11f0-146d-8fab10b1f04e
# ╠═b934cb60-4920-4596-9d7c-ecb468924232
# ╠═47bf2463-6930-4eb9-869a-ebbf9b73c3c3
# ╠═73e4226f-95a5-4f65-943f-a2ba95436a7a
# ╠═058c7e1b-a56a-4d50-ab87-c4c6224b4d21
# ╠═cc8a161a-a25f-49e5-b81b-9c145b0b3a39
# ╠═fb70097e-d7a7-4c51-947d-2b5f6d0f762d
# ╠═c14e2c85-551d-4446-8f53-8f06b6fdbc72
# ╠═fa970632-fe01-4a7c-9d3e-ef94fe4ee570
# ╠═d5728266-3ae7-424d-8478-d6b7b627bdbc
# ╠═bc4961dd-56b8-41e3-8da5-49d47194be93
# ╠═026307e6-a873-4755-a54d-bc9dba03999c
# ╠═0de74586-9147-4e0f-a298-80510621defe
# ╠═4c090621-2748-41f3-9086-c2108efccb47
# ╠═a038a262-6418-4718-b402-150fa99e9b2e
# ╠═cf02edad-3739-480f-8c17-908c4887e61c
# ╠═4f31fbc5-1d86-4775-bb60-9758f7809ddb
# ╠═1b4639cb-8a91-4602-90b6-4a92a0612720
# ╠═6996e675-6a2f-45e4-832e-d49d4b3b1d47
# ╠═e8f2d822-b8fc-4529-9334-550d208d5d0e
# ╠═175b595c-3cd9-46de-bbfb-d7454e0a30a1
# ╠═7fe859db-3853-4a93-82cd-ddb7e497fa52
# ╠═1188b2b5-95d9-484f-9f83-3b612b679052

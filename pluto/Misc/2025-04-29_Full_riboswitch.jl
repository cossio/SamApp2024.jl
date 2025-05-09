### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 9e5ab2f6-2505-11f0-2a08-2599aa3bd022
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 75226411-330d-4450-a6d2-0ca3d6774f65
using Statistics: mean, cov

# ╔═╡ b97d20b6-70e2-4194-b49c-3166761d6403
using LogExpFunctions: xlogx

# ╔═╡ 73398a3e-1c08-4416-ac15-9bbe282834ce
using LinearAlgebra: diagind

# ╔═╡ 48c28fd4-1b35-4632-9215-a24f3ecb1a80
using Unitful: ustrip

# ╔═╡ c89c23fe-3b38-415b-b11d-541e814d94fb
md"# Imports"

# ╔═╡ 2d70c8b7-8fed-47d4-9b88-c90dbb0b4cbc
import SamApp2024

# ╔═╡ 68544f95-f782-4560-b4c3-0aa7f643ef75
import BioSequences

# ╔═╡ 81a424bf-f09b-480a-99b3-528a561ecfad
import Makie, CairoMakie

# ╔═╡ f3ec55d7-80e5-467d-8aa6-7f7804dd5f58
import ViennaRNA

# ╔═╡ fef739cf-e49c-4a34-b987-85b90b356bc9
import Logomaker

# ╔═╡ ccd05903-4194-4d11-9025-1aac33df4fde
md"# Data"

# ╔═╡ 4083a737-315c-467c-9b0a-eb6e829e47dd
fullseqs = SamApp2024.artifact_load_20250428_aligned_full_riboswitch_sequences()

# ╔═╡ 6bd76c54-45cc-4158-b6c0-604503a2fb28
alnseqs = map(BioSequences.LongRNA{4}, [filter(!islowercase, replace(seq, '.' => "")) for seq = fullseqs])

# ╔═╡ 7f6c9a0b-7bf4-4e07-8226-38de3746420a
X = SamApp2024.onehot(map(BioSequences.LongRNA{4}, alnseqs));

# ╔═╡ 71a5a43e-bd92-40c3-8243-ed5c9532c64e
let fig = Makie.Figure()
	C = dropdims(sum(abs2, reshape(cov(reshape(X, :, size(X, 3)); dims=2), 5, 153, 5, 153); dims=(1,3)); dims=(1,3))
	C[diagind(C)] .= 0
	
	ax = Makie.Axis(fig[1,1]; width=700, height=700, xticks=0:10:150, yticks=0:10:150)
	Makie.heatmap!(ax, C)
	Makie.vlines!(ax, 108; color=:white, linestyle=:dash)
	Makie.hlines!(ax, 108; color=:white, linestyle=:dash)
	Makie.vlines!(ax, 116; color=:red, linestyle=:dash)
	Makie.vlines!(ax, 123; color=:red, linestyle=:dash)
	Makie.vlines!(ax, 131; color=:red, linestyle=:dash)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 56f9abc6-f0f4-4850-8dc5-f1df2158f32c
alnseqs[1][1:8]

# ╔═╡ 6715935d-e456-41c5-b985-c1d80f6bc583
alnseqs[1][116:123]

# ╔═╡ e1003c93-585e-43cf-add9-0daa4d2069d2
[occursin(string(s[1:8]), string(s[109:end])) for s = alnseqs]

# ╔═╡ fa6a0fc5-8e1f-4622-9952-17989b7d2d1a
[occursin(string(s[101:108]), string(s[109:end])) for s = alnseqs]

# ╔═╡ a0afc36c-fda8-4a6e-a284-0cf3553a74d3
seqlogo_color_scheme = Logomaker.color_scheme(
	'C' => "blue",
	'U' => "red",
	'A' => "green",
	'G' => "orange",
	'⊟' => "black"
)

# ╔═╡ 8262229d-1492-47f0-a6bb-ff024e8fc4f8
xlog2x(x) = xlogx(x) / log(oftype(x,2))

# ╔═╡ 7abd350a-2877-42a4-994e-c5ce62ac7f47
function seqlogo_entropic(p::AbstractMatrix; max_ylim=true)
    @assert size(p, 1) == 5 # nucleotides + gap
    w = p ./ sum(p; dims=1)
    H = sum(-xlog2x.(w); dims=1)
    @assert all(0 .≤ H .≤ log2(5))

    cons = w .* (log2(5) .- H)
    logo = Logomaker.Logo(cons, collect("ACGU⊟"); color_scheme=seqlogo_color_scheme)
    max_ylim && logo.ax.set_ylim(0, log2(5))
    logo.ax.set_ylabel("conservation (bits)")
    logo.ax.set_xlabel("site")

    return logo
end

# ╔═╡ daa6be95-b968-4002-8ad4-b1c8aa2d38b7
seqlogo_entropic(reshape(mean(X; dims=3), 5, :)).fig

# ╔═╡ fafcc597-a322-48a2-ad6d-eb9ef95b6527
wuss = SamApp2024.rfam_ss("RF00162"; inserts=false)

# ╔═╡ 84cb87f0-e81a-4799-aa19-57b9b5fb63c3
ss = SamApp2024.clean_wuss(wuss)

# ╔═╡ f1de1a29-8301-447f-adad-ec906d94724a
ss_without_P1 = join([i ∈ SamApp2024.RF00162_sites_annotated_secondary_structure().p1 ? '.' : c for (i,c) in enumerate(ss)]);

# ╔═╡ a0ee4304-50ff-440a-b363-2ea91a911384
# full seq, with P1 and terminator bound (OFF state)
ss_full_P1_and_terminator = ss * repeat('.', length(109:115)) * repeat('(', length(116:123)) * repeat('.', length(124:130)) * repeat(')', length(131:138)) * repeat('.', length(139:153))

# ╔═╡ 85895b80-369c-474e-a134-9b95b45b60a6
# full seq, with P1 but without terminator
ss_full_P1_and_no_terminator = ss * repeat('.', length(109:153))

# ╔═╡ d1244e1c-91f8-4378-8734-c52d3fa95385
# full seq, with P1 unbound and anti-terminator bound (ON state)
ss_full_noP1 = repeat('.', 8) * ss[9:100] * repeat('(', 8) *
	repeat('.', length(109:115)) * repeat(')', length(116:123)) * repeat('.', length(124:130)) * repeat('.', length(131:138)) * repeat('.', length(139:153))

# ╔═╡ a6744f6a-477d-47c4-9dcf-6466026893a5
Vienna_energies_aptamer_only_ss = [ustrip(ViennaRNA.energy(string(seq[1:108]), ss)) for seq = alnseqs];

# ╔═╡ 93c341ff-dde3-4dd9-b245-fc9db2d9a5ab
Vienna_energies_aptamer_only_noP1 = [ustrip(ViennaRNA.energy(string(seq[1:108]), ss_without_P1)) for seq = alnseqs];

# ╔═╡ 5eefdb8e-85bd-4116-9346-e9394de89475
Vienna_energies_full_P1_and_terminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_P1_and_terminator)) for seq = alnseqs];

# ╔═╡ 06236adb-22da-4e8e-80f6-f4398e1f5b5f
Vienna_energies_full_P1_and_no_terminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_P1_and_no_terminator)) for seq = alnseqs];

# ╔═╡ e66a8b94-73a8-42a4-b65f-4190d37fe172
Vienna_energies_full_noP1_antiterminator = [ustrip(ViennaRNA.energy(string(seq), ss_full_noP1)) for seq = alnseqs];

# ╔═╡ c5735f3f-8982-4eab-bc30-39b0b80a27c5
P1_deltaF_aptamer_only = Vienna_energies_aptamer_only_ss - Vienna_energies_aptamer_only_noP1

# ╔═╡ 4a8133b7-f404-4364-ae3a-ab7457364456
P1_deltaF_full = Vienna_energies_full_P1_and_terminator - Vienna_energies_full_noP1_antiterminator

# ╔═╡ eddb9d24-017e-40a6-ae17-00cc4e407e00
P1_deltaF_full_no_terminator = Vienna_energies_full_P1_and_no_terminator - Vienna_energies_full_noP1_antiterminator

# ╔═╡ ee17763e-7ab0-4234-b10f-cca58f12e814
let fig = Makie.Figure()
	_sz = 350
	
	ax = Makie.Axis(fig[1,1]; width=_sz, height=_sz, xlabel="ΔF(OFF) - ΔF(ON) for aptamer only", ylabel="ΔF(OFF) - ΔF(ON) for full riboswitch")
	Makie.scatter!(ax, P1_deltaF_aptamer_only, P1_deltaF_full)
	Makie.ablines!(ax, 0, 1; color=:red, linestyle=:dash)

	ax = Makie.Axis(fig[1,2]; width=_sz, height=_sz, xlabel="ΔF(OFF) - ΔF(ON) for aptamer only", ylabel="ΔF(OFF) - ΔF(ON) for full riboswitch (without terminator)")
	Makie.scatter!(ax, P1_deltaF_aptamer_only, P1_deltaF_full_no_terminator)
	Makie.ablines!(ax, 0, 1; color=:red, linestyle=:dash)
	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ b2b18a43-2a76-4f34-b041-0d123ea48e10
Makie.heatmap(ViennaRNA.bpp(string(alnseqs[1][1:108]))[1:8,101:108])

# ╔═╡ 0fbc4fa8-5903-4afb-9759-2ac24a637003
BPs_full = [ViennaRNA.bpp(string(seq)) for seq = alnseqs]

# ╔═╡ 7648e808-a0e0-485f-8bd8-5d898b5b3941
BPs_aptamer = [ViennaRNA.bpp(string(seq[1:108])) for seq = alnseqs]

# ╔═╡ 43564b0d-3a29-4f88-8d15-8120300660fa
let seq = string(alnseqs[30])
	println(seq[1:108])
	println(ViennaRNA.mfe(seq[1:108])[1])
	println(seq)
	println(ViennaRNA.mfe(seq)[1])
end

# ╔═╡ 888345b7-8661-47d5-bedd-5d59612071b6
let fig = Makie.Figure()
	bins = 0:0.1:1
	for i = 1:8
		ax = Makie.Axis(fig[fld1(i, 4), mod1(i, 4)]; width=200, height=200, xlabel="Base-pair prob. P1", ylabel="aptamers (freq.)", title="site = $i")
		
		Makie.hist!(ax, [bp[i, 109 - i] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i, 109 - i] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)

		Makie.vlines!(ax, mean([bp[i, 109 - i] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i, 109 - i] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)

		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 10)
		
		if i == 8
			Makie.axislegend(ax; position=:lt)
		end
	end
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ abfb5f30-5e6d-4c3f-8030-81b319787203
prod(mean([bp[i, 109 - i] for bp = BPs_full]) for i=1:8) / 
	prod(mean([bp[i, 109 - i] for bp = BPs_aptamer]) for i=1:8)

# ╔═╡ 858cf27e-ea29-4d4f-9205-c6de4ba49f33
let fig = Makie.Figure()
	bins = 0:0.1:1
	for (n, (i,j)) = enumerate(zip([13:17; 21:23], reverse([29:31; 38:42])))
		ax = Makie.Axis(fig[fld1(n, 4), mod1(n, 4)]; width=200, height=200, xlabel="Base-pair prob. (P2)", ylabel="aptamers (freq.)", title="site = $n")
		
		Makie.hist!(ax, [bp[i,j] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i,j] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)

		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)

		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 10)
		
		if n == 8
			Makie.axislegend(ax; position=:lt)
		end
	end
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ c783cf11-f5c2-4ef5-8dc4-cf9dc0e4b479
prod(mean([bp[i,j] for bp = BPs_full]) for (i,j) = zip([13:17; 21:23], reverse([29:31; 38:42]))) / 
	prod(mean([bp[i,j] for bp = BPs_aptamer]) for (i,j) = zip([13:17; 21:23], reverse([29:31; 38:42])))

# ╔═╡ df7fccd6-a5c4-47c4-afe4-a55a19ed57e3
let fig = Makie.Figure()
	bins = 0:0.1:1
	for i = 81:86
		j = 97 - i + 81
		n = i - 80
		ax = Makie.Axis(fig[fld1(n, 4), mod1(n,4)]; width=200, height=200, xlabel="Base-pair prob. (P4)", ylabel="aptamers (freq.)", title="sites = ($i,$j)")
		
		Makie.hist!(ax, [bp[i,j] for bp = BPs_full]; normalization=:pdf, bins, label="full", color=:gray)
		Makie.stephist!(ax, [bp[i,j] for bp = BPs_aptamer]; normalization=:pdf, bins, label="aptamer only", color=:blue, linewidth=2)

		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_full]); color=:black, linestyle=:dash, linewidth=4)
		Makie.vlines!(ax, mean([bp[i,j] for bp = BPs_aptamer]); color=:blue, linestyle=:dash, linewidth=4)

		Makie.xlims!(ax, 0, 1)
		Makie.ylims!(ax, 0, 10)
		
		if n == 4
			Makie.axislegend(ax; position=:lt)
		end
	end
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 64128cc4-44dc-4023-9600-1dcca76cfe31
prod(mean([bp[i,97-i+81] for bp = BPs_full]) for i=81:86) / prod(mean([bp[i,97-i+81] for bp = BPs_aptamer]) for i=81:86)

# ╔═╡ e4044118-c03d-4b20-ba8c-0284d7d250fb
SamApp2024.RF00162_sites_annotated_secondary_structure().p4

# ╔═╡ 21c629a1-5b4e-45ee-ae6b-18da275f3d0e
length([13:16; 21:23])

# ╔═╡ Cell order:
# ╠═c89c23fe-3b38-415b-b11d-541e814d94fb
# ╠═9e5ab2f6-2505-11f0-2a08-2599aa3bd022
# ╠═2d70c8b7-8fed-47d4-9b88-c90dbb0b4cbc
# ╠═68544f95-f782-4560-b4c3-0aa7f643ef75
# ╠═81a424bf-f09b-480a-99b3-528a561ecfad
# ╠═f3ec55d7-80e5-467d-8aa6-7f7804dd5f58
# ╠═fef739cf-e49c-4a34-b987-85b90b356bc9
# ╠═75226411-330d-4450-a6d2-0ca3d6774f65
# ╠═b97d20b6-70e2-4194-b49c-3166761d6403
# ╠═73398a3e-1c08-4416-ac15-9bbe282834ce
# ╠═48c28fd4-1b35-4632-9215-a24f3ecb1a80
# ╠═ccd05903-4194-4d11-9025-1aac33df4fde
# ╠═4083a737-315c-467c-9b0a-eb6e829e47dd
# ╠═6bd76c54-45cc-4158-b6c0-604503a2fb28
# ╠═7f6c9a0b-7bf4-4e07-8226-38de3746420a
# ╠═71a5a43e-bd92-40c3-8243-ed5c9532c64e
# ╠═56f9abc6-f0f4-4850-8dc5-f1df2158f32c
# ╠═6715935d-e456-41c5-b985-c1d80f6bc583
# ╠═e1003c93-585e-43cf-add9-0daa4d2069d2
# ╠═fa6a0fc5-8e1f-4622-9952-17989b7d2d1a
# ╠═a0afc36c-fda8-4a6e-a284-0cf3553a74d3
# ╠═8262229d-1492-47f0-a6bb-ff024e8fc4f8
# ╠═7abd350a-2877-42a4-994e-c5ce62ac7f47
# ╠═daa6be95-b968-4002-8ad4-b1c8aa2d38b7
# ╠═fafcc597-a322-48a2-ad6d-eb9ef95b6527
# ╠═84cb87f0-e81a-4799-aa19-57b9b5fb63c3
# ╠═f1de1a29-8301-447f-adad-ec906d94724a
# ╠═a0ee4304-50ff-440a-b363-2ea91a911384
# ╠═85895b80-369c-474e-a134-9b95b45b60a6
# ╠═d1244e1c-91f8-4378-8734-c52d3fa95385
# ╠═a6744f6a-477d-47c4-9dcf-6466026893a5
# ╠═93c341ff-dde3-4dd9-b245-fc9db2d9a5ab
# ╠═5eefdb8e-85bd-4116-9346-e9394de89475
# ╠═06236adb-22da-4e8e-80f6-f4398e1f5b5f
# ╠═e66a8b94-73a8-42a4-b65f-4190d37fe172
# ╠═c5735f3f-8982-4eab-bc30-39b0b80a27c5
# ╠═4a8133b7-f404-4364-ae3a-ab7457364456
# ╠═eddb9d24-017e-40a6-ae17-00cc4e407e00
# ╠═ee17763e-7ab0-4234-b10f-cca58f12e814
# ╠═b2b18a43-2a76-4f34-b041-0d123ea48e10
# ╠═0fbc4fa8-5903-4afb-9759-2ac24a637003
# ╠═7648e808-a0e0-485f-8bd8-5d898b5b3941
# ╠═43564b0d-3a29-4f88-8d15-8120300660fa
# ╠═888345b7-8661-47d5-bedd-5d59612071b6
# ╠═abfb5f30-5e6d-4c3f-8030-81b319787203
# ╠═858cf27e-ea29-4d4f-9205-c6de4ba49f33
# ╠═c783cf11-f5c2-4ef5-8dc4-cf9dc0e4b479
# ╠═df7fccd6-a5c4-47c4-afe4-a55a19ed57e3
# ╠═64128cc4-44dc-4023-9600-1dcca76cfe31
# ╠═e4044118-c03d-4b20-ba8c-0284d7d250fb
# ╠═21c629a1-5b4e-45ee-ae6b-18da275f3d0e

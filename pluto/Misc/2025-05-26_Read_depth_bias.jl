### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ bdb8966d-ee65-4cc1-8e03-95edc7248551
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ 5e8b92b5-66f0-4f28-8010-c103de66d1b4
using Distributions: Gamma

# ╔═╡ b4804a68-60fd-4655-a9a6-b0b82f4f70c4
using Statistics: mean, std

# ╔═╡ 468738a8-3a61-11f0-0206-291073c6db4d
md"# Imports"

# ╔═╡ 658362d5-eb92-473d-9f2b-ea83429e9363
import Makie, CairoMakie

# ╔═╡ 002d9917-53e7-462f-8b81-553e859c9d73
import SamApp2024

# ╔═╡ dc281c13-24ae-481a-918d-fb69e3392435
# load SHAPE data
shape_data_045 = SamApp2024.load_shapemapper_data_pierre_demux_20230920(; demux=true);

# ╔═╡ a41e0a5b-f029-4b1f-a386-570d5603391b
shape_data = SamApp2024.select_conditions_20231002(shape_data_045, filter(endswith("_rep0"), shape_data_045.conditions));

# ╔═╡ 2b5f7847-bec2-4cd1-b95d-d04dc5bec30e
shape_gamma_theta_M = (shape_data.shape_M_stderr.^2) ./ shape_data.shape_M;

# ╔═╡ 3899ac6f-e9ba-426a-a916-17eb40b54847
shape_gamma_theta_U = (shape_data.shape_U_stderr.^2) ./ shape_data.shape_U;

# ╔═╡ 8ac4b744-7c1c-466c-92a7-72bd315d7b83
shape_gamma_theta_D = (shape_data.shape_D_stderr.^2) ./ shape_data.shape_D;

# ╔═╡ f1bd692c-fc59-4c3c-8f15-021a25201e85
shape_gamma_alpha_M = (shape_data.shape_M ./ shape_data.shape_M_stderr).^2;

# ╔═╡ f48e3d7a-24fc-4ec0-a031-043b256f5541
shape_gamma_alpha_U = (shape_data.shape_U ./ shape_data.shape_U_stderr).^2;

# ╔═╡ 24957a09-2b7c-46cf-8e8c-dff49186fd4d
shape_gamma_alpha_D = (shape_data.shape_D ./ shape_data.shape_D_stderr).^2;

# ╔═╡ 562a46fb-d29e-43b2-99c3-ca5e0691765b
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=500, height=500)
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_M_depth)); normalization=:pdf, label="M")
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_U_depth)); normalization=:pdf, label="U")
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_D_depth)); normalization=:pdf, label="D")
	Makie.axislegend(ax; position=:rt)

	ax = Makie.Axis(fig[1,2]; width=500, height=500)
	bins = 0:0.001:0.1
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_M)); normalization=:pdf, bins)
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_U)); normalization=:pdf, bins)
	Makie.stephist!(ax, filter(!isnan, vec(shape_data.shape_D)); normalization=:pdf, bins)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ c7a6a735-8c5b-497c-91eb-8016b50fd8d2
let fig = Makie.Figure()
	_sz = 400
	
	pos_M = findall(!isnan, vec(shape_data.shape_M_depth)) ∩ findall(!isnan, vec(shape_data.shape_M))
	ax = Makie.Axis(fig[1,1]; width=_sz, height=_sz, xlabel="M_depth", ylabel="M")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_M_depth)[pos_M], vec(shape_data.shape_M)[pos_M]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 0.03)

	pos_U = findall(!isnan, vec(shape_data.shape_U_depth)) ∩ findall(!isnan, vec(shape_data.shape_U))
	ax = Makie.Axis(fig[1,2]; width=_sz, height=_sz, xlabel="U_depth", ylabel="U")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_U_depth)[pos_U], vec(shape_data.shape_U)[pos_U]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 0.03)

	pos_D = findall(!isnan, vec(shape_data.shape_D_depth)) ∩ findall(!isnan, vec(shape_data.shape_D))
	ax = Makie.Axis(fig[1,3]; width=_sz, height=_sz, xlabel="D_depth", ylabel="D")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_D_depth)[pos_D], vec(shape_data.shape_D)[pos_D]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 0.03)

	pos = findall(!isnan, vec(shape_data.shape_M_depth)) ∩ findall(!isnan, vec(shape_data.shape_U_depth))
	ax = Makie.Axis(fig[2,1]; width=_sz, height=_sz, xlabel="M_depth", ylabel="U_depth")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_M_depth)[pos], vec(shape_data.shape_U_depth)[pos]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 6e4)

	pos = findall(!isnan, vec(shape_data.shape_M_depth)) ∩ findall(!isnan, vec(shape_data.shape_D_depth))
	ax = Makie.Axis(fig[2,2]; width=_sz, height=_sz, xlabel="M_depth", ylabel="D_depth")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_M_depth)[pos], vec(shape_data.shape_D_depth)[pos]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 6e4)

	pos = findall(!isnan, vec(shape_data.shape_U_depth)) ∩ findall(!isnan, vec(shape_data.shape_D_depth))
	ax = Makie.Axis(fig[2,3]; width=_sz, height=_sz, xlabel="U_depth", ylabel="D_depth")
	Makie.datashader!(ax, Makie.Point2f.(vec(shape_data.shape_U_depth)[pos], vec(shape_data.shape_D_depth)[pos]))
	Makie.xlims!(ax, 0, 6e4)
	Makie.ylims!(ax, 0, 6e4)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ a97ff5ca-c6b7-4bb5-b5c4-09aafb85a8e3
L, N, C = size(shape_data.shape_reactivities)

# ╔═╡ 61417d18-0bd1-407f-8e33-fc8d1c3796c7
nsamples = 100

# ╔═╡ 98e0c738-904b-4778-aa0b-c4a294b80d77
begin
	R_list = Float32[]
	M_depth_list = Float32[]
	U_depth_list = Float32[]
	D_depth_list = Float32[]
	
	for c = 1:C, n = 1:N, i = 1:L
		α_M = shape_gamma_alpha_M[i,n,c]
		α_U = shape_gamma_alpha_U[i,n,c]
		α_D = shape_gamma_alpha_D[i,n,c]
	
		θ_M = shape_gamma_theta_M[i,n,c]
		θ_U = shape_gamma_theta_U[i,n,c]
		θ_D = shape_gamma_theta_D[i,n,c]
		
		if isnan(α_M) || isnan(α_U) || isnan(α_D) || isnan(θ_M) || isnan(θ_U) || isnan(θ_D)
			continue
		end
	
		M = rand(Gamma(α_M, θ_M), nsamples)
		U = rand(Gamma(α_U, θ_U), nsamples)
		D = rand(Gamma(α_D, θ_D), nsamples)
	
		R = (M - U) ./ D
	
		append!(R_list, R)
		append!(M_depth_list, fill(shape_data.shape_M_depth[i,n,c], nsamples))
		append!(U_depth_list, fill(shape_data.shape_U_depth[i,n,c], nsamples))
		append!(D_depth_list, fill(shape_data.shape_D_depth[i,n,c], nsamples))
	end
end

# ╔═╡ c25e20fe-d5f2-4cbd-91a6-40d153e0ab51
let fig = Makie.Figure()
	bins = -2:0.1:5
	ax = Makie.Axis(fig[1,1]; width=500, height=400, xlabel="Reactivities sampled from posterior P(r | r_in)", ylabel="Posterior prob.")
	
	lb, ub = 0e4, 0.25e4
	Makie.hist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, color=(:blue, 0.1))
	Makie.stephist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, label="Read depths ∈ [0, 2e4]", color=:blue, linewidth=1)

	lb, ub = 2e4, 4e4
	Makie.hist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, color=(:orange, 0.1))
	Makie.stephist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, label="Read depths ∈ [2e4, 4e4]", color=:orange, linewidth=1)

	lb, ub = 4e4, 6e4
	Makie.hist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, color=(:red, 0.1))
	Makie.stephist!(ax, R_list[(lb .< M_depth_list .< ub) .&& (lb .< U_depth_list .< ub) .&& (lb .< D_depth_list .< ub)]; normalization=:pdf, bins, label="Read depths ∈ [4e4, 6e4]", color=:red, linewidth=1)

	Makie.axislegend(ax; position=:rt)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═468738a8-3a61-11f0-0206-291073c6db4d
# ╠═bdb8966d-ee65-4cc1-8e03-95edc7248551
# ╠═658362d5-eb92-473d-9f2b-ea83429e9363
# ╠═002d9917-53e7-462f-8b81-553e859c9d73
# ╠═5e8b92b5-66f0-4f28-8010-c103de66d1b4
# ╠═b4804a68-60fd-4655-a9a6-b0b82f4f70c4
# ╠═dc281c13-24ae-481a-918d-fb69e3392435
# ╠═a41e0a5b-f029-4b1f-a386-570d5603391b
# ╠═2b5f7847-bec2-4cd1-b95d-d04dc5bec30e
# ╠═3899ac6f-e9ba-426a-a916-17eb40b54847
# ╠═8ac4b744-7c1c-466c-92a7-72bd315d7b83
# ╠═f1bd692c-fc59-4c3c-8f15-021a25201e85
# ╠═f48e3d7a-24fc-4ec0-a031-043b256f5541
# ╠═24957a09-2b7c-46cf-8e8c-dff49186fd4d
# ╠═562a46fb-d29e-43b2-99c3-ca5e0691765b
# ╠═c7a6a735-8c5b-497c-91eb-8016b50fd8d2
# ╠═a97ff5ca-c6b7-4bb5-b5c4-09aafb85a8e3
# ╠═61417d18-0bd1-407f-8e33-fc8d1c3796c7
# ╠═98e0c738-904b-4778-aa0b-c4a294b80d77
# ╠═c25e20fe-d5f2-4cbd-91a6-40d153e0ab51

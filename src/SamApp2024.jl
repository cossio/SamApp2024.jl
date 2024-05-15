module SamApp2024

import CSV
import FASTX
import LazyArtifacts
using Base: Fix2
using BioSequences: LongRNA
using BioSequences: LongSequence
using DataFrames: DataFrame
using LazyArtifacts: @artifact_str
using LinearAlgebra: Diagonal
using RestrictedBoltzmannMachines: batchmean
using RestrictedBoltzmannMachines: flat_w
using RestrictedBoltzmannMachines: Potts
using RestrictedBoltzmannMachines: RBM
using RestrictedBoltzmannMachines: var_h_from_v

include("effective_contacts.jl")
include("shape_500.jl")
include("hamming.jl")
include("rfam.jl")
include("onehot.jl")

end

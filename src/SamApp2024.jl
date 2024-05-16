module SamApp2024

import CSV
import FASTX
import LazyArtifacts
import HDF5
import Infernal
import Rfam
using Base: Fix2
using BioSequences: LongRNA
using BioSequences: LongSequence
using BioSequences: @rna_str
using DataFrames: DataFrame
using LazyArtifacts: @artifact_str
using LinearAlgebra: Diagonal
using RestrictedBoltzmannMachines: batchmean
using RestrictedBoltzmannMachines: flat_w
using RestrictedBoltzmannMachines: Potts
using RestrictedBoltzmannMachines: RBM
using RestrictedBoltzmannMachines: var_h_from_v
using RestrictedBoltzmannMachines: dReLU
using XLSX: readxlsx

include("effective_contacts.jl")
include("shape_500.jl")
include("hamming.jl")
include("rfam.jl")
include("onehot.jl")
include("artifacts.jl")
include("shape.jl")
include("shape_alignment.jl")
include("secondary_structure.jl")

end

module SamApp2024

import Artifacts
import CairoMakie
import CSV
import DCAUtils
import DelimitedFiles
import Distributions
import FASTX
import HDF5
import Infernal
import JLD2
import KernelDensity
import LazyArtifacts
import LoggingExtras
import Logomaker
import Makie
import NaNStatistics
import Rfam
import StatsBase
import StringDistances
import ViennaRNA
import ViennaRNA_jll
using Base: Fix2
using BioSequences: @rna_str
using BioSequences: LongRNA
using BioSequences: LongSequence
using DataFrames: DataFrame
using LazyArtifacts: @artifact_str
using LinearAlgebra: Diagonal
using RestrictedBoltzmannMachines: batchmean
using RestrictedBoltzmannMachines: dReLU
using RestrictedBoltzmannMachines: flat_w
using RestrictedBoltzmannMachines: Potts
using RestrictedBoltzmannMachines: RBM
using RestrictedBoltzmannMachines: var_h_from_v
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

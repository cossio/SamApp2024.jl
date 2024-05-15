module SamApp2024

import CSV

using Base: Fix2
using BioSequences: LongRNA
using LinearAlgebra: Diagonal
using RestrictedBoltzmannMachines: batchmean
using RestrictedBoltzmannMachines: flat_w
using RestrictedBoltzmannMachines: Potts
using RestrictedBoltzmannMachines: RBM
using RestrictedBoltzmannMachines: var_h_from_v

include("effective_contacts.jl")

end

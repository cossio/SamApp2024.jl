"""
    rbm2022()

Loads the RBM used to generate sequences, probed in the experiments in 2022.
"""
function rbm2022()
    rbm = HDF5.h5open(joinpath(artifact"RBM2022", "RBM_nh_100_l1b_0.01.hdf5"), "r") do h5py
        RBM(
            Potts(θ = HDF5.read(h5py, "v_fields")),
            dReLU(;
                θp = -HDF5.read(h5py, "h_theta_plus"), θn = HDF5.read(h5py, "h_theta_minus"),
                γp = HDF5.read(h5py, "h_gamma_plus"), γn = HDF5.read(h5py, "h_gamma_minus")
            ),
            HDF5.read(h5py, "weights")
        )
    end
    @assert size(rbm.visible) == (5, 108)
    @assert size(rbm.hidden) == (100,)
    @assert size(rbm.w) == (5, 108, 100)
    return rbm
end

function probed_aptamers_table_20221027()
    path = joinpath(artifact"probed_aptamers_table_20221027", "2022-10-27-tagged-aptamer-sequences-natural_artifical_pdb.tsv")
    aptamers_df = CSV.read(path, DataFrame)
    return aptamers_df
end

"""
    pierre20221107post(:natural)
    pierre20221107post(:synthetic)

My post-processing of the Excel file 2022-11-07_APSAM-I_data_analysis sent by Pierre.
There's one table for the `:natural` aptamers, and another for the `:synthetic` ones.
"""
function pierre20221107post(name::Symbol)
    if name === :natural
        return joinpath(artifact"pierre20221107post", "natural_aptamers.csv")
    elseif name === :synthetic
        return joinpath(artifact"pierre20221107post", "synthetic_aptamers.csv")
    else
        throw(ArgumentError("Argument should be `:natural` or `:synthetic`"))
    end
end

"""
    probed_artificial_sequences_2022_df()

Loads the Excel file joint_proposed_RF00162_annotated.xlsx as a DataFrame.
"""
function probed_artificial_sequences_2022_df()
    xls = readxlsx(probed_artificial_sequences_2022_excel())["joint_proposed_RF00162 (2)"]
    df = DataFrame([xls[:][1,n] => xls[:][2:end,n] for n in 1:size(xls[:], 2)]);
    df.RBM_energy = parse.(Float64, df.RBM_energy)
    df.INF_score = parse.(Float64, df.INF_score)
    return df
end

"""
    probed_artificial_sequences_2022_excel()

Excel table joint_proposed_RF00162_annotated.xlsx of the sequences we designed for probing in 2022.
RBM sequences were generated with `rbm2022()`. Infernal sequences with the model downloaded
from Rfam (which has the entropic noise added).
"""
function probed_artificial_sequences_2022_excel()
    return joinpath(artifact"Pierre20221107ShapeAnnotations", "joint_proposed_RF00162_annotated.xlsx")
end

#= Samples generated from the rbm2022 model (equilibrated). =#
function rbm2022samples()
    sampled_v = HDF5.h5open(joinpath(artifact"rbm2022-samples", "rbm_samples.hdf5"), "r") do hdf5
        HDF5.read_dataset(hdf5, "samples")
    end
    return BitArray(sampled_v)
end
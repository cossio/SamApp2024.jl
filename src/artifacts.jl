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

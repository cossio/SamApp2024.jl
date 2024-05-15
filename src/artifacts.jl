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

function load_dms_data_20250609_pdb_4kqy()
    dms_dir = artifact_dir_DMS_PDB_4KQY_20250607()
    conditions = filter(startswith("SAMAP-P1_DMS"), readdir(dms_dir))
    @assert issorted(conditions)

    positions_mapping_2022 = shape_positions_alignment_2022()
    sequence_length = 108
    number_of_sequences = 1

    RF00162_hits_afa_1410 = Infernal.cmalign(
        Infernal.cmfetch(Rfam.cm(; rfam_version="14.10"), "RF00162").out,
        Rfam.fasta_file("RF00162"; rfam_version="14.10");
        matchonly=true, outformat="afa"
     )
     RF00162_hits_dsc_1410 = FASTX.description.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out)))
     RF00162_hits_seqs_1410 = LongRNA{4}.(FASTX.sequence.(FASTX.FASTA.Reader(open(RF00162_hits_afa_1410.out))))

    aligned_sequence = [string(RF00162_hits_seqs_1410[only(i for (i, dsc) = enumerate(RF00162_hits_dsc_1410) if occursin("4KQY", dsc))])]

    shape_M = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_U = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_D = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_M_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_U_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_D_depth = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_reactivities = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_reactivities_err = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    shape_raw_reactivities = fill(NaN, sequence_length, number_of_sequences, length(conditions))
    shape_raw_reactivities_err = fill(NaN, sequence_length, number_of_sequences, length(conditions))

    for (c, cond) = enumerate(conditions), n = 1:number_of_sequences
        profile_file = joinpath(dms_dir, cond, "$(cond)_APSAM-PDB10-P1_profile.txt")
        profile_df = CSV.read(profile_file, DataFrame)
        profile_df_sequence = join(profile_df.Sequence)
        for i = 1:108
            mapped_position = positions_mapping_2022.pdb[2,i]
            if ismissing(mapped_position)
                continue
            else
                @assert only(aligned_sequence)[i] == profile_df_sequence[mapped_position]

                shape_M[i, n, c] = profile_df.Modified_rate[mapped_position]
                shape_U[i, n, c] = profile_df.Untreated_rate[mapped_position]
                shape_D[i, n, c] = profile_df.Denatured_rate[mapped_position]

                shape_M_depth[i, n, c] = profile_df.Modified_effective_depth[mapped_position]
                shape_U_depth[i, n, c] = profile_df.Untreated_effective_depth[mapped_position]
                shape_D_depth[i, n, c] = profile_df.Denatured_effective_depth[mapped_position]

                shape_reactivities[i, n, c] = profile_df.HQ_profile[mapped_position]
                shape_reactivities_err[i, n, c] = profile_df.HQ_stderr[mapped_position]

                shape_raw_reactivities[i, n, c] = profile_df.Reactivity_profile[mapped_position]
                shape_raw_reactivities_err[i, n, c] = profile_df.Std_err[mapped_position]
            end
        end
    end

    shape_M_stderr = sqrt.(shape_M ./ shape_M_depth)
    shape_U_stderr = sqrt.(shape_U ./ shape_U_depth)
    shape_D_stderr = sqrt.(shape_D ./ shape_D_depth)

    return (;
        shape_M, shape_M_depth, shape_M_stderr,
        shape_U, shape_U_depth, shape_U_stderr,
        shape_D, shape_D_depth, shape_D_stderr,
        shape_reactivities, shape_reactivities_err,
        shape_raw_reactivities, shape_raw_reactivities_err,
        conditions, aligned_sequence
    )
end

function infernal_align_fasta_to_cm(fasta_path::AbstractString, cm_path::AbstractString)
    # trimmed (no inserts) aligned fasta
    aln = Infernal.cmalign(cm_path, fasta_path; matchonly=true, outformat="AFA")

    # these are already aligned and without inserts
    sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(aln.out)))
    @assert allequal(map(length, sequences))

    return LongRNA{4}.(sequences)
end

function infernal_score_sequences(cm_path::AbstractString, sequences::AbstractVector; informat=nothing, notrunc=true, glob=true)
    mktemp() do fasta_path, io
        FASTX.FASTA.Writer(io) do writer
            for (i, sequence) = enumerate(sequences)
                record = FASTX.FASTA.Record(string(i), sequence)
                write(writer, record)
            end
        end
        result = Infernal.cmalign(cm_path, fasta_path; glob, informat, notrunc)
        return Infernal.cmalign_parse_sfile(result.sfile)
    end
end

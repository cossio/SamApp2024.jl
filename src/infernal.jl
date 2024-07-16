function align_fasta_to_cm(fasta_path::AbstractString, cm_path::AbstractString)
    # trimmed (no inserts) aligned fasta
    aln = Infernal.cmalign(cm_path, fasta_path; matchonly=true, outformat="AFA")

    # these are already aligned and without inserts
    sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(aln.out)))
    @assert allequal(map(length, sequences))

    return LongRNA{4}.(sequences)
end

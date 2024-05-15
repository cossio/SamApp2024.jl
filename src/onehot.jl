#=
Define the order in which amino-acids and nucleotides are encoded as integers.
=#

function alphabet(::Type{<:LongAA}; gap::Bool = true)
    AAs = aa"ACDEFGHIKLMNPQRSTVWY"
    if gap
        return AAs * aa"-"
    else
        return AAs
    end
end

function alphabet(::Type{<:LongRNA}; gap::Bool = true)
    NTs = rna"ACGU"
    if gap
        return NTs * rna"-"
    else
        return NTs
    end
end

function alphabet(::Type{<:LongDNA}; gap::Bool = true)
    NTs = dna"ACGT"
    if gap
        return NTs * rna"-"
    else
        return NTs
    end
end

function alphabet(s::Union{<:LongAA, <:LongRNA{4}, <:LongDNA{4}}; gap::Bool = true)
    return alphabet(typeof(s); gap)
end

#=
Encodes sequences as Potts arrays, where the integer entry indicates the sequence letter.
=#

function potts(X::BitMatrix)
    return vec(Int8.(first.(Tuple.(argmax(X; dims=1)))))
end

function potts(X::BitArray{3})
    return reshape(Int8.(first.(Tuple.(argmax(X; dims=1)))), size(X,2), size(X,3))
end

function potts(s::Union{LongAA, LongAA, LongRNA{4}, LongDNA{4}, AbstractVector{<:LongAA}, AbstractVector{<:LongRNA}, AbstractVector{<:LongDNA}})
    potts(onehot(s))
end

function aaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongAA([alphabet(LongAA; gap)[i] for i in P])
end

function rnaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongRNA{4}([alphabet(LongRNA; gap)[i] for i in P])
end

function dnaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongDNA{4}([alphabet(LongDNA; gap)[i] for i in P])
end

function aaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [aaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function rnaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [rnaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function dnaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [dnaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function onehot(seq::Union{LongAA, LongRNA{4}, LongDNA{4}}; gap::Bool=true)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(alphabet(seq; gap))
end

function onehot(seqs::Union{AbstractVector{<:LongAA}, AbstractVector{<:LongRNA}, AbstractVector{<:LongDNA}}; gap::Bool=true)
    L = only(unique(length.(seqs))) # all sequences must have same length
    return reshape(mapreduce(s -> onehot(s; gap), hcat, seqs), :, L, length(seqs))
end

function aaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 21
        gap = true
    elseif size(X, 1) == 20
        gap = false
    else
        error("Expected 20 or 21 rows; got $(size(X, 1))")
    end

    return aaseq(potts(X); gap)
end

function rnaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 5
        gap = true
    elseif size(X, 1) == 4
        gap = false
    else
        error("Expected 4 or 5 rows; got $(size(X, 1))")
    end

    return rnaseq(potts(X); gap)
end

function dnaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 5
        gap = true
    elseif size(X, 1) == 4
        gap = false
    else
        error("Expected 4 or 5 rows; got $(size(X, 1))")
    end

    return dnaseq(potts(X); gap)
end

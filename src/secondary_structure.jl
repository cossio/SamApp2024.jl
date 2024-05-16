"""
    RF00162_sites_paired()

Returns the lists of paired, unpaired, and pseudoknot sites in the consensus secondary
structure of RF00162.
"""
function RF00162_sites_paired()
    wuss = RF00162_wuss(; insertions=false)
    bps = findall(∈("()<>"), wuss)
    nps = findall(∉("()<>Aa"), wuss)
    pks = findall(∈("Aa"), wuss)
    return (; bps, nps, pks)
end

"""
    RF00162_wuss(; insertions=false)

Fetch the WUSS secondary structure of RF00162 from Rfam.
"""
function RF00162_wuss(; insertions=false)
    wuss_full = stockholm_ss(Infernal.esl_afetch(Rfam.seed(), "RF00162").out)
    if insertions
        return wuss_full
    else
        # remove insertions
        return filter(!=('.'), wuss_full)
    end
end

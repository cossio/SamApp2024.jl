import FASTX
import Infernal
import Rfam
import SamApp2024
using BioSequences: LongRNA

_tmpfasta = tempname()

FASTX.FASTA.Writer(open(_tmpfasta, "w")) do writer
    for (n, seq) in enumerate(SamApp2024.rnaseq(SamApp2024.rbm2022samples()))
        ismissing(seq) && continue
        write(writer, FASTX.FASTA.Record("RBM-$n", filter(!=('-'), string(seq))))
    end
end

Rfam_cm = Infernal.cmfetch(Rfam.cm(), "RF00162")
Rfam_cm_emitted_sequences_afa = Infernal.cmemit(Rfam_cm.out; N=5000, aligned=true, outformat="AFA")
Rfam_cm_emitted_sequences = FASTX.sequence.(FASTX.FASTA.Reader(open(Rfam_cm_emitted_sequences_afa.out)))
Rfam_cm_emitted_sequences = [filter(!=('.'), filter(!islowercase, seq)) for seq = Rfam_cm_emitted_sequences]
Rfam_cm_emitted_sequences = LongRNA{4}.(Rfam_cm_emitted_sequences)

_tmpfasta_cm = tempname()

FASTX.FASTA.Writer(open(_tmpfasta_cm, "w")) do writer
    for (n, seq) in enumerate(Rfam_cm_emitted_sequences)
        write(writer, FASTX.FASTA.Record("rCM-$n", filter(!=('-'), string(seq))))
    end
end

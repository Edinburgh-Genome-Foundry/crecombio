import crecombio

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Translocation
seq0 = SeqRecord(
        Seq(
            "AAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAA"
        )
    )
seq1 = SeqRecord(
        Seq(
            "TTTTTTTTTTTTTTTTTTTTGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTT"
        )
    )

recombined_seqs = crecombio.recombine_two_sequences([seq0, seq1])

print(recombined_seqs[0].seq)
# AAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTT
print(recombined_seqs[1].seq)
# TTTTTTTTTTTTTTTTTTTTGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAA

# Insertion
excision_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "AAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
target_seq = SeqRecord(
        Seq(
            "GGGGGGGGGGGGGGGGGGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTTTTTTT"
        )
    )
recombined_seqs = crecombio.recombine_two_sequences([excision_seq, target_seq])
print(recombined_seqs[0].seq)
# GGGGGGGGGGGGGGGGGGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTTTTTTT
print(recombined_seqs[1].seq)
# GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC

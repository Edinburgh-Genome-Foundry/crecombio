from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import crecombio


def test_recombine():
    excision_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "AAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    results = crecombio.recombine_one_sequence([excision_seq])

    assert (
        str(results[0][0].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCC"
        "CCCCCCCCCCCCC"
    )

    inversion_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "GGGGGGGGGGGGGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    results = crecombio.recombine_one_sequence([inversion_seq])
    assert (
        str(results[0][0].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCT"
        "TTTTTTTTTTTTGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )

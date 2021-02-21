import pytest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import crecombio


def test_recombine_one_sequence():
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


def test_recombine():
    excision_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "GGGGGGGGGGGGGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    crecombio.recombine([excision_seq])  # 1 seq
    crecombio.recombine([excision_seq, excision_seq])  # 2 seq
    with pytest.raises(Exception):
        crecombio.recombine(
            [excision_seq, excision_seq, excision_seq]
        )  # not 1 or 2 seq

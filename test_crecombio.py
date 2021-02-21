import pytest

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import crecombio


def test_count_number_of_sites():
    seq_with_three_sites = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "AAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            "GGGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCC"
        )
    )
    site = crecombio.SITES["Flp"]["seq"]
    matches, rc_matches = crecombio.count_number_of_sites(seq_with_three_sites, site)
    assert matches == [26, 82]
    assert rc_matches == [145]


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


def test_recombine_two_sequences():
    seq0 = SeqRecord(
        Seq("AAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAA")
    )
    seq1 = SeqRecord(
        Seq("TTTTTTTTTTTTTTTTTTTTGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTT")
    )
    recombined_seqs = crecombio.recombine_two_sequences([seq0, seq1])
    assert (
        str(recombined_seqs[0].seq)
        == "AAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTT"
    )
    assert (
        str(recombined_seqs[1].seq)
        == "TTTTTTTTTTTTTTTTTTTTGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAA"
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

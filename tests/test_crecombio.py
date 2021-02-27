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
    # Not 2 FRT sites:
    no_frt_seq = SeqRecord(Seq("ATCG"))
    with pytest.raises(Exception):
        crecombio.recombine_one_sequence([no_frt_seq])

    # FRTs in same direction:
    excision_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "AAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    results = crecombio.recombine_one_sequence([excision_seq])

    assert (
        str(results[0].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCC"
        "CCCCCCCCCCCCC"
    )

    # FRTs in opposite direction:
    inversion_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "GGGGGGGGGGGGGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    results = crecombio.recombine_one_sequence([inversion_seq])
    assert (
        str(results[0].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCT"
        "TTTTTTTTTTTTGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
    )


def test_recombine_two_sequences():
    # Test translocation:
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

    # Test insertion:
    excision_seq = SeqRecord(
        Seq(
            "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAA"
            "AAAAAAAAAGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        )
    )
    target_seq = SeqRecord(
        Seq("GGGGGGGGGGGGGGGGGGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTTTTTTT")
    )
    # Excision seq first, target seq second:
    recombined_seqs = crecombio.recombine_two_sequences([excision_seq, target_seq])
    assert (
        str(recombined_seqs[0].seq)
        == "GGGGGGGGGGGGGGGGGGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAAAAAAAAAAAGA"
        "AGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTTTTTTT"
    )
    assert (
        str(recombined_seqs[1].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCC"
        "CCCCCCCCCCCCC"
    )

    # Target seq first, excision seq second:
    recombined_seqs = crecombio.recombine_two_sequences([target_seq, excision_seq])
    assert (
        str(recombined_seqs[0].seq)
        == "GGGGGGGGGGGGGGGGGGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCAAAAAAAAAAAAAAAAAAAAAAGA"
        "AGTTCCTATTCTCTAGAAAGTATAGGAACTTCTTTTTTTTTTTTTTTTT"
    )
    assert (
        str(recombined_seqs[1].seq)
        == "GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCC"
        "CCCCCCCCCCCCC"
    )

    # 2x2 FRT sites:
    with pytest.raises(Exception):
        crecombio.recombine_two_sequences([excision_seq, excision_seq])

    # Incorrect number of FRT sites:
    with pytest.raises(Exception):
        crecombio.recombine_two_sequences(
            [SeqRecord(Seq("ATCG")), SeqRecord(Seq("ATCG"))]
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

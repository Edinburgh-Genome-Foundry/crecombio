import re

from Bio.Seq import Seq

SITES = {
    "Flp": {
        "seq": Seq("GAAGTTCCTATTCtctagaaaGtATAGGAACTTC".upper()),
        "5p": Seq("GAAGTTCCTATTC"),
        "3p": Seq("tctagaaaGtATAGGAACTTC".upper()),
    },
    "Cre": {"seq": Seq("ATAACTTCGTATAATGTATGCTATACGAAGTTAT")},
    "Cre_lox_511": {"seq": Seq("ATAACTTCGTATAATGTATACTATACGAAGTTAT")},
    "Cre_lox_5171": {"seq": Seq("ATAACTTCGTATAATGTGTACTATACGAAGTTAT")},
    "Cre_lox_2272": {"seq": Seq("ATAACTTCGTATAAAGTATCCTATACGAAGTTAT")},
    "Cre_M2": {"seq": Seq("ATAACTTCGTATAAGAAACCATATACGAAGTTAT")},
    "Cre_M3": {"seq": Seq("ATAACTTCGTATATAATACCATATACGAAGTTAT")},
    "Cre_M7": {"seq": Seq("ATAACTTCGTATAAGATAGAATATACGAAGTTAT")},
    "Cre_M11": {"seq": Seq("ATAACTTCGTATACGATACCATATACGAAGTTAT")},
}


def recombine(sequences, recombinase="Flp"):
    """Simulate recombination and return sequences.


    **Parameters**

    **sequences**
    > List of one or two Biopython SeqRecord objects (`list`).

    **recombinases**
    > Recombinase name (`str`). Default Flp.
    """
    if len(sequences) == 1:
        return recombine_one_sequence(sequences, recombinase)
    if len(sequences) == 2:
        return recombine_two_sequences(sequences, recombinase)
    else:
        raise Exception("Only 1- or 2-sequence recombination is supported")


def recombine_one_sequence(sequences, recombinase="Flp"):
    """Simulate recombination of one sequence."""
    recombined_sequences = []
    for sequence in sequences:
        site = SITES[recombinase]["seq"]
        matches, rc_matches = count_number_of_sites(sequence, site)
        if (len(matches) + len(rc_matches)) != 2:
            raise Exception("The sequence must contain 2 recombinase sites")

        # Excision
        if len(matches) == 2:
            new_sequence = sequence[: matches[0]] + sequence[matches[1] :]
            recombined_sequences += [new_sequence]
        # Excision on reverse strand
        if len(rc_matches) == 2:
            new_sequence = sequence[: matches[0]] + sequence[matches[1] :]
            recombined_sequences += [new_sequence]
        # Inversion
        if len(matches) == 1 and len(rc_matches) == 1:
            if matches[0] > rc_matches[0]:
                raise Exception(
                    "Orientation incorrect: reverse complement of recombination site is"
                    " before the site."
                )
                # extend to circular sequences later
            inside_sequence = sequence[
                (matches[0] + len(site)) : rc_matches[0]
            ]  # added length to not include the site
            inverted_sequence = inside_sequence.reverse_complement()
            left_part = sequence[: matches[0] + len(site)]
            right_part = sequence[rc_matches[0] :]
            new_sequence = left_part + inverted_sequence + right_part
            recombined_sequences += [new_sequence]

    return recombined_sequences


def recombine_two_sequences(sequences, recombinase="Flp"):
    """Simulate recombination of two sequences."""
    seq0 = sequences[0]
    seq1 = sequences[1]
    site = SITES[recombinase]["seq"]
    seq0_matches, seq0_rc_matches = count_number_of_sites(seq0, site)
    seq1_matches, seq1_rc_matches = count_number_of_sites(seq1, site)
    if len(seq0_rc_matches) != 0 or len(seq1_rc_matches) != 0:
        print(
            "Simulating recombination of reverse complement sites are not implemented "
            "yet and these sites are ignored"
        )

    recombined_sequences = []
    if len(seq0_matches) == 1 and len(seq1_matches) == 1:
        # Translocate
        seq0_left_part = seq0[: seq0_matches[0]]
        seq0_right_part = seq0[seq0_matches[0] :]  # includes the site
        seq1_left_part = seq1[: seq1_matches[0]]
        seq1_right_part = seq1[seq1_matches[0] :]  # includes the site
        seq2 = seq0_left_part + seq1_right_part
        seq3 = seq1_left_part + seq0_right_part
        recombined_sequences += [seq2, seq3]

    elif len(seq0_matches) == 1 and len(seq1_matches) == 2:
        # Insert from seq1 to seq0
        donor = seq1
        acceptor = seq0
        donor_matches = seq1_matches
        acceptor_matches = seq0_matches
        new_donor = donor[: donor_matches[0]] + donor[donor_matches[1] :]
        excised_seq = donor[donor_matches[0] : donor_matches[1]]
        new_acceptor = (
            acceptor[: acceptor_matches[0]]
            + excised_seq
            + acceptor[acceptor_matches[0] :]
        )
        recombined_sequences += [new_acceptor, new_donor]

    elif len(seq0_matches) == 2 and len(seq1_matches) == 1:
        # Insert from seq0 to seq1
        donor = seq0
        acceptor = seq1
        donor_matches = seq0_matches
        acceptor_matches = seq1_matches
        # Same code as above:
        new_donor = donor[: donor_matches[0]] + donor[donor_matches[1] :]
        excised_seq = donor[donor_matches[0] : donor_matches[1]]
        new_acceptor = (
            acceptor[: acceptor_matches[0]]
            + excised_seq
            + acceptor[acceptor_matches[0] :]
        )
        recombined_sequences += [new_acceptor, new_donor]

    elif len(seq0_matches) == 2 and len(seq1_matches) == 2:
        raise Exception("Both sequences have 2 recombination sites")

    else:
        raise Exception(
            "The number of recombination sites must be 1 or 2 in each sequence"
        )

    return recombined_sequences


def count_number_of_sites(sequence, site):
    """Count the number of recombination sites in a sequence.


    **Parameters**

    **sequence**
    > A Biopython SeqRecord object.

    **site**
    > The site sequence (`str`).
    """
    matches = [
        m.start() for m in re.finditer(re.escape(str(site)), str(sequence.seq).upper())
    ]

    rc_matches = [
        m.start()
        for m in re.finditer(
            re.escape(str(site.reverse_complement())), str(sequence.seq).upper()
        )
    ]

    return matches, rc_matches

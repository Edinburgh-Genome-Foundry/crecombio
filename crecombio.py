import re

from Bio.Seq import Seq

SITES = {
    "Flp": {
        "seq": Seq("GAAGTTCCTATTCtctagaaaGtATAGGAACTTC".upper()),
        "5p": Seq("GAAGTTCCTATTC"),
        "3p": Seq("tctagaaaGtATAGGAACTTC".upper()),
    }
}


def recombine(sequences, recombinase="Flp"):
    """Simulate recombination and return sequences.


    **Parameters**

    **sequences**
    > List of Biopython SeqRecord objects (`list`).

    **recombinases**
    > Recombinase (`str`). Default Flp.
    """
    if len(sequences) != 1:
        raise Exception("Multisequence recombination not supported yet")

    recombined_sequences = []
    for sequence in sequences:
        site = SITES[recombinase]["seq"]
        matches = [
            m.start()
            for m in re.finditer(re.escape(str(site)), str(sequence.seq).upper())
        ]
        print(matches)

        rc_matches = [
            m.start()
            for m in re.finditer(
                re.escape(str(site.reverse_complement())), str(sequence.seq).upper()
            )
        ]

        if (len(matches) + len(rc_matches)) != 2:
            raise Exception("The sequence must contain 2 recombinase sites")

        # Excision
        if len(matches) == 2:
            new_sequence = sequence[: matches[0]] + sequence[matches[1] :]
            recombined_sequences += [[new_sequence]]
        # Excision on reverse strand
        if len(rc_matches) == 2:
            new_sequence = sequence[: matches[0]] + sequence[matches[1] :]
            recombined_sequences += [[new_sequence]]
        # Inversion
        if len(matches) == 1 and len(rc_matches) == 1:
            if matches[0] > rc_matches[0]:
                raise Exception(
                    "Orientation incorrect: reverse complement of site is before the site."
                )
                # extend to circular sequences later
            inside_sequence = sequence[
                (matches[0] + len(site)) : rc_matches[0]
            ]  # added length to not include the site
            inverted_sequence = inside_sequence.reverse_complement()
            left_part = sequence[: matches[0] + len(site)]
            right_part = sequence[rc_matches[0] :]
            new_sequence = left_part + inverted_sequence + right_part
            recombined_sequences += [[new_sequence]]

    return recombined_sequences

<p align="center">
<img alt="Crecombio logo" title="GeneAlloy" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Crecombio/main/images/crecombio.png" width="150">
</p>

# Crecombio

[![Build Status](https://github.com/Edinburgh-Genome-Foundry/crecombio/actions/workflows/build.yml/badge.svg)](https://github.com/Edinburgh-Genome-Foundry/crecombio/actions/workflows/build.yml)
[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/crecombio/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/crecombio?branch=main)

A simple **Cre**, Flp and other site-specific re**combi**nation simulat**o**r.

## Background

Site-specific *recombinase proteins* rearrange one or two DNA molecules, through recombination
at certain recognized sequences. In *Flp-FRT recombination*, a tyrosine family member recombinase protein, *flippase*, recognizes *FRT* (flippase recognition target) sites and performs inversion, excision, insertion or translocation of DNA in a reversible manner, depending on the configuration of the FRT sites. Another, directional, system utilizes *serine integrases* to catalyze rearrangements at DNA sequences called attachment *(att)* sites.

Site-specific recombinations can be categorised in many ways:

- directionality: reversible or not
- recombination enzyme used
- recombination sequence used
- number of DNA molecules present in the reaction and number of recombination sites in each molecule
- purpose: inversion, excision, insertion and translocation

For simulating homologous recombination and other assemblies, use [DNA Cauldron](https://github.com/Edinburgh-Genome-Foundry/dnacauldron) and EGF CUBA [Simulate Golden Gate Assemblies](https://cuba.genomefoundry.org/simulate_gg_assemblies) / [Simulate multi-method assemblies](https://cuba.genomefoundry.org/simulate_multi_method_assemblies).

Crecombio is currently intended for simple 1- or 2-molecule recombination simulations. For more complicated procedures, such as the *serine integrase recombinational assembly* (SIRA), described in *Merrick et al.* (Serine Integrases: Advancing Synthetic Biology. [ACS Synth. Biol. 2018, 7, 299−310](https://pubs.acs.org/doi/10.1021/acssynbio.7b00308)), use it in an iterative manner. Alternatively DNA Cauldron's classes can be used with custom-defined enzymes to simulate the cleavage and recombination.

**Work in progress:**

- Simulating *att* recombinations
- Searching nonspecific (ambiguous) sites
- Handling circular sequences
- Summary of simulation in a PDF report

<p align="left">
<img alt="Flp/FRT recombination" title="Flp/FRT recombination" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Crecombio/main/images/crecombio_plot.jpg" width="1000">
</p>

## Install

```
pip install crecombio
```

## Usage

```python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import crecombio
excision_seq = SeqRecord(Seq("GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCtctagaaaGtATAGGAACTTCAAAAAAAAAAAAAAAAAAAAAAGAAGTTCCTATTCtctagaaaGtATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC".upper()))
results = crecombio.recombine([excision_seq])
print(results[0].seq)
# GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC

inversion_seq = SeqRecord(Seq("GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCtctagaaaGtATAGGAACTTCAAAAAAAAAAAAAGGGGGGGGGGGGGAAGTTCCTATaCtttctagaGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC".upper()))
results = crecombio.recombine([inversion_seq])
print(results[0].seq)
# GACTGATGTGACGTGTGACAGCTGACGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCCCCCCCCCCCCTTTTTTTTTTTTTGAAGTTCCTATACTTTCTAGAGAATAGGAACTTCCCCCCCCCCCCCCCCCCCCCCCCCCCC
```

See the ``examples`` folder for two-sequence recombination examples.

## License = MIT

Crecombio is [free software](https://www.gnu.org/philosophy/free-sw.en.html), which means the users have the freedom to run, copy, distribute, study, change and improve the software.

Crecombio was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/) by [Peter Vegh](https://github.com/veghp) and is released under the MIT license.

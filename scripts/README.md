# Supplemental scripts

## `draw_biogeochemical_cycles.R`

R scritps edited from [R script](https://github.com/AnantharamanLab/METABOLIC/blob/master/draw_biogeochemical_cycles.R) to create biogeochemical diagram cycles. I added Acetogenesis and Nitric Oxide Dismutase and several colors to identify major functions.

To run it, just use the command `Rscript draw_biogeochemical_cycles.R bigecyhmm_output_folder/diagram_input_folder/ diagram_output TRUE` on the diagram folder (`diagram_input_folder`) of the results folder of bigecyhmm.

## `create_hmm.py`

Little python script to create HMM file from a protein fasta files. The input protein sequences are associated with a specific genes.

It requires the following package:

- [biopython](https://github.com/biopython/biopython)
- [pyhmmer](https://github.com/althonos/pyhmmer)
- [pyfamsa](https://github.com/althonos/pyfamsa)
- [pytrimal](https://github.com/althonos/pytrimal)

A multiple alignments of the protein is performed with `pyfamsa`. Then this alignement is trimmed with `pytrimal`. Finally a HMM file is created with `pyhmmer`. `biopython` is used to read the input file.

To use it ru nthe following command: `python3 create_hmm.py -i proteins_sequences.faa`.

This will create a `proteins_sequences.hmm` file.
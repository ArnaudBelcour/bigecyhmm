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

Added HMMs are:

- `aspA` (for Adenosine-5′-Phosphosulfate Reductase), sequences extracted from the article of [Friedrich 2002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC134748/).
- `nod` (for nitric oxide dismutase), sequences extracted from the article of [Zhu et al. 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6636425/) and from [Schmitz et al. 2023](https://journals.asm.org/doi/10.1128/msphere.00571-22).
- `fthfs` (for Formate-tetrahydrofolate ligase), sequences extracted from UniPort on the 12 March 2024 using the request `(protein_name%3AFTHFS)&facets=reviewed%3Atrue`.
- `mcrG` with KO `K00402`.
- `mtaB` with KO `K04480`.
- `mtbA` with KO `K14082`.
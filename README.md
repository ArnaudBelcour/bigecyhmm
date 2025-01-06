# bigecyhmm: Biogeochemical cycle HMMs search

This is a package to search for genes associated with biogeochemical cycles in protein sequence fasta files. The HMMs come from METABOLIC article, KEGG, PFAM, TIGR.

## Dependencies

bigecyhmm is developed to be as minimalist as possible. It requires:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [Pillow](https://github.com/python-pillow/Pillow): to create biogeochemical cycle diagrams.

The HMMs used are stored inside the package as a zip file ([hmm_files.zip](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/bigecyhmm/hmm_databases)). It makes this python package a little heavy (around 15 Mb) but in this way, you do not have to download other files and can directly use it.

## Installation

It can be installed with pip by cloning the repository:

```sh
git clone https://github.com/ArnaudBelcour/bigecyhmm.git

cd bigecyhmm

pip install -e .

```

## Run bigecyhmm

You can used the tools with two calls:

- by giving as input a protein fasta file:

```sh
bigecyhmm -i protein_sequence.faa -o output_dir
```

- by giving as input a folder containing multiple fasta files:

```sh
bigecyhmm -i protein_sequences_folder -o output_dir
```

There is one option:

* `-c` to indicate the number of core used. It is only useful if you have multiple protein fasta files as the added cores will be used to run another HMM search on a different protein fasta files. 

## Output

It gives as output:

- a folder `hmm_results`: one tsv files showing the hits for each protein fasta file.
- `function_presence.tsv` a tsv file showing the presence/absence of generic functions associated with the HMMs that matched.
- a folder `diagram_input`, the necessary input to create Carbon, Nitrogen, Sulfur and other cycles with the [R script](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/scripts/draw_biogeochemical_cycles.R) modified from the [METABOLIC repository](https://github.com/AnantharamanLab/METABOLIC) using the following command: `Rscript draw_biogeochemical_cycles.R bigecyhmm_output_folder/diagram_input_folder/ diagram_output TRUE`. This script requires the diagram package that could be installed in R with `install.packages('diagram')`.
- a folder `diagram_figures` contains biogeochemical diagram figures drawn from template situated in `bigecyhmm/templates`.


## bigecyhmm_visualisation

There is a second command associated with bigecyhmm (`bigecyhmm_visualisation`), to create visualisation of the results.

To create the associated figures, there are other dependencies:

- seaborn
- pandas
- plotly
- kaleido

Four inputs are expected:

- `--esmecata`: esmecata output folder associated with the run (as the visualisation works on esmecata results).
- `--bigecyhmm`: bigecyhmm output folder associated with the run.
- `--abundance-file`: abundance file indicating the abundance for each organisms selected by EsMeCaTa.
- `-o`: an output folder.


## Citation

If you have used bigecyhmm in an article, please cite:

- this github repository for bigecyhmm.

- PyHMMER for the search on the HMMs:

Martin Larralde and Georg Zeller. PyHMMER: a python library binding to HMMER for efficient sequence analysis. Bioinformatics, 39(5):btad214, May 2023.  https://doi.org/10.1093/bioinformatics/btad214

- HMMer website for the search on the HMMs:

HMMER. http://hmmer.org. Accessed: 2022-10-19.

- the following articles for the creation of the custom HMMs:

Zhou, Z., Tran, P.Q., Breister, A.M. et al. METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks. Microbiome 10, 33 (2022). https://doi.org/10.1186/s40168-021-01213-8

Anantharaman, K., Brown, C., Hug, L. et al. Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system. Nat Commun 7, 13219 (2016). https://doi.org/10.1038/ncomms13219

- the following article for KOfam HMMs:

Takuya Aramaki, Romain Blanc-Mathieu, Hisashi Endo, Koichi Ohkubo, Minoru Kanehisa, Susumu Goto, Hiroyuki Ogata, KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold, Bioinformatics, Volume 36, Issue 7, April 2020, Pages 2251–2252, https://doi.org/10.1093/bioinformatics/btz859

- the following article for TIGRfam HMMs:

Jeremy D. Selengut, Daniel H. Haft, Tanja Davidsen, Anurhada Ganapathy, Michelle Gwinn-Giglio, William C. Nelson, Alexander R. Richter, Owen White, TIGRFAMs and Genome Properties: tools for the assignment of molecular function and biological process in prokaryotic genomes, Nucleic Acids Research, Volume 35, Issue suppl_1, 1 January 2007, Pages D260–D264, https://doi.org/10.1093/nar/gkl1043

- the following article for Pfam HMMs:

Robert D. Finn, Alex Bateman, Jody Clements, Penelope Coggill, Ruth Y. Eberhardt, Sean R. Eddy, Andreas Heger, Kirstie Hetherington, Liisa Holm, Jaina Mistry, Erik L. L. Sonnhammer, John Tate, Marco Punta, Pfam: the protein families database, Nucleic Acids Research, Volume 42, Issue D1, 1 January 2014, Pages D222–D230, https://doi.org/10.1093/nar/gkt1223
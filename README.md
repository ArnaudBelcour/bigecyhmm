# bigecyhmm: Biogeochemical cycle HMMs search

This is a package to search for genes associated with biogeaochemical cycles in protein sequence fasta files. The HMMs come from METABOLIC article, KEGG, PFAM, TIGR.

## Depedencies

I develop bigecyhmm to be as minimalist as possible so it requires only:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform the HMM search.

The searched HMMs are contained inside the package as a zip file ([]()).

## Installation

It can be isntalled by cloning the repository:

```sh
git clone https://github.com/ArnaudBelcour/bigecyhmm.git

cd bigecyhmm
```

And then use a pip install command:

```sh
pip install -e . --config-settings editable_mode=compat

```

## Run bigecyhmm

You can used the tools with two calls:

- by giving as input a protein fasta file:

```sh
bigecyhmm -i protein_sequence.faa -o output_dir -c 10
```

- by giving as input a folder containing multiple fasta files:

```sh
bigecyhmm -i protein_sequences_folder -o output_dir -c 10
```

## Output

It gives as output:

- one tsv files showing the hits for each protein fasta file.
- a tsv file showing the presence/absence of generic functions associated with the HMMs that matched.


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
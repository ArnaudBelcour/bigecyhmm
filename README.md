[![PyPI version](https://img.shields.io/pypi/v/bigecyhmm.svg)](https://pypi.org/project/bigecyhmm/)

# bigecyhmm: Biogeochemical cycle HMMs search

This is a package to search for genes associated with biogeochemical cycles in protein sequence fasta files. The HMMs come from METABOLIC article, KEGG, PFAM, TIGR.

## Table of contents
- [bigecyhmm: Biogeochemical cycle HMMs search](#bigecyhmm-biogeochemical-cycle-hmms-search)
  - [Table of contents](#table-of-contents)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [Run bigecyhmm](#run-bigecyhmm)
  - [Output](#output)
  - [bigecyhmm\_visualisation](#bigecyhmm_visualisation)
    - [Function occurrence and abundance](#function-occurrence-and-abundance)
    - [Output](#output-1)
  - [Citation](#citation)

## Dependencies

bigecyhmm is developed to be as minimalist as possible. It requires:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [Pillow](https://github.com/python-pillow/Pillow): to create biogeochemical cycle diagrams.

The HMMs used are stored inside the package as a zip file ([hmm_files.zip](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/bigecyhmm/hmm_databases)). It makes this python package a little heavy (around 15 Mb) but in this way, you do not have to download other files and can directly use it.

## Installation

It can be installed from PyPI:

`pip install bigecyhmm`

Or it can be installed with pip by cloning the repository:

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
- `bigecyhmm.log`: log file.
- `bigecyhmm_metadata.json`: bigecyhmm metadata (Python version used, package version used).
- `function_presence.tsv`: occurrence of the functions in the different input protein files.
- `pathway_presence.tsv`: occurrence of the major metabolic pathways in the different inputs files.
- `pathway_presence_hmms.tsv`: HMMs with matches for the major metabolic pathways in the different inputs files.
- `Total.R_input.txt`: ratio of the occurrence of major metabolic pathways in the all communities.

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

### Function occurrence and abundance

For visualisation, two values are used to represent the functions.
First, the **occurrence** corresponding to the number of organisms having this function dividing by the total number of organisms in the community. If you give an `abundance file`, a second value is used, the **abundance** (computed for each sample in the abundance file). The abundance of a function is the sum of the abundance of organisms having it divided by the sum of abundance of all organisms in the sample.

For example, if we look at the function `Formate oxidation fdoG` in a community. If 20 organisms in this community have this function on a community having a total of 80 organisms, the **occurrence** of this function is 0.25 (20 / 80). Then, let's say that these 20 organisms ahve a summed abundance of 600 and the total abundance of all organisms in the community is 1200, then the **abundance** of the function is 0.5 (600 / 1200).

### Output

Several output are created by bigecyhmm_visualisation.

````
output_folder
├── function_abundance
│   └── function_participation
│       └── sample_1.tsv
│       └── ...
│   └── heatmap_abundance_samples.png
│   └── heatmap_abundance_samples.tsv
│   └── polar_plot_merged.png
├── function_occurrence
│   └── heatmap_occurrence.png
│   └── heatmap_occurrence.tsv
│   └── polar_plot_merged.png
├── bigecyhmm_visualisation.log
├── bigecyhmm_visualisation_metadata.json
````

`function_abundance` is a folder containing all visualisation associated with abundance values. It contains:

- `function_participation`: a folder containing one tabulated file per sample from the abundance file. For each sample, it gives the function abundance associated with each organism in the community.

- `heatmap_abundance_samples.png`: a heatmap showing the abundance for all the HMMs searched by bigecyhmm in the different samples.
- `heatmap_abundance_samples.tsv`: the tabulated file associated with the creation of the `heatmap_abundance_samples.png` file.
- `polar_plot_merged.png`: a polar plot showing the abundance of major functions in the samples.

`function_occurrence` is a folder containing all visualisation associated with occurrence values. It contains:

- `heatmap_occurrence.png`: a heatmap showing the occurrence for all the HMMs searched by bigecyhmm in the community (all the input protein files).
- `heatmap_occurrence.tsv`: the tabulated file associated with the creation of the `heatmap_occurrence.png` file.
- `polar_plot_merged.png`: a polar plot showing the occurrence of major functions in the samples.

`bigecyhmm_visualisation.log` is a log file.

`bigecyhmm_visualisation_metadata.json` is a metadata file giving information on the version of the package used.

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
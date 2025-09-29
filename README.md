[![PyPI version](https://img.shields.io/pypi/v/bigecyhmm.svg)](https://pypi.org/project/bigecyhmm/) [![](pictures/doi_tabigecy.svg)](https://doi.org/10.1093/bioinformatics/btaf230)

# bigecyhmm: Biogeochemical cycle HMMs search

Bigecyhmm is a Python package to search for genes associated with biogeochemical cycles in protein sequence fasta files. It begins as a self-contained, lightweight reimplementation of a subtask performed in [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) but has since grown. Bigecyhmm default behaviour searches for enzymes associated with carbon, sulfur, nitrogen and phosphorus cycles using HMMs from METABOLIC article, KEGG, PFAM, TIGR. It can be also used with a custom database and then will output network representation of the cycle.

## 0 Table of contents
- [bigecyhmm: Biogeochemical cycle HMMs search](#bigecyhmm-biogeochemical-cycle-hmms-search)
  - [0 Table of contents](#0-table-of-contents)
  - [1 Dependencies](#1-dependencies)
  - [2 Installation](#2-installation)
  - [3 bigecyhmm](#3-bigecyhmm)
    - [3.1 Usage](#31-usage)
    - [3.2 Output](#32-output)
  - [4 bigecyhmm\_visualisation](#4-bigecyhmm_visualisation)
    - [4.1 Function occurrence and abundance](#41-function-occurrence-and-abundance)
    - [4.2 Output of bigecyhmm\_visualisation](#42-output-of-bigecyhmm_visualisation)
  - [5 Custom usage](#5-custom-usage)
    - [5.1 Bigecyhmm internal database](#51-bigecyhmm-internal-database)
      - [5.1.1 Contribution to bigecyhmm internal database](#511-contribution-to-bigecyhmm-internal-database)
      - [5.1.2 Modifying bigecyhmm internal database](#512-modifying-bigecyhmm-internal-database)
    - [5.2 `bigecyhmm_custom`: using custom database](#52-bigecyhmm_custom-using-custom-database)
      - [5.2.1 Requirements](#521-requirements)
      - [5.2.2 Inputs](#522-inputs)
      - [5.2.3 Outputs](#523-outputs)
  - [6 Citation](#6-citation)

## 1 Dependencies

bigecyhmm is developed to be as minimalist as possible. It requires:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [Pillow](https://github.com/python-pillow/Pillow): to create biogeochemical cycle diagrams.

The HMMs used are stored inside the package as a zip file ([hmm_files.zip](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/bigecyhmm/hmm_databases)). It makes this python package a little heavy (around 19 Mb) but in this way, you do not have to download other files and can directly use it.

For `bigecyhmm_visualisation`, you also needs the following packages:

- [pandas](https://pypi.org/project/pandas/): to read the input files.
- [seaborn](https://github.com/mwaskom/seaborn) and [matplotlib](https://github.com/matplotlib/matplotlib): to create most of the figures.

For `bigecyhmm_custom`, you also needs the following package:

- [networkx](https://github.com/networkx/networkx): to handle custom biogeochemical cycle as a graph.
- [matplotlib](https://github.com/matplotlib/matplotlib): to create automatically (bad) visualisation of the cycle.

## 2 Installation

It can be installed from PyPI:

`pip install bigecyhmm`

Or it can be installed with pip by cloning the repository:

```sh
git clone https://github.com/ArnaudBelcour/bigecyhmm.git

cd bigecyhmm

pip install -e .
```

For `bigecyhmm_visualisation`, you also needs to run:

`pip install pandas seaborn`

For `bigecyhmm_custom`, you also needs to run:

`pip install networkx matplotlib`

## 3 bigecyhmm

### 3.1 Usage

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

* `-c` to indicate the number of core used. It is only useful if you have multiple protein fasta files as the added cores will be used to run another HMM search on a different protein fasta file.

### 3.2 Output

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

## 4 bigecyhmm_visualisation

There is a second command associated with bigecyhmm (`bigecyhmm_visualisation`), to create visualisation of the results.

To create the associated figures, there are other dependencies:

- [pandas](https://pypi.org/project/pandas/): to read the input files.
- [seaborn](https://github.com/mwaskom/seaborn) and [matplotlib](https://github.com/matplotlib/matplotlib): to create most of the figures.

Two subcommands are available for `bigecyhmm_visualisation`:

- `bigecyhmm_visualisation esmecata`: to create visualisation from EsMeCaTa and bigecyhmm outputs folder (with optionally an abundance file).
- `bigecyhmm_visualisation genomes`: to create visualisation from bigecyhmm output folder (with optionally an abundance file).

There are four parameters:

- `--esmecata`: EsMeCaTa output folder associated with the run (as the visualisation works on esmecata results). Only required for `bigecyhmm_visualisation esmecata`.
- `--bigecyhmm`: bigecyhmm output folder associated with the run. Required for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `--abundance-file`: abundance file indicating the abundance for each organisms selected by EsMeCaTa. Optional for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `-o`: an output folder. Required for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.

### 4.1 Function occurrence and abundance

For visualisation, two values are used to represent the functions.
First, the **occurrence** corresponding to the number of organisms having this function dividing by the total number of organisms in the community. If you give an `abundance file`, a second value is used, the **abundance** (computed for each sample in the abundance file). The abundance of a function is the sum of the abundance of organisms having it divided by the sum of abundance of all organisms in the sample.

For example, if we look at the function `Formate oxidation fdoG` in a community. If 20 organisms in this community have this function on a community having a total of 80 organisms, the **occurrence** of this function is 0.25 (20 / 80). Then, let's say that these 20 organisms have a summed abundance of 600 and the total abundance of all organisms in the community is 1200, then the **abundance** of the function is 0.5 (600 / 1200).

### 4.2 Output of bigecyhmm_visualisation

Several output are created by bigecyhmm_visualisation.

```
output_folder
├── function_abundance
│   ├── cycle_diagrams_abundance
│   |   └── sample_1_carbon_cycle.png
│   |   └── sample_1_nitrogen_cycle.png
│   |   └── ...
│   ├── function_participation
│   |   └── sample_1.tsv
│   |   └── ...
│   ├── cycle_participation
│   |   └── sample_1.tsv
│   |   └── ...
│   └── barplot_esmecata_found_taxon_sample.png
│   └── barplot_esmecata_found_organism_sample.tsv
│   └── cycle_abundance_sample.tsv
│   └── cycle_abundance_sample_melted.tsv
│   └── cycle_abundance_sample_raw.tsv
│   └── function_abundance_sample.tsv
│   └── heatmap_abundance_samples.png
│   └── hmm_functional_profile.tsv
│   └── polar_plot_abundance_sample_1.png
│   └── polar_plot_abundance_sample_XXX.png
├── function_occurrence
│   └── cycle_occurence.tsv
│   └── diagram_carbon_cycle.png
│   └── diagram_nitrogen_cycle.png
│   └── diagram_sulfur_cycle.png
│   └── diagram_other_cycle.png
│   └── function_occurrence.tsv
│   └── function_occurrence_in_organism.tsv
│   └── heatmap_occurrence.png
│   └── pathway_presence_in_organism.tsv
├── bigecyhmm_visualisation.log
├── bigecyhmm_visualisation_metadata.json
```

`function_abundance` is a folder containing all visualisation associated with abundance values. It contains:

- `cycle_diagrams_abundance`: a folder containing 4 cycle diagrams (carbon, sulfur, nitrogen and other) from METABOLIC per sample from the abundance file. For each sample, it gives the abundance and the relative abundance of the major function.
- `function_participation`: a folder containing one tabulated file per sample from the abundance file. For each sample, it gives the function abundance associated with each organism in the community.
- `cycle_participation`: a folder containing one tabulated file per sample from the abundance file. For each sample, it gives the cycle abundance associated with each organism in the community.
- `barplot_esmecata_found_taxon_sample.png`: a barplot displaying the coverage of EsMeCaTa according to the abundances from samples. Each bar corresponds to a sample, the y-axis shows the relative abundances of the organisms in the sample. The color indicates which taxonomic rank has been used by EsMeCaTa to predict the consensus proteomes. If EsMeCaTa was not able to predict a consensus proteomes, it is displayed in category `Not found`. With this figure, you can have an idea if there is enough predictions for the different samples in the dataset and at which taxonomic ranks these predictiosn have been made. Thus allowing the estimation of the quality of the predictions: predictions are better if they are closer to lower taxonomic ranks (genus family). `barplot_esmecata_found_organism_sample.tsv` is the input file used to create the figure.
- `function_abundance_sample.tsv`: a tabulated file containing the relative abundance of each function according to the abundance of the associated organisms in the different sample. Rows correspond to the functions and columns correspond to the samples. It is used to create the `heatmap_abundance_samples.png` file. The file `hmm_functional_profile.tsv` contains the absolute abundance of the functions.
- `heatmap_abundance_samples.png`: a heatmap showing the abundance for all the HMMs searched by bigecyhmm in the different samples.
- `cycle_abundance_sample.tsv`: a tabulated file showing the relative abundance of major functions in biogeochemical cycles according to the organisms. Rows correspond to the major functions and columns correspond to the samples. The file `cycle_abundance_sample_melted.tsv` is a melted version of this file. The file `cycle_abundance_sample_raw.tsv` contains the absolute abundance of the functions.
- `polar_plot_abundance_samples_XXXX.png`: a polar plot showing the abundance of major functions in the sample `XXXX`.

`function_occurrence` is a folder containing all visualisation associated with occurrence values. It contains:

- `cycle_occurence.tsv`: a tabulated file showing the occurrence of major functions in biogeochemical cycles. Rows correspond to the major function and the column corresponds to the community.
- `diagram_*.png`: diagram representing a biogeochemical cycles (carbon, nitrogen, sulfur, other) from METABOLIC. It shows the number of organisms with predicted major functions and the relative occurrence of these functions.
- `function_occurrence.tsv`: a tabulated file containing the ratio for each function. Rows correspond to the function and the column corresponds to the community. It is used to create the `heatmap_occurrence.png` file.
- `function_occurrence_in_organism.tsv`: a tabulated file containing the occurrence of function in each organism of the samples.
- `heatmap_occurrence.png`: a heatmap showing the occurrence for all the HMMs searched by bigecyhmm in the community (all the input protein files).
- `pathway_presence_in_organism.tsv`: a tabulated file containing the occurrence of cycle funcitons in each organism of the samples.
- `swarmplot_function_ratio_community.png`: a swarmplot showing the occurrence of major functions in the samples.

`bigecyhmm_visualisation.log` is a log file.

`bigecyhmm_visualisation_metadata.json` is a metadata file giving information on the version of the package used.

## 5 Custom usage

### 5.1 Bigecyhmm internal database

#### 5.1.1 Contribution to bigecyhmm internal database

If you are interested in specific functions associated with cycles present in bigecyhmm (carbon, sulfur, nitrogen, phosphorus) and want to propose an addition, you can create an issue or a Pull Request.
Depending on the additions or modifications, it will be taken into account. Keep in mind that bigecyhmm's goal is to limit itself to a small internal database.
If you want to completely add another cycle, please refer to the next subsection.

#### 5.1.2 Modifying bigecyhmm internal database

You can also edit the database to add your own functions. To do so, you can either clone this repository or make a fork. Then install bigecyhmm using `pip install -e .` inside bigecyhmm folder (where the file `pyproject.toml` is located). You can modify the internal database in different ways. There are files to potentially modify:

- `hmm_databases/hmm_files.zip`: a zip file containing the HMM files used to screen the associated genes. Uncompress it, add the HMM files you want in it and then compress it.
- `hmm_databases/hmm_table_template.tsv`: a tabulated file containing the association between functions and HMMs. For each HMM you add, you have to add a line in this file. There are two mandatory columns (1) `Hmm file` (name of the HMM file in `hmm_files.zip`) and (2) `Hmm detecting threshold` (threshold used to filter matches).
- `hmm_databases/cycle_pathways.tsv`: a tabulated file linking major functions to HMMs. This file is linked to the creation of the diagrams. If you want to modify this file and propagate the change to the diagram, you must (1) edit diagram templates located at `templates/*` and (2) edit `diagram_cycles.py`, especially function called `create_carbon_cycle` (and the one for the other cycles). In this function several lines are associated with the major function: `data_step_01 = diagram_data['C-S-01:Organic carbon oxidation']` extracts function abundance from predictions, `imgdraw.text((800,80), 'Step1: Organic carbon\n oxidation\n{0}: {1}\n{2}: {3}%'.format(first_term, data_step_01[0], second_term, data_step_01[1]), (0,0,0), font=font)` puts the prediction on the template. Modifying the template requires to also modifies these scripts.

**Adding new HMM**

If you want to add a new HMM for the search, just modify `hmm_databases/hmm_files.zip` and `hmm_databases/hmm_table_template.tsv`.

**Adding new pathway/diagram**

If you want to modify or create a diagram: you have to put the new associated HMMs in `hmm_databases/hmm_files.zip` and `hmm_databases/hmm_table_template.tsv`. Then modify `hmm_databases/cycle_pathways.tsv` by adding the new pathways associated with HMMs.

If you have multiple HMMs for the same pathway, you can separate them with a `, `. If you have two HMMs that are required at the same time, you have to separated them with a `; `. For example, `soxZ.hmm, soxA.hmm; soxC.hmm, soxD.hmm` means that to have the associated function you must have either *soxZ* or *soxZ* **AND** either *soxC* or *soxD*.

It is also possible to say that a function should not be associated with a HMM by prefixing `NO|` to the HMM filename. For example, `soxZ.hmm, soxA.hmm; NO|soxC.hmm, NO|soxD.hmm` means that to have the associated function you must have either *soxZ* or *soxZ* **AND** **NOT** *soxC* or *soxD*.

Then you also have to modify the diagram template in `templates/*` (you can modify the svg and then extract the new template in png). Finally, you will have to modify `diagram_cycles.py` to correctly place the new predictions on the template.

To do so, you have to change the coordinates of the text in `diagram_cycles.py`. For example, in the line

```python
  data_step_01 = diagram_data['C-S-01:Organic carbon oxidation']
  imgdraw.text((800,80), 'Step1: Organic carbon\n oxidation\n{0}: {1}\n{2}: {3}%'.format(first_term, data_step_01[0], second_term, data_step_01[1]), (0,0,0), font=font)
```

`(800,80)` corresponds to the coordinates of the text on the figure, by adjusting it you can move the text. First number is associated with x-axis and second number is associated with y-axis. For x-axis, 0 begins at the left of the figure with higher numbers going towards the right. For y-axis, 0 begins at the top of the figure with higher numbers going towards the bottom. `(0,0,0)` corresponds to the color of the text.


### 5.2 `bigecyhmm_custom`: using custom database

**Warning**: This is a prototype.

It is possible to create a completely custom database that is linked to a specific biogeochemical cycles (or metabolic networks) using `bigecyhmm_custom`.

#### 5.2.1 Requirements

This command requires three packages:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [networkx](https://github.com/networkx/networkx): to handle biogeochemical cycle as a graph.
- [matplotlib](https://github.com/matplotlib/matplotlib): to create automatically (bad) visualisation of the cycle.

#### 5.2.2 Inputs

This command line expects two arguments:

- `-i`: an input protein sequence fasta file/folder.
- `-d`: a file/folder containing the custom databases. `bigecyhmm_custom` will iterate other the file/folder to search for every `.json` files. If it finds one, it will search for associated `.tsv` and `.zip` (files with the same name and at the same location but with either a tsv or zip extension). The three expected files are listed below:
  - a `json` file representing the biogeochemical cycle as a bipartite graph with nodes representing `metabolite` and `function`. Example can be found in the test folder, such as [carbon cycle json file](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/custom_db/carbon_cycle.json). The `hmm` field in the `function node` in the json is mandatory to indicate the HMMs associated with the functions of the cycle. The HMMs are represented as a string with `, ` separating HMMs as a `OR` relation (meaning these HMMs are redundant) and `; ` as a `AND` relation (meaning that both HMMs are required).
  - a `zip` file containing the HMM profiles (`.hmm` files) such as the one used by bigecyhmm ([hmm_files.zip](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/hmm_files.zip)). If no file is present in the folder, bigecyhmm will use its internal HMM database. You can search for HMM in [KEGG Ortholog database](https://www.genome.jp/kegg/ko.html), [Protein Family Models from NIH](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/), [PFAM](https://www.ebi.ac.uk/interpro/download/Pfam/) or [EggNOG](http://eggnog5.embl.de/#/app/home). It is also possible to build them, an example can be found with [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/examples/msa_to_hmm.html#Build-an-HMM-from-a-multiple-sequence-alignment).
  - a `tsv` file containing the threshold for the different HMMs. If no file is present in the folder, bigecyhmm will use its internal template file for threshold. An example can be found in bigecyhmm internal database ([hmm_table_template.tsv](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/hmm_table_template.tsv)) or in the test folder ([hmm_table_template.tsv](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/mini_custom_db/hmm_table_template.tsv)).

An example with mini database is present in the [test folder](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/test/input_data/mini_custom_db).

Here are several examples of inputs:

- Only a json file, bigecyhmm will use its internal HMM database to search for HMM files from the json file (associated argument `-d custom_db_cycle.json`):
```
custom_db_cycle.json
```

- A folder with one json file and tsv/zip files (associated argument `-d custom_db_cycle`):
```
custom_db_cycle
├── custom_db_cycle.json
├── custom_db_cycle.tsv
├── custom_db_cycle.zip
```

- A folder with several json files (associated argument `-d custom_db_cycle`):
```
custom_db_cycle
├── carbon_cycle.json
├── carbon_cycle.tsv
├── carbon_cycle.zip
├── nitrogen_cycle.json
├── nitrogen_cycle.tsv
├── nitrogen_cycle.zip
├── sulfur_cycle.json
├── sulfur_cycle.tsv
├── sulfur_cycle.zip
```

Usage example:
```
bigecyhmm_custom -i protein_sequences.faa -d custom_db -o output_folder
```

It can take five optional arguments:

- `--abundance-file`: an abundance file containing the abundance of the organisms associated with the protein sequences given as input in different samples.
- `--measure-file`: a measurement file containing the measures of metabolites of the biogeochemical cycle in different samples.
- `--esmecata`: by giving an esmecata output folder, `bigecyhmm_custom` maps taxon_id to organism names to associate organism abundance with esmecata predicitons.
- `-m`: JSON file containing gene associated with protein motifs to check for predictions. This verification comes from the [METABOLIC article](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8#Sec2) (you can find information about it, in the section `Motif validation`). The protein motif corresponds to a regex associated with amnio-acids or `X` (the latter being any amino-acid). The idea of this verification is to check if an expected amino-acid motif is present in the sequence matching the associated HMM. You can see an example file in the test folder ([motif.json](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/motif.json)). The name of the gene corresponds to the name of its HMM. If no file is given, it will be using the default ones from METABOLIC (you can find it [here](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/__init__.py#L39) as a dicitonary).
- `-p`: JSON file containing association between two genes to check for predictions. This verification comes from the [METABOLIC article](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8#Sec2) (you can find information about it, in the section `Motif validation`). It ensures that a sequence is properly associated with a specific HMM and not to anotehr yet similar HMM. An example file can be found in the test folfer ([motif_pair.json](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/motif_pair.json)). It contains association between two gene names. The HMM search results of the sequence against these two gnee profiles are compared to find the one with a better score. The name of the gene corresponds to the name of its HMM. If no file is given, it will be using the default ones from METABOLIC (you can find it [here](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/__init__.py#L50) as a dicitonary).

#### 5.2.3 Outputs

`bigecyhmm_custom` creates inside the output folder on folder per input custom json file.
It outputs similar files than bigecyhmm classical output except for the cycle visualisation.
As it is more difficult to provide a direct visualisation from a custom database `bigecyhmm_custom` relies on network representation to create these visualisations.
To do so, it creates network file (`cycle_diagram_bipartite_occurrence.graphml`) as output. These files can be used in network software (such as [Cytoscape](https://cytoscape.org/), or [pyvis](https://github.com/WestHealth/pyvis)) to generate visualisation. It also tries to create a visualisation with networkx and matpltolib but they are not very good.

If you have given abundance and measure files, a second network file (`cycle_diagram_bipartite_abundance.graphml`) is created where function nodes are associated with the summed abundance of organisms in the different samples and metabolite node are associated with their measures in the different samples.

## 6 Citation

If you have used bigecyhmm in an article, please cite:

Arnaud Belcour, Loris Megy, Sylvain Stephant, Caroline Michel, Sétareh Rad, Petra Bombach, Nicole Dopffel, Hidde de Jong and Delphine Ropers. Predicting coarse-grained representations of biogeochemical cycles from metabarcoding data *Bioinformatics*, Volume 41, Issue Supplement_1, July 2025, Pages i49–i57, https://doi.org/10.1093/bioinformatics/btaf230

Bigecyhmm relies on:

- PyHMMER for the search on the HMMs:

Martin Larralde and Georg Zeller. PyHMMER: a python library binding to HMMER for efficient sequence analysis. Bioinformatics, 39(5):btad214, 2023.  https://doi.org/10.1093/bioinformatics/btad214

- HMMer website for the search on the HMMs:

HMMER. http://hmmer.org. Accessed: 2022-10-19.

- the following articles for the creation of the custom HMMs:

Zhou, Z., Tran, P.Q., Breister, A.M. et al. METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks. Microbiome 10, 33, 2022. https://doi.org/10.1186/s40168-021-01213-8

Anantharaman, K., Brown, C., Hug, L. et al. Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system. Nat Commun 7, 13219, 2016. https://doi.org/10.1038/ncomms13219

- the following article for KOfam HMMs:

Takuya Aramaki, Romain Blanc-Mathieu, Hisashi Endo, Koichi Ohkubo, Minoru Kanehisa, Susumu Goto, Hiroyuki Ogata, KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold, Bioinformatics, Volume 36, Issue 7, 2020, Pages 2251–2252, https://doi.org/10.1093/bioinformatics/btz859

- the following article for TIGRfam HMMs:

Jeremy D. Selengut, Daniel H. Haft, Tanja Davidsen, Anurhada Ganapathy, Michelle Gwinn-Giglio, William C. Nelson, Alexander R. Richter, Owen White, TIGRFAMs and Genome Properties: tools for the assignment of molecular function and biological process in prokaryotic genomes, Nucleic Acids Research, Volume 35, Issue suppl_1, 2007, Pages D260–D264, https://doi.org/10.1093/nar/gkl1043

- the following article for Pfam HMMs:

Robert D. Finn, Alex Bateman, Jody Clements, Penelope Coggill, Ruth Y. Eberhardt, Sean R. Eddy, Andreas Heger, Kirstie Hetherington, Liisa Holm, Jaina Mistry, Erik L. L. Sonnhammer, John Tate, Marco Punta, Pfam: the protein families database, Nucleic Acids Research, Volume 42, Issue D1, 2014, Pages D222–D230, https://doi.org/10.1093/nar/gkt1223

- the following articles for phosphorus cycle:

Boden, J.S., Zhong, J., Anderson, R.E. et al. Timing the evolution of phosphorus-cycling enzymes through geological time using phylogenomics. Nature Communications, 15, 3703 (2024). https://doi.org/10.1038/s41467-024-47914-0

Siles, J. A., Starke, R., Martinovic, T., Fernandes, M. L. P., Orgiazzi, A., & Bastida, F. Distribution of phosphorus cycling genes across land uses and microbial taxonomic groups based on metagenome and genome mining. Soil Biology and Biochemistry, 174, 108826, 2022. https://doi.org/10.1016/j.soilbio.2022.108826


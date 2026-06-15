[![PyPI version](https://img.shields.io/pypi/v/bigecyhmm.svg)](https://pypi.org/project/bigecyhmm/) [![](pictures/doi_tabigecy.svg)](https://doi.org/10.1093/bioinformatics/btaf230)

# bigecyhmm: Biogeochemical cycle HMMs search

Bigecyhmm is a Python package to search for genes associated with biogeochemical cycles in protein sequence fasta files. It begins as a self-contained, lightweight reimplementation of a subtask performed in [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) but has since grown. Bigecyhmm, by default, searches for enzymes associated with carbon, sulfur, nitrogen and phosphorus cycles using HMMs from METABOLIC article, KEGG, PFAM, TIGR. It can be also used with a custom database and then will output network representation of the cycle.

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
        - [5.2.2.1 Inputs: custom tsv file](#5221-inputs-custom-tsv-file)
        - [5.2.2.2 Inputs: custom json file](#5222-inputs-custom-json-file)
        - [5.2.2.3 Inputs: folder](#5223-inputs-folder)
        - [5.2.2.3 Inputs: internal custom database](#5223-inputs-internal-custom-database)
      - [5.2.3 Usage](#523-usage)
      - [5.2.4 Outputs](#524-outputs)
  - [6 Citation](#6-citation)
  - [License](#license)

## 1 Dependencies

`bigecyhmm` is developed to be as minimalist as possible. It requires:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [Pillow](https://github.com/python-pillow/Pillow): to create biogeochemical cycle diagrams.

The HMMs used are stored inside the package as a folder ([hmm_files](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/bigecyhmm/hmm_databases)). It makes this python package a little heavy (around 19 Mb when compressed) but in this way, you do not have to download other files and can directly use it.

For `bigecyhmm_custom` (use custom database instead of the default one), you also needs the following package:

- [networkx](https://github.com/networkx/networkx): to handle custom biogeochemical cycle as a graph.
- [pygraphviz](https://github.com/pygraphviz/pygraphviz): to render layout of bipartite graph.
- [matplotlib](https://github.com/matplotlib/matplotlib): to create automatically (bad) visualisation of the cycle.

For `bigecyhmm_visualisation` (linking metabolic predictions to organism abundance), you also needs the following packages:

- [pandas](https://pypi.org/project/pandas/): to read the input files.
- [seaborn](https://github.com/mwaskom/seaborn) and [matplotlib](https://github.com/matplotlib/matplotlib): to create most of the figures.
- [scipy](https://github.com/scipy/scipy) and [statsmodels](https://github.com/statsmodels/statsmodels/) for statistical tests when giving a group file.
- [esmecata](github.com/AuReMe/esmecata) to handle esmecata results.
- If use with results form `bigecyhmm_custom`:
   - [networkx](https://github.com/networkx/networkx): to handle custom biogeochemical cycle as a graph.
   - [pygraphviz](https://github.com/pygraphviz/pygraphviz): to render layout of bipartite graph.


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

`pip install pandas seaborn esmecata scipy statsmodels esmecata`

For `bigecyhmm_custom`, you also needs to run:

`pip install networkx matplotlib pygraphviz`

## 3 bigecyhmm

### 3.1 Usage

Bigecyhmm can be used on two types of input:

- by giving a protein fasta file:

```sh
bigecyhmm -i protein_sequence.faa -o output_dir
```

- by giving a folder containing multiple fasta files:

```
protein_sequences_folder
├── protein_sequences_1.faa
├── protein_sequences_2.faa
├── ...
```

```sh
bigecyhmm -i protein_sequences_folder -o output_dir
```

There is one option:

* `-c` to indicate the number of core used. It is only useful if you have multiple protein fasta files as the added cores will be used to run another HMM search on a different protein fasta file.

### 3.2 Output

It gives as output:

```
output_folder
├── diagram_figures
│   ├── carbon_cycle.png
│   ├── nitrogen_cycle.png
│   ├── other_cycle.png
│   ├── phosphorus_cycle.png
│   ├── sulfur_cycle.png
├── diagram_input
│   ├── org_1.txt
│   ├── org_2.txt
│   ├── org_3.txt
├── hmm_results
│   ├── org_1.tsv
│   ├── org_2.tsv
│   ├── org_3.tsv
├── bigecyhmm.log
├── bigecyhmm_metadata.json
├── function_presence.tsv
├── mapping_pathway_to_function_name.tsv
├── pathway_presence.tsv
├── pathway_presence_hmms.tsv
├── R_input.txt
```

- a folder `hmm_results`: one tsv files showing the hits for each protein fasta file.
- `function_presence.tsv` a tsv file showing the presence/absence of generic functions associated with the HMMs that matched.
- a folder `diagram_input`, the necessary input to create Carbon, Nitrogen, Sulfur and other cycles with the [R script](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/scripts/draw_biogeochemical_cycles.R) modified from the [METABOLIC repository](https://github.com/AnantharamanLab/METABOLIC) using the following command: `Rscript draw_biogeochemical_cycles.R bigecyhmm_output_folder/diagram_input_folder/ diagram_output TRUE`. This script requires the diagram package that could be installed in R with `install.packages('diagram')`.
- a folder `diagram_figures` contains biogeochemical diagram figures drawn from template situated in `bigecyhmm/templates`.
- `bigecyhmm.log`: log file.
- `bigecyhmm_metadata.json`: bigecyhmm metadata (Python version used, package version used).
- `function_presence.tsv`: occurrence of the functions in the different input protein files.
- `mapping_pathway_to_function_name.tsv`: linking pathway name to more specific function name.
- `pathway_presence.tsv`: occurrence of the major metabolic pathways in the different inputs files.
- `pathway_presence_hmms.tsv`: HMMs with matches for the major metabolic pathways in the different inputs files.
- `Total.R_input.txt`: ratio of the occurrence of major metabolic pathways in the all communities.

## 4 bigecyhmm_visualisation

There is a second command associated with bigecyhmm (`bigecyhmm_visualisation`), to create visualisation of the results.

To create the associated figures, there are other dependencies:

- [pandas](https://pypi.org/project/pandas/): to read the input files.
- [seaborn](https://github.com/mwaskom/seaborn) and [matplotlib](https://github.com/matplotlib/matplotlib): to create most of the figures.
- [networkx](https://github.com/networkx/networkx): to handle biogeochemical cycle as a graph if hanlding results from `bigecyhmm_custom`.
- [pygraphviz](https://github.com/pygraphviz/pygraphviz): to render layout of bipartite graph.
- [scipy](https://github.com/scipy/scipy) and [statsmodels](https://github.com/statsmodels/statsmodels/) for statistical tests when giving a group file.

Two subcommands are available for `bigecyhmm_visualisation`:

- `bigecyhmm_visualisation esmecata`: to create visualisation from EsMeCaTa and bigecyhmm outputs folder (with optionally an abundance file).
- `bigecyhmm_visualisation genomes`: to create visualisation from bigecyhmm output folder (with optionally an abundance file).

There are several parameters:

- `--esmecata`: EsMeCaTa output folder associated with the run (as the visualisation works on EsMeCaTa results). Only required for `bigecyhmm_visualisation esmecata`.
- `--bigecyhmm`: bigecyhmm output folder associated with the run. Required for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `-o`: an output folder. Required for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `--abundance-file`: abundance file indicating the abundance for each organisms selected by EsMeCaTa. Optional for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `--measure-file`: abundance file indicating the abundance for each metabolites (for bipartite graph). Optional for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `--group-file`: tabulated file indicating the group for each sample. Optional for both `bigecyhmm_visualisation esmecata` and `bigecyhmm_visualisation genomes`.
- `--background-file`: background image used for the donut plot. If no background images are given, bigecyhmm uses (1) a template one or (2) generate a bipartite graph representation (if results are from `bigecyhmm_custom`). Optional for `bigecyhmm_visualisation esmecata`.

`--group-file` expects a tabulated file like this (you have [an example](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/tests/input_data/group_sample.tsv) in test folder):

| sample   | group |
|----------|-------|
| sample_1 | A     |
| sample_2 | A     |
| sample_3 | B     |

It will be used to generate donut plot with different colours associated with each group (and then computes statistics).

### 4.1 Function occurrence and abundance

For visualisation, two values are used to represent the functions.
First, the **occurrence** corresponding to the number of organisms having this function dividing by the total number of organisms in the community. If you give an `abundance file`, a second value is used, the **abundance** (computed for each sample in the abundance file). The abundance of a function is the sum of the abundance of organisms having it divided by the sum of abundance of all organisms in the sample.

For example, if we look at the function `Formate oxidation fdoG` in a community. If 20 organisms in this community have this function on a community having a total of 80 organisms, the **occurrence** of this function is 0.25 (20 / 80). Then, let's say that these 20 organisms have a summed abundance of 600 and the total abundance of all organisms in the community is 1200, then the **abundance** of the function is 0.5 (600 / 1200).

### 4.2 Output of bigecyhmm_visualisation

Several output are created by bigecyhmm_visualisation:

```
output_folder
├── function_abundance
│   ├── cycle_diagrams_abundance
│   |   └── sample_1_carbon_cycle.png
│   |   └── sample_1_nitrogen_cycle.png
│   |   └── ...
│   ├── cycle_participation
│   |   └── sample_1.tsv
│   |   └── ...
│   ├── cycle_taxa_abundance
│   |   └── function_1.png
│   |   └── function_1.tsv
│   |   └── ...
│   ├── function_participation
│   |   └── sample_1.tsv
│   |   └── ...
│   ├── plot_donut
│   |   └── cleaned_data.tsv
│   |   └── combined_donut_table.png
│   |   └── group_medians_donut.png
│   |   └── group_stats_table.png
│   |   └── group_stats.tsv
│   ├── plot_donut_graph (if running bigecyhmm_custom)
│   |   └── ...
│   ├── polar_plot_abundance
│   |   └── sample_1.png
│   |   └── ...
│   └── barplot_esmecata_found_taxon_sample.png
│   └── barplot_esmecata_found_organism_sample.tsv
│   └── barplot_esmecata_missing_organism_sample.tsv
│   └── cycle_abundance_sample.tsv
│   └── cycle_abundance_sample_melted.tsv
│   └── cycle_abundance_sample_raw.tsv
│   └── cycle_pathways_bubble_plot.png
│   └── function_abundance_sample.tsv
│   └── function_abundance_sample_raw.tsv
│   └── heatmap_abundance_samples.png
│   └── heatmap_abundance_samples.svg
│   └── hmm_functional_profile.tsv
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

If you have used `bigecyhmm_custom`, the output folder will be a little different as bigecyhmm will create one subfolder for each custom database used and will not generate `cycle_diagrams_abundance` folder:

```
output_folder
├── custom_database_1
│   ├── function_abundance
│   ├── function_occurrence
│   ├── bigecyhmm_visualisation_metadata.json
│   ├── ...
├── custom_database_2
│   ├── function_abundance
│   ├── function_occurrence
│   ├── bigecyhmm_visualisation_metadata.json
│   ├── ...
├── bigecyhmm_visualisation.log
```

`function_abundance` is a folder containing all visualisation associated with abundance values. It contains:

- `cycle_diagrams_abundance`: a folder containing 4 cycle diagrams (carbon, sulfur, nitrogen and other) from METABOLIC per sample from the abundance file. For each sample, it gives the abundance and the relative abundance of the major function.
- `cycle_participation`: a folder containing one tabulated file per sample from the abundance file. For each sample, it gives the cycle abundance associated with each organism in the community.
- `cycle_taxa_abundance`: a folder containing one png and one tabulated file per metabolic function. They show the organims abundance identified to possess the associated metabolic functions in the samples.
- `function_participation`: a folder containing one tabulated file per sample from the abundance file. For each sample, it gives the function abundance associated with each organism in the community.
- `plot_donut`: a folder containing a donut plot showing metabolic abundance in samples according to group and results of statistical tests (Kruskal-Wallis with Benjamini-Hochberg correction). A background image is used to represent the studied metabolism. 
- `plot_donut_graph` (only with `bigecyhmm_custom`): a folder containing a donut plot showing metabolic abundance in samples according to group and results of statistical tests (Kruskal-Wallis with Benjamini-Hochberg correction). A bipartite graph image is generated as a background image showing a representation of the studied metabolism.
- `polar_plot_abundance`: a folder containing polar plots showing the abundance of major functions in each sample.
- `barplot_esmecata_found_taxon_sample.png`: a barplot displaying the coverage of EsMeCaTa according to the abundances from samples. Each bar corresponds to a sample, the y-axis shows the relative abundances of the organisms in the sample. The color indicates which taxonomic rank has been used by EsMeCaTa to predict the consensus proteomes. If EsMeCaTa was not able to predict a consensus proteomes, it is displayed in category `Not found`. With this figure, you can have an idea if there is enough predictions for the different samples in the dataset and at which taxonomic ranks these predictions have been made. Thus allowing the estimation of the quality of the predictions: predictions are better if they are closer to lower taxonomic ranks (genus family). `barplot_esmecata_found_organism_sample.tsv` is the input file used to create the figure.
- `barplot_esmecata_missing_organism_sample.tsv`: a tabulated file showing organisms without predictions from EsMeCaTa.
- `function_abundance_sample.tsv`: a tabulated file containing the relative abundance of each function according to the abundance of the associated organisms in the different sample. Rows correspond to the functions and columns correspond to the samples. It is used to create the `heatmap_abundance_samples.png` file. The file `function_abundance_sample_raw.tsv` contains the absolute abundance of the functions.
- `heatmap_abundance_samples.png`: a heatmap showing the abundance for all the HMMs searched by bigecyhmm in the different samples.
- `cycle_abundance_sample.tsv`: a tabulated file showing the relative abundance of major functions in biogeochemical cycles according to the organisms. Rows correspond to the major functions and columns correspond to the samples. The file `cycle_abundance_sample_melted.tsv` is a melted version of this file. The file `cycle_abundance_sample_raw.tsv` contains the absolute abundance of the functions.

`function_occurrence` is a folder containing all visualisation associated with occurrence values. It contains:

- `cycle_occurence.tsv`: a tabulated file showing the occurrence of major functions in biogeochemical cycles. Rows correspond to the major function and the column corresponds to the community.
- `diagram_*.png`: diagram representing a biogeochemical cycles (carbon, nitrogen, sulfur, other) from METABOLIC. It shows the number of organisms with predicted major functions and the relative occurrence of these functions.
- `function_occurrence.tsv`: a tabulated file containing the ratio for each function. Rows correspond to the function and the column corresponds to the community. It is used to create the `heatmap_occurrence.png` file.
- `function_occurrence_in_organism.tsv`: a tabulated file containing the occurrence of function in each organism of the samples.
- `heatmap_occurrence.png`: a heatmap showing the occurrence for all the HMMs searched by bigecyhmm in the community (all the input protein files).
- `pathway_presence_in_organism.tsv`: a tabulated file containing the occurrence of cycle functions in each organism of the samples.

`bigecyhmm_visualisation.log` is a log file.

`bigecyhmm_visualisation_metadata.json` is a metadata file giving information on the version of the package used.

## 5 Custom usage

If the metabolic functions present in the internal database of bigecyhmm do not cover some metabolic functions of interest, you have several possibilities:
- propose to add some functions (`5.1.1 Contribution to bigecyhmm internal database`)
- modify the internal database (`5.1.2 Modifying bigecyhmm internal database`)
- generate your own custom database (`5.2 bigecyhmm_custom: using custom database`)

### 5.1 Bigecyhmm internal database

#### 5.1.1 Contribution to bigecyhmm internal database

If you are interested in specific functions associated with cycles present in bigecyhmm (carbon, sulfur, nitrogen, phosphorus) and want to propose an addition, you can create an issue or a Pull Request.
Depending on the additions or modifications, it will be taken into account. Keep in mind that bigecyhmm's goal is to limit itself to a small internal database.
If you want to completely add another cycle, please refer to the next subsection.

#### 5.1.2 Modifying bigecyhmm internal database

You can also edit the database to add your own functions. To do so, you can either clone this repository or make a fork. Then install bigecyhmm by running `pip install -e .` inside bigecyhmm folder (where the file `pyproject.toml` is located). You can modify the internal database in different ways. The following files can be modified:

- `hmm_databases/hmm_files`: a folder containing the HMM files used to screen the associated genes. Add the HMM files you want in it.
- `hmm_databases/hmm_table_template.tsv`: a tabulated file containing the association between functions and HMMs. For each HMM you add, you have to add a line in this file. There are two mandatory columns (1) `Hmm file` (name of the HMM file in `hmm_files`) and (2) `Hmm detecting threshold` (threshold used to filter matches).
- `hmm_databases/cycle_pathways.tsv`: a tabulated file linking major functions to HMMs. This file is linked to the creation of the diagrams. If you want to modify this file and propagate the change to the diagram, you must (1) edit diagram templates located at `templates/*` and (2) edit `diagram_cycles.py`, especially function called `create_carbon_cycle` (and the one for the other cycles). In this function several lines are associated with the major function: `data_step_01 = diagram_data['C-S-01:Organic carbon oxidation']` extracts function abundance from predictions, `imgdraw.text((800,80), 'Step1: Organic carbon\n oxidation\n{0}: {1}\n{2}: {3}%'.format(first_term, data_step_01[0], second_term, data_step_01[1]), (0,0,0), font=font)` puts the prediction on the template. Modifying the template requires to also modifies these scripts.

**Adding new HMM**

If you want to add a new HMM for the search, just modify `hmm_databases/hmm_files` and `hmm_databases/hmm_table_template.tsv`.

**Adding new pathway/diagram**

If you want to modify or create a diagram: you have to put the new associated HMMs in `hmm_databases/hmm_files` and `hmm_databases/hmm_table_template.tsv`. Then modify `hmm_databases/cycle_pathways.tsv` by adding the new pathways associated with HMMs. As shown below, this requires to modify script and figure template files, it is easier to rely on the `bigecyhmm_custom` when adding new pathways (see `5.2 bigecyhmm_custom: using custom database`).

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

It is possible to create a custom database that is linked to a specific biogeochemical cycles (or metabolic networks) and then run it with `bigecyhmm_custom`.

#### 5.2.1 Requirements

This command requires three packages:

- [PyHMMER](https://github.com/althonos/pyhmmer): to perform HMM search.
- [networkx](https://github.com/networkx/networkx): to handle biogeochemical cycle as a graph.
- [matplotlib](https://github.com/matplotlib/matplotlib): to create automatically (bad) visualisation of the cycle.

#### 5.2.2 Inputs

This command expects two mandatory arguments:

- `-i`: an input protein sequence fasta file/folder.
- `-d`: a file containing the custom databases (potentially a folder). This file can be either a `.json` file ([carbon cycle json file](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/custom_db/carbon_cycle.json)) or a `.tsv` file ([carbon_cycle_od file](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/custom_db_one_file/carbon_cycle_od.tsv)). If it finds one, it will search for associated folder/file with the same name containing HMMs. Furthermore, if a `.json` file is found, it will search for a `.tsv` corresponding to the `hmm_template`.

Here are several examples of inputs:

- A tsv file and an HMM folder (associated argument `-d custom_db_cycle.tsv`):
```
custom_db_cycle.tsv
custom_db_cycle
```

- use internal custom dabatase `-d internal_hydrogen_table`:

- Only a json file, bigecyhmm will use its internal HMM database to search for HMM files from the json file (associated argument `-d custom_db_cycle.json`):
```
custom_db_cycle.json
```

- A folder with one json file and tsv file/folder (associated argument `-d custom_db_cycle`):
```
custom_db_cycle
├── custom_db_cycle.json
├── custom_db_cycle.tsv
├── custom_db_cycle
```

- A folder with several json files (associated argument `-d custom_db_cycle`):
```
custom_db_cycle
├── carbon_cycle.json
├── carbon_cycle.tsv
├── carbon_cycle
├── nitrogen_cycle.json
├── nitrogen_cycle.tsv
├── nitrogen_cycle
├── sulfur_cycle.json
├── sulfur_cycle.tsv
├── sulfur_cycle
```

##### 5.2.2.1 Inputs: custom tsv file

This format tries to merge both pathway definition and HMM association in one file. It expects a folder (with the same name) containing associated HMMs, such as:

```
├── custom_database.tsv     <- custom tabulated file
├── custom_database         <- folder containing HMM file
```

The input file looks like this:

| Graph_No | ID                      | Type     | Implementation      | HMM_threshold | enzyme_long                               | function_input | function_output | HMM_source                                 | enzyme_ref             | kegg_ortholog       |
|----------|-------------------------|----------|---------------------|---------------|-------------------------------------------|----------------|-----------------|--------------------------------------------|------------------------|---------------------|
| 1.0      | Hydrogen generation     | FUNCTION | hydA1 OR frhA       | NA            | NA                                        | e-             | Hydrogen        | NA                                         | NA                     | NA                  |
| 1.1      | hydA1                   | HMM      | fe.hmm              | 188\|domain   | hydrogenase [FeFe] A1                     |                |                 | doi: 10.3389/fbinf.2022.918853             | DOI: 10.1038/srep34212 | K00532              |
| 1.2      | frhA                    | HMM      | nife-group-3abd.hmm | 246\|domain   | hydrogenase [NiFe] 3a                     |                |                 | doi: 10.3389/fbinf.2022.918853             | DOI: 10.1038/srep34212 | KEGG TOO UNSPECIFIC |
| 4.0      | Methanogen.             | FUNCTION | mcrA                | NA            | NA                                        | Hydrogen, CO2  | Methane         | NA                                         | NA                     | NA                  |
| 4.1      | mcrA                    | HMM      | TIGR03256.hmm       | 314.45\|full  | methyl-coenzyme M reductase alpha subunit |                |                 | https://doi.org/10.1186/s40168-021-01213-8 |                        | K00399              |

The first column (`Graph_No`) corresponds to the number of the metabolic pathway/function in the database (number ending by `.0` corresponds to the function, the ones ending by another number corresponds to a HMM/gene).

The second column (`ID`) indicates the ID/name of the function or the abbreviated name of the gene.

The third column (`Type`) shows if the row is either a metabolic pathway/function (`FUNCTION`) or a gene/HMM (`HMM`).

The fourth column (`Implementation`) is different if the row is a function or a HMM:
- for `FUNCTION`, it corresponds to the boolean expression (such as `hydA1 OR frhA `) showcasing the association of gene/HMM to predict the function. This corresponds to column `HMMs` of `cycle_pathways.tsv` from the internal database of bigecyhmm.
- for `HMM`, it is the name of the HMM file in the HMM folder (such as `fe.hmm`) given with this `.tsv` file.

The fifth column (`HMM_threshold`) is only for row of `Type HMM`. It indicates the threshold to use to filder the detection (such as `188|domain`).

The sixth column (`enzyme_long`) is only for row of `Type HMM` and indicates the name of the gene/enzyme (such as `hydrogenase [FeFe]`).

The seventh and eigth columns (`function_input` and `function_output`) is only for row of `Type FUNCTION` and describes the metabolite given as input to the function. If there are multiple metabolites they are separated bya  `, ` (such as `Hydrogen, CO2`). This is important for the generation of the birpartite graph as this will generate edge between function and metabolite nodes.

Columns after the eight are optional and used to store comments.

There is a full example in bigecyhmm internal database folder ([custom_hydrogen_central_cycles.tsv](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/bigecyhmm/hmm_databases/custom_hydrogen_central_cycles.tsv)).

##### 5.2.2.2 Inputs: custom json file

It is also possible to give a json file as input, such as:

```
├── custom_database.json    <- custom json file
├── custom_database.tsv     <- file similar to hmm_template
├── custom_database         <- folder containing HMM file
```

The three expected files/folders are:
  - a `json` file representing the biogeochemical cycle as a bipartite graph with nodes representing `metabolite` and `function`. Example can be found in the test folder, such as [carbon cycle json file](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/custom_db/carbon_cycle.json). The `hmm` field in the `function node` in the json is mandatory to indicate the HMMs associated with the functions of the cycle. The HMMs are represented as a string with ` OR ` separating HMMs as a `OR` relation (meaning these HMMs are redundant) and ` AND ` as a `AND` relation (meaning that both HMMs are required).
  - a folder containing the HMM profiles (`.hmm` files) such as the one used by bigecyhmm ([hmm_files](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/hmm_files)). If no file is present in the folder, bigecyhmm will use its internal HMM database. You can search for HMM in [KEGG Ortholog database](https://www.genome.jp/kegg/ko.html), [Protein Family Models from NIH](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/), [PFAM](https://www.ebi.ac.uk/interpro/download/Pfam/) or [EggNOG](http://eggnog5.embl.de/#/app/home). It is also possible to build them, an example can be found with [pyhmmer](https://pyhmmer.readthedocs.io/en/stable/examples/msa_to_hmm.html#Build-an-HMM-from-a-multiple-sequence-alignment).

  - a `tsv` file containing the threshold for the different HMMs, similar to hmm_template from bigecyhmm database. If no file is present in the folder, bigecyhmm will use its internal template file for threshold. An example can be found in bigecyhmm internal database ([hmm_table_template.tsv](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/hmm_table_template.tsv)) or in the test folder ([hmm_table_template.tsv](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/mini_custom_db/hmm_table_template.tsv)).

##### 5.2.2.3 Inputs: folder

It is also possible to give an input folder containing multiple custom databases:

```
input_folder
├── custom_database.json
├── custom_database.tsv
├── custom_database
├── other_custom_database.tsv
├── other_custom_database
```

An example with mini database is present in the [test folder](https://github.com/ArnaudBelcour/bigecyhmm/tree/main/test/input_data/mini_custom_db).

##### 5.2.2.3 Inputs: internal custom database

Several custom database (using either json or tsv files) are already present in bigecyhmm and can be called by using a str (such as `-d internal_hydrogen_table`):

- `internal_hydrogen_table`: a tsv custom database focus on hydrogen metabolism ((custom_hydrogen_central_cycles.tsv)[https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/custom_hydrogen_central_cycles.tsv]).
- multiple json files that are similar to the biogeochemical cycles present in bigecyhmm default database:
  - `internal_carbon`: carbon-centered cycle ((custom_carbon_cycle.json)[https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/hmm_databases/custom_carbon_cycle.json]).
  - `internal_sulfur`: sulfur-centered cycle.
  - `internal_nitrogen`: nitrogen-centered cycle.
  - `internal_phosphorus`: phosphorus-centered cycle.
  - `internal_other`: other cycle (arsenate, selenate).
  - `internal_all`: combinations of all the previous cycles.

#### 5.2.3 Usage

Usage example:
```
bigecyhmm_custom -i protein_sequences.faa -d custom_db -o output_folder
```

It can take four optional arguments:

- `-c`: number of cores for multiprocessing.
- `--esmecata`: by giving an EsMeCaTa output folder, `bigecyhmm_custom` maps taxon_id to organism names to associate organism abundance with EsMeCaTa predictions.
- `-m`: JSON file containing gene associated with protein motifs to check for predictions. This verification comes from the [METABOLIC article](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8#Sec2) (you can find information about it, in the section `Motif validation`). The protein motif corresponds to a regex associated with amnio-acids or `X` (the latter being any amino-acid). The idea of this verification is to check if an expected amino-acid motif is present in the sequence matching the associated HMM. You can see an example file in the test folder ([motif.json](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/motif.json)). The name of the gene corresponds to the name of its HMM. If no file is given, it will be using the default ones from METABOLIC (you can find it [here](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/__init__.py#L39) as a dicitonary).
- `-p`: JSON file containing association between two genes to check for predictions. This verification comes from the [METABOLIC article](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01213-8#Sec2) (you can find information about it, in the section `Motif validation`). It ensures that a sequence is properly associated with a specific HMM and not to anotehr yet similar HMM. An example file can be found in the test folfer ([motif_pair.json](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/test/input_data/motif_pair.json)). It contains association between two gene names. The HMM search results of the sequence against these two gnee profiles are compared to find the one with a better score. The name of the gene corresponds to the name of its HMM. If no file is given, it will be using the default ones from METABOLIC (you can find it [here](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/bigecyhmm/__init__.py#L50) as a dicitonary).

#### 5.2.4 Outputs

`bigecyhmm_custom` creates inside the output folder one folder per input custom database file. It generates similar files than bigecyhmm default outputs except for the cycle visualisation. Similar to output folder from bigecyhmm, output folder of `bigecyhmm_custom` can be used as input folder for `bigecyhmm_visualisation`.

```
output_folder
├── input_file_name
├── hmm_results
|   │   ├── org_1.tsv
|   │   ├── org_2.tsv
|   │   ├── org_3.tsv
|   ├── bigecyhmm.log
|   ├── bigecyhmm_metadata.json
|   ├── function_presence.tsv
|   ├── mapping_pathway_to_function_name.tsv
|   ├── pathway_presence.tsv
|   ├── pathway_presence_hmms.tsv
|   ├── R_input.txt
├── database
│   ├── hmm_template_file.tsv
│   ├── input_graph.graphml
│   ├── pathway_template_file.tsv
├── bigecyhmm_custom.log
├── bigecyhmm_custom_metadata.json
```

As there are no diagram templates, `diagram_figures` and `diagram_input` are not generated. `bigecyhmm_custom` relies on network representation to generate graph visualisations.
To do so, it creates network file (`input_graph.graphml`) as output. These files can be used in network software (such as [Cytoscape](https://cytoscape.org/), or [pyvis](https://github.com/WestHealth/pyvis)) to generate visualisation.

Output folder of `bigecyhmm_custom` can be used as input to `bigecyhmm_visualisation`. This will generate similar output files except for the diagram that can not be generated. But a visualisation will be made for the donut plot either: (1) using itnernal template for custom internal database, (2) automatically generated bipartite graph representing the metabolism or (3) you can use your own background image with option `--background-file`.

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

## License

This software is licensed under the GNU GPL-3.0-or-later, see the [LICENSE](https://github.com/ArnaudBelcour/bigecyhmm/blob/main/LICENSE) file for details.

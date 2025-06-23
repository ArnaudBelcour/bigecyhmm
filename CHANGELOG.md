# Changelog

# bigecyhmm v0.1.7 (2025-05-22)

## Add

* Create two files showing occurrence of functions in organism in `bigecyhmm_visualisation` (helpful to link function to organism when using results from EsMeCaTa).
* HMM file `K07699.hmm` to database associated with sporulation regulation gene `spo0A`.

## Modify

* Optimise function abundance computation by using matrix multiplication.
* Update readme according to the additions/changes.

## Fix

* Fix several typos.
* Fix issue with Methanogenesis (two wrongs KOs associated with it).

# bigecyhmm v0.1.6 (2025-04-14)

## Add

* Draft of phosphorous cycles (with genes `PhnD`, `PhnJ`, `PstS`, `PtxD`, `PitA`, `PNaS`, `PhnZ`, `PhnX`, `PhnA`, `PhnM`, `PhnW`, `PepM`, `Ppd`, `HtxB`, `HtxA` and `ptxB`).
* Template figure for phosphorous cycles (from [https://doi.org/10.1038/s41467-024-47914-0](https://doi.org/10.1038/s41467-024-47914-0)).
* Three test HMMs (`HtxB`, `HtxA` and `ptxB` created from [https://doi.org/10.1038/s41467-024-47914-0](https://doi.org/10.1038/s41467-024-47914-0)).
* HMMs for phosphorus cycle from [https://doi.org/10.1016/j.soilbio.2022.108826](https://doi.org/10.1016/j.soilbio.2022.108826): `phoA`, `phoD`, `phoX`, `phoN`, `aphA`, `olpA`, `appA`, `phnG`, `phnH`, `phnI`, `phnK`, `phnL`, `ppa`, `ppx`, `ppk1`, `gcd`, `pqqC`, `phoB`, `phoR`, `phoU`.
* HMMs for carbonic anhydrase: `K01672`, `K01673`, `K01674`, `K01743`, `K01726`, `K18245` and `K18246`.
* HMM for Thiosuflate oxidation second path (S4intermediate): `tsdA`.
* Creation of barplot showing the coverage of EsMeCaTa on each sample of a dataset when using `bigecyhmm_visualisation` with abundance file.
* Test for computing abundance of taxonomic rank (used to create the barplot showing the coverage of EsMeCaTa predicitons on the samples).
* `Aerobic respiration` with Cytochrome c oxidase in other cycle to highlight aerobic taxa. Update associated template.
* Creation of HMM functional profile showing the abundance of each HMM in the sample of the dataset.
* Motif validation step from METABOLIC.
* Prototype to use custom biogeochemical diagram (or metabolic network) and HMM databases with `bigecyhmm_custom`.
* Add examples of custom json database files in internal databases: `custom_carbon_cycle.json`, `custom_phosphorus_cycle.json`, `custom_sulfur_cycle.json`, `custom_other_cycle.json` and `custom_nitrogen_cycle.json`.
* Add `create_html_visualisation_pyvis.py` in scripts folder to test visualisation with `bigecyhmm_custom`.

## Modify

* Move global variables in `__init__.py`.
* Update readme according to the additions/changes.
* Remove several internal files as they were deprecated or not used anymore.

## Fix

* Fix issue of function with zero abundance/occurrence not present in output files.
* Fix several issues with threshold not being applied correctly.

# bigecyhmm v0.1.5 (2025-01-29)

## Add

* Complete refactoring of `bigecyhmm_visualisation` with two new subcommands `bigecyhmm_visualisation esmecata` for input from EsMeCaTa or `bigecyhmm_visualisation genomes` if the inputs are genomes.
* Handle abundance data with `bigecyhmm_visualisation`.
* In `bigecyhmm_visualisation`, add the possibility to create heatmap, polar plot and diagrams showing abundance of specific function categories. Also creates tabulated files indicating the ratio of organisms associated with functions for each sample.
* Add creation of metadata file for `bigecyhmm_visualisation`.
* Numerous tests for the visualisation.
* HMM associated with formate fermentation (formate lyase, hycE, `K15830`).

## Modify

* Rename pathway `C-S-10:Acetogenesis` into `C-S-10:Acetogenesis WL` and make it focus on Wood–Ljungdahl pathway.
* Update readme to better explain outputs and `bigecyhmm_visualisation`.

## Fix

* Fix issue when computing abundance where abundances were not normalised by total abundance.

# bigecyhmm v0.1.4 (2025-01-08)

## Fix

* Try to fix issue in pyproject.toml with pypi publishing.

# bigecyhmm v0.1.3 (2025-01-06)

## Add

* `bigecyhmm_visualisation` to create visualisation from esmecata and bigecyhmm outputs.

# bigecyhmm v0.1.2 (2024-12-09)

## Add

* `pathway_presence.tsv` showing the presence of pathways in each taxon of the communities.
* file associated with hydrogen consumption `hydrogen_consumption.tsv`.
* test for carbon cycle.
* missing HMMs (and other HMMs for hydrogenases).

## Modify

* move `Total.R_input.txt` outside of `diagram_input` folder.

## Fix

* typo in metadata.

# bigecyhmm v0.1.1 (2024-11-01)

## Add

* presence of formula for pathways.

## Fix

* issue in heterogeneity between HMMs presence in pathways, template and compressed file.
* issue in pyproject.toml.

# bigecyhmm v0.1.0 (2024-10-22)

## Add

* use of threshold with pyhmmer.

## Modify

* information for several HMMs: `acsB`, `mcrC`, `mcrG`, `mtaB`, `mtbA`, `rubisco_form_II_III`, `rubisco_form_III`, `rubisco_form_IV`, `fthfs`, `fhs`, `nod`, `apsA`, `apsK`, `arsC (trx)`, `cysN`.

# bigecyhmm v0.0.4 (2024-03-13)

## Add

* scripts to create diagram figures in Python using pillow.
* pillow as dependency.
* templates folder for biogeochemical diagram creation.

# bigecyhmm v0.0.3 (2024-03-12)

## Add

* `scripts` folder with script `create_hmm.py` trying to create HMMs and `draw_biogeochemical_cycles.R` to draw geochemical cycle diagrams (modified from METABOLIC script).
* `R_diagram_pathways.tsv` in database folder.
* creation of metadata file.

## Modify

* update HMMs database with: `aspA`, `nod`, `fthfs`, `mcrG`, `mtaB` and `mtbA`.
* set HMMs score threshold to 40 and evalue to 1e-5 to extract functions.

# bigecyhmm v0.0.2 (2024-03-11)

## Add

* diagram files creation.
* `R_diagram_pathways.tsv` in database folder.

## Modify

* multiprocessing: previously done using pyhmmer multiprocessing but it is much faster to use python multiprocessing on the different input fasta files and give 1 CPU to pyhmmer.

# bigecyhmm v0.0.1 (2024-03-08)

First test release for bigecyhmm.

## Add

* bigecyhmm command line and python functions.
* HMMs used in `hmm_databases` folder.
* Readme.
* CHANGELOG.
* License.

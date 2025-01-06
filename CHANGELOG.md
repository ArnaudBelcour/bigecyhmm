# Changelog

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

# Copyright (C) 2024-2025 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

__version__ = '0.1.6'

# Create global constants containing path to package internal database.
import os
ROOT = os.path.dirname(__file__)
PATHWAY_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'cycle_pathways.tsv')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
HMM_COMPRESS_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
PHENOTYPE_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'phenotypes.tsv')
HYDROGEN_CONSUMPTION_FILE = os.path.join(ROOT, 'hmm_databases', 'hydrogen_consumption.tsv')
TEMPLATE_CARBON_CYCLE = os.path.join(ROOT, 'templates', 'template_carbon_cycle_total.png')
TEMPLATE_NITROGEN_CYCLE = os.path.join(ROOT, 'templates', 'template_nitrogen_cycle_total.png')
TEMPLATE_SULFUR_CYCLE = os.path.join(ROOT, 'templates', 'template_sulfur_cycle_total.png')
TEMPLATE_OTHER_CYCLE = os.path.join(ROOT, 'templates', 'template_other_cycle_total.png')
TEMPLATE_PHOSPHORUS_CYCLE = os.path.join(ROOT, 'templates', 'template_phosphorus_cycle.png')
TEMPLATE_PHOSPHORUS_GENE_CYCLE = os.path.join(ROOT, 'templates', 'template_phosphorus_genes.png')
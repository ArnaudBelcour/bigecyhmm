# Copyright (C) 2024 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
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

import csv
import os

ROOT = os.path.dirname(__file__)
HMM_COMPRESS_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
DIAGRAM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'R_diagram_pathways.tsv')

def get_diagram_pathways_hmms():
    pathway_hmms = {}
    sorted_pathways = []
    with open(DIAGRAM_TEMPLATE_FILE, 'r') as open_r_pathways:
        csvreader = csv.reader(open_r_pathways, delimiter = '\t')
        for line in csvreader:
            sorted_pathways.append(line[0])
            if ';' in line[1]:
                hmm_combinations = [combination.split(',') for combination in line[1].split(';')] 
                pathway_hmms[line[0]] = hmm_combinations
            else:
                pathway_hmms[line[0]] = [line[1].split(',')]

    sorted_pathways = sorted(sorted_pathways)

    return pathway_hmms, sorted_pathways


def get_organism_matching_hmms(hmm_folder):
    org_hmms = {}
    for org_hmm_file in os.listdir(hmm_folder):
        org_hmm_file_basename = os.path.splitext(org_hmm_file)[0]
        org_hmm_file_path = os.path.join(hmm_folder, org_hmm_file)
        with open(org_hmm_file_path, 'r') as open_hmm_file:
            csvreader = csv.DictReader(open_hmm_file, delimiter = '\t')
            hmms = []
            next(csvreader)
            for line in csvreader:
                if float(line['evalue']) < 1e-5 and float(line['score']) >= 40:
                    hmms.append(line['HMM'].replace('_full', ''))
            org_hmms[org_hmm_file_basename] = hmms

    return org_hmms

def check_diagram_pathways(sorted_pathways, org_hmms, pathway_hmms):
    all_pathways = {pathway: 0 for pathway in sorted_pathways}
    org_pathways = {}
    for org in org_hmms:
        for pathway in sorted_pathways:
            pathway_check = []
            for hmm_combination in pathway_hmms[pathway]:
                negative_hmms = [hmm for hmm in hmm_combination if 'NO' in hmm]
                if len(negative_hmms) > 0:
                    hmm_combination = [hmm.replace('NO|', '') for hmm in hmm_combination]
                intersection_hmms = set(hmm_combination).intersection(org_hmms[org])
                if len(intersection_hmms) > 0:
                    if len(negative_hmms) == 0:
                        pathway_check.append(True)
                    else:
                        pathway_check.append(False)
                else:
                    pathway_check.append(False)


            if all(pathway_check) is True:
                if org not in org_pathways:
                    org_pathways[org] = {}
                org_pathways[org][pathway] = 1
                if pathway not in all_pathways:
                    all_pathways[pathway] = 1
                else:
                    all_pathways[pathway] += 1
            else:
                if org not in org_pathways:
                    org_pathways[org] = {}
                org_pathways[org][pathway] = 0

    return all_pathways, org_pathways

def create_input_diagram(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    pathway_hmms, sorted_pathways = get_diagram_pathways_hmms()
    org_hmms = get_organism_matching_hmms(input_folder)
    all_pathways, org_pathways = check_diagram_pathways(sorted_pathways, org_hmms, pathway_hmms)

    for org in org_pathways:
        org_file = os.path.join(output_folder, org+'.R_input.txt')
        with open(org_file, 'w') as open_output_file:
            csvwriter = csv.writer(open_output_file, delimiter='\t')
            for pathway in org_pathways[org]:
                csvwriter.writerow([pathway, org_pathways[org][pathway]])

    total_file = os.path.join(output_folder, 'Total.R_input.txt')
    with open(total_file, 'w') as open_total_file:
        csvwriter = csv.writer(open_total_file, delimiter='\t')
        for pathway in all_pathways:
            csvwriter.writerow([pathway, all_pathways[pathway], all_pathways[pathway] / len(org_hmms)])


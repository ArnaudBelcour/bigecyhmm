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

import logging
from PIL import Image, ImageDraw, ImageFont

from bigecyhmm.utils import parse_result_files

ROOT = os.path.dirname(__file__)
HMM_COMPRESS_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
DIAGRAM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'R_diagram_pathways.tsv')
TEMPLATE_CARBON_CYCLE = os.path.join(ROOT, 'templates', 'template_carbon_cycle_total.png')
TEMPLATE_NITROGEN_CYCLE = os.path.join(ROOT, 'templates', 'template_nitrogen_cycle_total.png')
TEMPLATE_SULFUR_CYCLE = os.path.join(ROOT, 'templates', 'template_sulfur_cycle_total.png')
TEMPLATE_OTHER_CYCLE = os.path.join(ROOT, 'templates', 'template_other_cycle_total.png')

logger = logging.getLogger(__name__)


def get_diagram_pathways_hmms():
    """From DIAGRAM_TEMPLATE_FILE extract HMMs associated with cycles of the diagrams.

    Returns:
        pathway_hmms (dict): dictionary with functions as key and list of list of HMMS as value
        sorted_pathways (list): ordered list of functions
    """
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


def check_diagram_pathways(sorted_pathways, org_hmms, pathway_hmms):
    """Compute the presence of functions of biogeochemical cycles in the dataset.

    Args:
        sorted_pathways (list): ordered list of functions
        org_hmms (dict): dictionary with organism as key and list of hit HMMs as value
        pathway_hmms (dict): dictionary with functions as key and list of list of HMMS as value

    Returns:
        all_pathways (dict): pathway as key and number of organisms having it as value
        org_pathways (dict): organism as key and subdict with pathway presence as value
    """
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
    """Create input files for the creation of the biogeochemical cycle diagram.
    This function creates input for this R script: https://github.com/AnantharamanLab/METABOLIC/blob/master/draw_biogeochemical_cycles.R

    Args:
        input_folder (str): path to HMM search results folder (one tsv file per organism)
        output_folder (str): path to output folder containing input files for diagram creation
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    pathway_hmms, sorted_pathways = get_diagram_pathways_hmms()
    org_hmms = parse_result_files(input_folder)
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


def parse_diagram_folder(input_diagram_folder):
    """Parse functions in Total.R_input.txt.

    Args:
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function

    Returns:
        diagram_data (dict): functions as key and (nb genomes containing in it, percentage coverage) as value
    """
    diagram_data_file = os.path.join(input_diagram_folder, 'Total.R_input.txt')

    diagram_data = {}
    with open(diagram_data_file, 'r') as open_diagram_data_file:
        csvreader = csv.reader(open_diagram_data_file, delimiter='\t')
        for line in csvreader:
            nb_genomes =  line[1]
            percentage_coverage = round(float(line[2]) * 100, 1)
            diagram_data[line[0]] = [nb_genomes, percentage_coverage]

    return diagram_data


def create_carbon_cycle(output_file, input_diagram_folder):
    """From png TEMPLATE_CARBON_CYCLE and input_diagram_folder file, create carbon cycle figure.

    Args:
        output_file (str): path to output file
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function
    """
    img = Image.open(TEMPLATE_CARBON_CYCLE, 'r')
    imgdraw = ImageDraw.Draw(img)
    font = ImageFont.load_default(20)
    diagram_data = parse_diagram_folder(input_diagram_folder)

    data_step_01 = diagram_data['C-S-01:Organic carbon oxidation']
    data_step_02 = diagram_data['C-S-02:Carbon fixation']
    data_step_03 = diagram_data['C-S-03:Ethanol oxidation']
    data_step_04 = diagram_data['C-S-04:Acetate oxidation']
    data_step_05 = diagram_data['C-S-05:Hydrogen generation']
    data_step_06 = diagram_data['C-S-06:Fermentation']
    data_step_07 = diagram_data['C-S-07:Methanogenesis']
    data_step_08 = diagram_data['C-S-08:Methanotrophy']
    data_step_09 = diagram_data['C-S-09:Hydrogen oxidation']
    data_step_10 = diagram_data['C-S-10:Acetogenesis']

    imgdraw.text((800,80), 'Step1: Organic carbon\n oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_01[0], data_step_01[1]), (0,0,0), font=font)
    imgdraw.text((100,70), 'Step2: Carbon fixation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_02[0], data_step_02[1]), (139,137,137), font=font)
    imgdraw.text((750,320), 'Step3: Ethanol oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_03[0], data_step_03[1]), (0,0,0), font=font)
    imgdraw.text((150,400), 'Step4: Acetate oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_04[0], data_step_04[1]), (0,0,0), font=font)
    imgdraw.text((530,225), 'Step5: Hydrogen generation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_05[0], data_step_05[1]), (139,117,0), font=font)
    imgdraw.text((375,150), 'Step6: Fermentation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_06[0], data_step_06[1]), (139,117,0), font=font)
    imgdraw.text((350,450), 'Step7: Methanogenesis\nGenomes: {0}\nCoverage: {1}%'.format(data_step_07[0], data_step_07[1]), (93,71,139), font=font)
    imgdraw.text((300,650), 'Step8: Methanotrophy\nGenomes: {0}\nCoverage: {1}%'.format(data_step_08[0], data_step_08[1]), (205,186,150), font=font)
    imgdraw.text((575,400), 'Step9: Hydrogen oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_09[0], data_step_09[1]), (238,162,173), font=font)
    imgdraw.text((275,300), 'Step10: Acetogenesis\nGenomes: {0}\nCoverage: {1}%'.format(data_step_10[0], data_step_10[1]), (0,134,139), font=font)

    img = img.resize((2112, 1632), Image.Resampling.LANCZOS)
    img.save(output_file, dpi=(300, 300), quality=100)


def create_nitrogen_cycle(output_file, input_diagram_folder):
    """From png TEMPLATE_NITROGEN_CYCLE and input_diagram_folder file, create nitrogen cycle figure.

    Args:
        output_file (str): path to output file
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function
    """
    img = Image.open(TEMPLATE_NITROGEN_CYCLE, 'r')
    imgdraw = ImageDraw.Draw(img)
    font = ImageFont.load_default(20)
    diagram_data = parse_diagram_folder(input_diagram_folder)

    data_step_01 = diagram_data['N-S-01:Nitrogen fixation']
    data_step_02 = diagram_data['N-S-02:Ammonia oxidation']
    data_step_03 = diagram_data['N-S-03:Nitrite oxidation']
    data_step_04 = diagram_data['N-S-04:Nitrate reduction']
    data_step_05 = diagram_data['N-S-05:Nitrite reduction']
    data_step_06 = diagram_data['N-S-06:Nitric oxide reduction']
    data_step_07 = diagram_data['N-S-07:Nitrous oxide reduction']
    data_step_08 = diagram_data['N-S-08:Nitrite ammonification']
    data_step_09 = diagram_data['N-S-09:Anammox']
    data_step_10 = diagram_data['N-S-10:Nitric oxide dismutase']

    imgdraw.text((700,120), 'Step1: Nitrogen fixation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_01[0], data_step_01[1]), (205,16,118), font=font)
    imgdraw.text((800,360), 'Step2: Ammonia oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_02[0], data_step_02[1]), (0,205,205), font=font)
    imgdraw.text((650,650), 'Step3: Nitrite oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_03[0], data_step_03[1]), (139,69,0), font=font)
    imgdraw.text((250,600), 'Step4: Nitrate reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_04[0], data_step_04[1]), (16,78,139), font=font)
    imgdraw.text((50,425), 'Step5: Nitrite reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_05[0], data_step_05[1]), (16,78,139), font=font)
    imgdraw.text((50,300), 'Step6: Nitric oxide reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_06[0], data_step_06[1]), (16,78,139), font=font)
    imgdraw.text((225,120), 'Step7: Nitrous oxide reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_07[0], data_step_07[1]), (16,78,139), font=font)
    imgdraw.text((500,425), 'Step8: Nitrite ammonification\nGenomes: {0}\nCoverage: {1}%'.format(data_step_08[0], data_step_08[1]), (95,158,160), font=font)
    imgdraw.text((500,275), 'Step9: Anammox\nGenomes: {0}\nCoverage: {1}%'.format(data_step_09[0], data_step_09[1]), (102,205,0), font=font)
    imgdraw.text((400,200), 'Step10: Nitric oxide dismutase\nGenomes: {0}\nCoverage: {1}%'.format(data_step_10[0], data_step_10[1]), (154,50,205), font=font)

    img = img.resize((2112, 1632), Image.Resampling.LANCZOS)
    img.save(output_file, dpi=(300, 300), quality=100)


def create_sulfur_cycle(output_file, input_diagram_folder):
    """From png TEMPLATE_SULFUR_CYCLE and input_diagram_folder file, create sulfur cycle figure.

    Args:
        output_file (str): path to output file
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function
    """
    img = Image.open(TEMPLATE_SULFUR_CYCLE, 'r')
    imgdraw = ImageDraw.Draw(img)
    font = ImageFont.load_default(20)
    diagram_data = parse_diagram_folder(input_diagram_folder)

    data_step_01 = diagram_data['S-S-01:Sulfide oxidation']
    data_step_02 = diagram_data['S-S-02:Sulfur reduction']
    data_step_03 = diagram_data['S-S-03:Sulfur oxidation']
    data_step_04 = diagram_data['S-S-04:Sulfite oxidation']
    data_step_05 = diagram_data['S-S-05:Sulfate reduction']
    data_step_06 = diagram_data['S-S-06:Sulfite reduction']
    data_step_07 = diagram_data['S-S-07:Thiosulfate oxidation']
    data_step_08 = diagram_data['S-S-08:Thiosulfate disproportionation 1']
    data_step_09 = diagram_data['S-S-09:Thiosulfate disproportionation 2']

    imgdraw.text((700,80), 'Step1: Sulfide oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_01[0], data_step_01[1]), (238,118,0), font=font)
    imgdraw.text((600,200), 'Step2: Sulfur reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_02[0], data_step_02[1]), (122,197,205), font=font)
    imgdraw.text((850,360), 'Step3: Sulfur oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_03[0], data_step_03[1]), (154,50,205), font=font)
    imgdraw.text((650,650), 'Step4: Sulfite oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_04[0], data_step_04[1]), (162,205,90), font=font)
    imgdraw.text((100,550), 'Step5: Sulfate reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_05[0], data_step_05[1]), (139,69,19), font=font)
    imgdraw.text((150,150), 'Step6: Sulfite reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_06[0], data_step_06[1]), (139,69,19), font=font)
    imgdraw.text((375,500), 'Step7: Thiosulfate oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_07[0], data_step_07[1]), (0,104,139), font=font)
    imgdraw.text((400,250), 'Step8: Thiosulfate \ndisproportionation 1\nGenomes: {0}\nCoverage: {1}%'.format(data_step_08[0], data_step_08[1]), (0,104,139), font=font)
    imgdraw.text((625,400), 'Step9: Thiosulfate \ndisproportionation 2\nGenomes: {0}\nCoverage: {1}%'.format(data_step_09[0], data_step_09[1]), (0,104,139), font=font)

    img = img.resize((2112, 1632), Image.Resampling.LANCZOS)
    img.save(output_file, dpi=(300, 300), quality=100)


def create_other_cycle(output_file, input_diagram_folder):
    """From png TEMPLATE_OTHER_CYCLE and input_diagram_folder file, create other cycle figure.

    Args:
        output_file (str): path to output file
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function
    """
    img = Image.open(TEMPLATE_OTHER_CYCLE, 'r')
    imgdraw = ImageDraw.Draw(img)
    font = ImageFont.load_default(20)
    diagram_data = parse_diagram_folder(input_diagram_folder)

    data_step_01 = diagram_data['O-S-01:Iron reduction']
    data_step_02 = diagram_data['O-S-02:Iron oxidation']
    data_step_03 = diagram_data['O-S-03:Arsenate reduction']
    data_step_04 = diagram_data['O-S-04:Arsenite oxidation']
    data_step_05 = diagram_data['O-S-05:Selenate reduction']

    imgdraw.text((100,175), 'Step1: Iron reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_01[0], data_step_01[1]), (0,100,0), font=font)
    imgdraw.text((375,175), 'Step2: Iron oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_02[0], data_step_02[1]), (0,100,0), font=font)
    imgdraw.text((10,575), 'Step3: Arsenate reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_03[0], data_step_03[1]), (205,102,0), font=font)
    imgdraw.text((330,575), 'Step4: Arsenite oxidation\nGenomes: {0}\nCoverage: {1}%'.format(data_step_04[0], data_step_04[1]), (205,102,0), font=font)
    imgdraw.text((800,575), 'Step5: Selenate reduction\nGenomes: {0}\nCoverage: {1}%'.format(data_step_05[0], data_step_05[1]), (0,0,0), font=font)

    img = img.resize((2112, 1632), Image.Resampling.LANCZOS)
    img.save(output_file, dpi=(300, 300), quality=100)


def create_diagram_figures(hmm_diagram_folder, biogeochemical_diagram_folder):
    """From png TEMPLATE_OTHER_CYCLE and input_diagram_folder file, create other cycle figure.

    Args:
        input_diagram_folder (str): path to input diagram files created by create_input_diagram function
        biogeochemical_diagram_folder (str): path to output folder containing diagram cycles
    """
    logger.info('Creating biogeochemical cycle figures.')

    if not os.path.exists(biogeochemical_diagram_folder):
        os.mkdir(biogeochemical_diagram_folder)

    carbon_cycle_file = os.path.join(biogeochemical_diagram_folder, 'carbon_cycle.png')
    create_carbon_cycle(carbon_cycle_file, hmm_diagram_folder)
    nitrogen_cycle_file = os.path.join(biogeochemical_diagram_folder, 'nitrogen_cycle.png')
    create_nitrogen_cycle(nitrogen_cycle_file, hmm_diagram_folder)
    sulfur_cycle_file = os.path.join(biogeochemical_diagram_folder, 'sulfur_cycle.png')
    create_sulfur_cycle(sulfur_cycle_file, hmm_diagram_folder)
    other_cycle_file = os.path.join(biogeochemical_diagram_folder, 'other_cycle.png')
    create_other_cycle(other_cycle_file, hmm_diagram_folder)
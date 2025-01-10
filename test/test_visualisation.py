import os
import csv
import subprocess
import shutil

from bigecyhmm.visualisation import compute_relative_abundance_per_tax_id, read_abundance_file, read_esmecata_proteome_file

def test_compute_relative_abundance_per_tax_id():
    sample_abundance = {'sample_1': {'org_1': 100, 'org_2': 100, 'org_3': 0},
                        'sample_2': {'org_1': 0, 'org_2': 200, 'org_3': 600},
                        'sample_3': {'org_1': 0, 'org_2': 120, 'org_3': 400}}

    sample_tot_abundance = {'sample_1': 200,
                  'sample_2': 800,
                  'sample_3': 520}

    observation_names_tax_id_names = {'org_1': 'tax_id_name_1', 'org_2': 'tax_id_name_2', 'org_3': 'tax_id_name_2'}

    abundance_data = compute_relative_abundance_per_tax_id(sample_abundance, sample_tot_abundance, observation_names_tax_id_names)
    expected_abundance_data = {'sample_1': {'tax_id_name_1': 0.5, 'tax_id_name_2': 0.5},
                               'sample_2': {'tax_id_name_1': 0, 'tax_id_name_2': 1},
                               'sample_3': {'tax_id_name_1': 0, 'tax_id_name_2': 1}}

    for sample in expected_abundance_data:
        for tax_id_name in expected_abundance_data[sample]:
            assert expected_abundance_data[sample][tax_id_name] == abundance_data[sample][tax_id_name]


def test_compute_relative_abundance_per_tax_id_file():
    proteome_tax_id_file = os.path.join('input_data', 'proteome_tax_id.tsv')
    abundance_file_path = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')

    sample_abundance, sample_tot_abundance = read_abundance_file(abundance_file_path)
    observation_names_tax_id_names = read_esmecata_proteome_file(proteome_tax_id_file)
    abundance_data = compute_relative_abundance_per_tax_id(sample_abundance, sample_tot_abundance, observation_names_tax_id_names)

    expected_abundance_data = {'sample_1': {'tax_id_name_1': 0.5, 'tax_id_name_2': 0.5},
                               'sample_2': {'tax_id_name_1': 0, 'tax_id_name_2': 1},
                               'sample_3': {'tax_id_name_1': 0, 'tax_id_name_2': 1}}

    for sample in expected_abundance_data:
        for tax_id_name in expected_abundance_data[sample]:
            assert expected_abundance_data[sample][tax_id_name] == abundance_data[sample][tax_id_name]

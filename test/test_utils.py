import os
import csv
import zipfile

from bigecyhmm.utils import read_measures_file

def test_read_measures_file():
    abundance_file_path = os.path.join('input_data', 'proteome_tax_id_abundance.tsv')
    sample_abundance, sample_tot_abundance = read_measures_file(abundance_file_path)

    expected_sample_abundance = {'sample_1': {'org_1': 100, 'org_2': 100, 'org_3': 0},
                        'sample_2': {'org_1': 0, 'org_2': 200, 'org_3': 600},
                        'sample_3': {'org_1': 0, 'org_2': 120, 'org_3': 400}}

    expected_sample_tot_abundance = {'sample_1': 200,
                  'sample_2': 800,
                  'sample_3': 520}

    for sample in expected_sample_abundance:
        for tax_id_name in expected_sample_abundance[sample]:
            assert expected_sample_abundance[sample][tax_id_name] == sample_abundance[sample][tax_id_name]

    for sample in expected_sample_tot_abundance:
        assert expected_sample_tot_abundance[sample] == sample_tot_abundance[sample]
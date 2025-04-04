import os

from bigecyhmm.utils import read_measures_file, read_esmecata_proteome_file


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


def test_read_esmecata_proteome_file():
    proteome_tax_id_file = os.path.join('input_data', 'esmecata_output_folder', '0_proteomes', 'proteome_tax_id.tsv')
    observation_names_tax_id_names, observation_names_tax_ranks = read_esmecata_proteome_file(proteome_tax_id_file)

    expected_observation_names_tax_id_names = {'org_1': 'tax_id_name_1', 'org_2': 'tax_id_name_2', 'org_3': 'tax_id_name_2'}
    expected_observation_names_tax_ranks = {'org_1': 'genus', 'org_2': 'family', 'org_3': 'order'}

    for org in expected_observation_names_tax_id_names:
        assert expected_observation_names_tax_id_names[org] == observation_names_tax_id_names[org]

    for org in expected_observation_names_tax_ranks:
        assert expected_observation_names_tax_ranks[org] == observation_names_tax_ranks[org]

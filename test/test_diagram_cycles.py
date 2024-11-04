import os
import csv
import zipfile

from bigecyhmm.diagram_cycles import check_diagram_pathways

def test_check_diagram_pathways():
    sorted_pathways = ['S-S-09:Thiosulfate disproportionation 2']
    org_hmms = {'org_1': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm'],
                'org_2': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm', 'soxC.hmm', 'soxD.hmm'],
                'org_3': ['soxC.hmm', 'soxD.hmm']}
    pathway_hmms = {'S-S-09:Thiosulfate disproportionation 2': [['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm'], ['NO|soxC.hmm', 'NO|soxD.hmm']]}

    expected_org_pathways = {'org_1': {'S-S-09:Thiosulfate disproportionation 2': 1},
                             'org_2': {'S-S-09:Thiosulfate disproportionation 2': 0},
                             'org_3': {'S-S-09:Thiosulfate disproportionation 2': 0}}
    all_pathways, org_pathways, org_pathways_hmms = check_diagram_pathways(sorted_pathways, org_hmms, pathway_hmms)

    for org in org_pathways:
        for pathway in org_pathways[org]:
            assert org_pathways[org][pathway] == expected_org_pathways[org][pathway]
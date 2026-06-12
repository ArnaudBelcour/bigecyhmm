import os
import csv
import zipfile

from bigecyhmm.diagram_cycles import check_diagram_pathways, check_boolean_expression

def test_check_diagram_pathways():
    sorted_pathways = ['S-S-09:Thiosulfate disproportionation 2']
    org_hmms = {'org_1': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm'],
                'org_2': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm', 'soxC.hmm', 'soxD.hmm'],
                'org_3': ['soxC.hmm', 'soxD.hmm']}
    pathway_hmms = {'S-S-09:Thiosulfate disproportionation 2': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm', 'soxC.hmm', 'soxD.hmm']}

    expected_org_pathways = {'org_1': {'S-S-09:Thiosulfate disproportionation 2': 1},
                             'org_2': {'S-S-09:Thiosulfate disproportionation 2': 0},
                             'org_3': {'S-S-09:Thiosulfate disproportionation 2': 0}}
    pathway_expression = {'S-S-09:Thiosulfate disproportionation 2': "(soxX.hmm or soxY.hmm or soxZ.hmm or soxA.hmm) and (not soxC.hmm and not soxD.hmm)"}
    all_pathways, org_pathways, org_pathways_hmms = check_diagram_pathways(sorted_pathways, pathway_expression, org_hmms, pathway_hmms)

    for org in org_pathways:
        for pathway in org_pathways[org]:
            assert org_pathways[org][pathway] == expected_org_pathways[org][pathway]

def test_check_boolean_expression():
    sorted_pathways = ['S-S-09:Thiosulfate disproportionation 2']
    org_hmms = {'org_1': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm'],
                'org_2': ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm', 'soxC.hmm', 'soxD.hmm'],
                'org_3': ['soxC.hmm', 'soxD.hmm']}
    pathway_hmms = ['soxX.hmm', 'soxY.hmm', 'soxZ.hmm', 'soxA.hmm', 'soxC.hmm', 'soxD.hmm']

    expected_org_pathways = {'org_1': True,
                             'org_2': False,
                             'org_3': False}
    hmm_boolean_expression = 	"(soxX.hmm or soxY.hmm or soxZ.hmm or soxA.hmm) and (not soxC.hmm and not soxD.hmm)"
    pathway_presences = {}
    for org in org_hmms:
        pathway_presence = org_pathways = check_boolean_expression(hmm_boolean_expression, org_hmms[org], pathway_hmms)
        pathway_presences[org] = pathway_presence

    for org in expected_org_pathways:
        assert pathway_presences[org] == expected_org_pathways[org]

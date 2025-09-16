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

__version__ = '0.1.8'

# Create global constants containing path to package internal database.
import os
ROOT = os.path.dirname(__file__)
PATHWAY_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'cycle_pathways.tsv')
HMM_TEMPLATE_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_table_template.tsv')
HMM_COMPRESSED_FILE = os.path.join(ROOT, 'hmm_databases', 'hmm_files.zip')
TEMPLATE_CARBON_CYCLE = os.path.join(ROOT, 'templates', 'template_carbon_cycle_total.png')
TEMPLATE_NITROGEN_CYCLE = os.path.join(ROOT, 'templates', 'template_nitrogen_cycle_total.png')
TEMPLATE_SULFUR_CYCLE = os.path.join(ROOT, 'templates', 'template_sulfur_cycle_total.png')
TEMPLATE_OTHER_CYCLE = os.path.join(ROOT, 'templates', 'template_other_cycle_total.png')
TEMPLATE_PHOSPHORUS_CYCLE = os.path.join(ROOT, 'templates', 'template_phosphorus_cycle.png')
TEMPLATE_PHOSPHORUS_GENE_CYCLE = os.path.join(ROOT, 'templates', 'template_phosphorus_genes.png')

CUSTOM_CARBON_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'custom_carbon_cycle.json')
CUSTOM_SULFUR_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'custom_sulfur_cycle.json')
CUSTOM_NITROGEN_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'custom_nitrogen_cycle.json')
CUSTOM_PHOSPHORUS_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'custom_phosphorus_cycle.json')
CUSTOM_OTHER_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'custom_other_cycle.json')
CUSTOM_HYDROGENOTROPHIC_CYCLE_NETWORK = os.path.join(ROOT, 'hmm_databases', 'hydrogenotrophic_salt_caverns.json')

# Motif validation using regex search in sequence.
MOTIF = {'dsrC': r'GPXKXXCXXXGXPXPXXCX',
        'tusE': r'G[PV]XKXX[ARNDQEGHILKMFPSTWYV]XXXGXPXPXXCX',
        'rubisco_form_I': r'TXKPKLGLXXXNYXRXXXEXLXGGL',
        'rubisco_form_III': r'QXGGGXXXHPXGXXAGAXAXRXXXXA',
        'dsrA': r'XXXDXGXXGXXXRXXXXCXGXXXCX',
        'dsrB': r'XXXXPXXXXXXXXXCXXXCXXXX',
        'soxB': r'[YEVRQAHMNTIFG][PRKQNHMFLTGSY][EDGKNS][VACTSLG][EDQNP][VIAFLM][VASGC][FMLT][STY][PA][AG][VFY][R][W][G][TPAGNYLVS][TSCAV][ILVKAMT][LPMIV][PAGEMS][GNED][QDEYASH][AEPTMNDVSGKQW][IVL][TLRM][WRMQLFVKYI][DEGA][HDRNASLWQK][LIVM][YWLMTHAIFS][ANDETGSH][YVQEMANFWH][TCVLM][GSAC][FMICVTL][TSGN][YD][PGST][EAWYQNMDHKS][LATVCSMNI][YTGCFS][LRVKTAISPM][FTRSNQMKGLAD][YENDSPQTRGAWKM][LMVIYRF][RTSEKDNL][GA][AEMNGQTSKRD][QFMTDKRHEYVAG][ILMV][KHARLQN][AVTNDKIEMGLSQ][VILMHT][LMIWFV][E][DGQSKE][IVAL][AGC][SDEQ][NAK][VLI][F][TNQHVSD][SPAKEQTLRVDN][DN][PA][FYL][YLFRIMQK][QIRH][QSG][GE][GQ][D][VM][SVIL][R][VLTIMA][FGEHYTAW][GNA][LIMAVF][RGQEHDTSANIK][YWF][VRTDASKEHMN][LICFVM][DETHAKNRSQM][PVIL][DTNGSAVKMLRYEH][AKNEQLRG][PRGKTESADNI][TQMILNSVAKGRF][GNYSH][ESQKNHRADG][R]',
        'soxC': r'[DERGSHQ][VAICTS][LIFMV][VLI][AVCG][YLFWI][AYKGFSRTN][QAM][N][G][E][ARMPSH][LI][RYM][PKRAVDM][EQSAG][QNH][G][YF][P][VLIAM][R][LAVIM][VILMFC][VLIANM][P][G][WVYCFLI][EQ][GA][SVNI][ITMSLVA][QNWHCS][VIT][K][WYHF][LIV][RHKQN][RQN][ILVM][LQEYKGNMHD][VLFAI][TAMGSHVYDIKL][DEARTNSQ][LQGMKVEARTH][PAERKS][AVYWFMLESTGI][MEANQYGIHVWF][ASTHCLQV][KRYF][DNESGQWF][E][TAVS][SAVIKTLGR][EKHRGYNL][Y][TIVNSM][DEVSNMAGTIQL][VLTMPSAIHQ][TLMIQVKARY][APEGKQDRS][DSNGTE][GKS][RQTLKCIFMHSDVEN][VASHLYWIFTQ][LRQIKYWMTEHV][AQKIRMGENL][FWYHM][TSFAVDNHQ][WFSLYMNG][VIDYPEALSTFMRH][MQLNIFC][EDNRGH][PACVST][QKNRDE][S]}',
}

# Motif validation by checking that it is not better associated with another HMM.
MOTIF_PAIR = {'dsrE': 'tusD',
    'tusD': 'dsrE',
    'dsrH': 'tusB',
    'tusB': 'dsrH',
    'dsrF': 'tusC',
    'tusC': 'dsrF',
    'amoA': 'pmoA',
    'amoB': 'pmoB',
    'amoC': 'pmoC',
    'pmoA': 'amoA',
    'pmoB': 'amoB',
    'pmoC': 'amoC'
}
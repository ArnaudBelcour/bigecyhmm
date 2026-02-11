#This script was made to allow a user of the custom version of bigecyhmm to provide a single document (easier for user) with Functions and HMMs
#and split it into the three files needed for bigecyhmm_custom to run: 
# cycle_pathways_custom.tsv , Hmm_table_template_custom.tsv and a nx.Bigraph object with nodes for every function and substance, which can be saved as a graphml file. 

import csv
import sys
from pathlib import Path
import re
import networkx as nx


def read_input(path: Path):    #expects a TSV file with header 
    with path.open(newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        rows = list(reader)
    return rows[1:] if len(rows) > 1 else []


def parse_rows(rows):
    cycle_rows = []  
    hmm_rows = []    
    func_to_path = {}  
    id_to_hmms = {}    

    last_path = None
    for row in rows:  #pad to 9 columns, and assign variables for input columns (input = i, 0-based column numbers. e.g. i0 = input col 0)
        row = [c.strip() for c in row] + [""] * max(0, 9 - len(row))
        i0, i1, i2, i3, i4, i5, i8 = row[0], row[1], row[2], row[3], row[4], row[5], row[8]
        typ = (i2 or "").upper()

        if 'FUNCTION' in typ:
            func_major = i0.split('.', 1)[0] if i0 else '' #remove dot from ID to get major function ID
            add_piece = i3  #expression stored in FUNCTION rows
            func_to_path[func_major] = i1
            if i1 != last_path:
                cycle_rows.append([i1, add_piece, func_major])
                last_path = i1
            else:
                if add_piece:
                    #separate concatenated pieces by a single space 
                    cycle_rows[-1][1] = (cycle_rows[-1][1] + ' ' + add_piece).strip()

        elif 'HMM' in typ:
            func_major = i0.split('.', 1)[0] if i0 else '' #remove dot from ID to get major function ID
            row11 = [''] * 11
            row11[0] = func_major
            row11[1] = func_to_path.get(func_major, '')
            row11[3] = i1
            row11[4] = i5
            row11[5] = i3
            row11[6] = i8
            row11[10] = i4
            hmm_rows.append(row11)

            if i1:
                id_to_hmms.setdefault(func_major, {})
                if i1 not in id_to_hmms[func_major]:
                    id_to_hmms[func_major][i1] = i3

    return cycle_rows, hmm_rows, id_to_hmms, func_to_path


def translate_expression(expr: str, func_major: str, id_to_hmms: dict) -> str:
    #translate FUNCTION expressions: normalize textual operators to lowercase
    #and replace ID tokens with HMM filenames within the same func_major.
    if not expr:
        return expr

    #transform text operators to lowercase
    expr = re.sub(r"\bAND\b", "and", expr, flags=re.IGNORECASE)
    expr = re.sub(r"\bOR\b", "or", expr, flags=re.IGNORECASE)
    expr = re.sub(r"\bNOT\b", "not", expr, flags=re.IGNORECASE)

    #preserve parentheses and whitespace; process other tokens replacing IDs
    parts = re.split(r'(\(|\))', expr)
    out_parts = []
    hmms_for_func = id_to_hmms.get(func_major, {})

    for p in parts:
        if p in {'(', ')'}:
            out_parts.append(p)
            continue

        #split by whitespace but keep separators
        subparts = re.split(r'(\s+)', p)
        newsub = []
        for s in subparts:
            if s.isspace() or s == '':
                newsub.append(s)
                continue
            key = s.strip()
            #leave operator words, make sure they are lowercase (should be the case through regex earlier), and replace ID by HMM if available. 
            if key.lower() in {'and', 'or', 'not'}:
                newsub.append(key.lower())
            else:
                newsub.append(hmms_for_func.get(key, s))
        out_parts.append(''.join(newsub))

    return ''.join(out_parts).strip()


def write_outputs(cycle_rows, hmm_rows, out_cycle: Path, out_hmm: Path):
    with out_cycle.open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['Pathways', 'Hmms', 'function'])
        for r in cycle_rows:
            w.writerow(r)

    with out_hmm.open('w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        headers = [''] * 11
        headers[0] = 'function'
        headers[1] = 'Pathways'
        headers[3] = 'Gene abbreviation'
        headers[4] = 'Gene name'
        headers[5] = 'Hmm file'
        headers[6] = 'Corresponding KO'
        headers[10] = 'Hmm detecting threshold'
        w.writerow(headers)
        for r in hmm_rows:
            w.writerow(r)


def build_bipartite_graph_from_rows(rows, func_column=1, type_column=2, inputs_column=6, outputs_column=7):

    G = nx.DiGraph()

    for row in rows:
        #pad to the required number of columns
        row = [c.strip() for c in row] + [""] * max(0, max(func_column, type_column, inputs_column, outputs_column) + 1 - len(row))
        typ = (row[type_column] or "").upper()
        if 'FUNCTION' not in typ:
            continue

        func_id = (row[func_column] or '').strip()
        if not func_id:
            continue

        #add function node if not already present
        if not G.has_node(func_id):
            G.add_node(func_id, type='Function')

        #parse function inputs, add nodes and edges which are not already present
        raw_inputs = row[inputs_column] or ''
        if raw_inputs:
            for s in [x.strip() for x in raw_inputs.split(',') if x.strip()]:
                subst_id = s
                if not G.has_node(subst_id):
                    G.add_node(subst_id, type='Substance')
                if not G.has_edge(subst_id, func_id):
                    G.add_edge(subst_id, func_id)

        #parse function outputs, add nodes and edges which are not already present
        raw_outputs = row[outputs_column] or ''
        if raw_outputs:
            for s in [x.strip() for x in raw_outputs.split(',') if x.strip()]:
                subst_id = s
                if not G.has_node(subst_id):
                    G.add_node(subst_id, type='Substance')
                if not G.has_edge(func_id, subst_id):
                    G.add_edge(func_id, subst_id)

    return G


def generate_custom_db_from_tsv_one_file(custom_database_input, output_folder):
    inp = Path(custom_database_input)
    out_folder = Path(output_folder)
    if not inp.exists():
        raise FileNotFoundError(f"Input not found: {inp}")

    database_folder = out_folder / 'database'
    database_folder.mkdir(parents=True, exist_ok=True)

    rows = read_input(inp)
    cycle_rows, hmm_rows, id_to_hmms, func_to_path = parse_rows(rows)

    #translate expressions in cycle_rows using per-func mappings
    for cr in cycle_rows:
        cr[1] = translate_expression(cr[1], cr[2], id_to_hmms)

    pathway_template_path = database_folder / 'pathway_template_file.tsv'
    hmm_template_path = database_folder / 'hmm_template_file.tsv'

    write_outputs(cycle_rows, hmm_rows, pathway_template_path, hmm_template_path)

    #build bipartite graph and write as GraphML into database folder
    G = build_bipartite_graph_from_rows(rows)

    #make sure attributes are strings and add type and name label for easier graphml visualization
    for n, d in G.nodes(data=True):
        dtype = d.get('type', 'Substance')
        d['label'] = n
        for k, v in list(d.items()):
            if isinstance(v, (list, dict)):
                d[k] = str(v)

    graphml_path = database_folder / 'input_graph.graphml'
    nx.write_graphml(G, str(graphml_path), encoding='utf-8', named_key_ids=True)

    return str(hmm_template_path), str(pathway_template_path), G


if __name__ == '__main__':  #Usage: python custom_tsv_parser.py <input_tsv> [output_folder]
    if len(sys.argv) < 2:
        sys.exit(1)
    inp = Path(sys.argv[1])
    outp = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('.')
    if not inp.exists():
        print(f"Input not found: {inp}")
        sys.exit(1)
    try:
        hmm_t, path_t, G = generate_custom_db_from_tsv_one_file(inp, outp)
        print(f"Wrote: {path_t} and {hmm_t} and {outp}/database/input_graph.graphml")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
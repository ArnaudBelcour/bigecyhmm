from pyvis.network import Network
import networkx as nx

def create_html_visualisation(input_graphml_filepath, output_html_filepath):
    nx_graph = nx.read_graphml(input_graphml_filepath)
    nt = Network(height="750px", width="100%", directed=True, notebook=False, select_menu=True, filter_menu=True)
    for node in nx_graph.nodes:
        print(node, nx_graph.nodes[node])

    function_nodes = [node for node in nx_graph.nodes if nx_graph.nodes[node]['type']=='Function']
    nt.add_nodes(function_nodes, shape=['box']*len(function_nodes))
    metabolite_nodes = [node for node in nx_graph.nodes if nx_graph.nodes[node]['type']=='Metabolite']
    nt.add_nodes(metabolite_nodes, shape=['ellipse']*len(metabolite_nodes))
    nt.add_edges(nx_graph.edges)
    nt.show_buttons(filter_=['physics'])
    nt.write_html(output_html_filepath)


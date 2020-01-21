from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.util_network import NetworkUtil as util
from py2cytoscape.data.style import StyleUtil as s_util
import py2cytoscape.cytoscapejs as renderer

from IPython.display import Image
import networkx as nx
from networkx.drawing import nx_pydot as pyd
import igraph as ig
import numpy as np
import pandas as pd


class CytoscapeIntegration:
    def __init__(self, node_df, edge_df):
        self.node_df = node_df
        self.edge_df = edge_df

    def cytoscape_successful(self):
        cytoscape_successful = True

        # Create client
        cy = CyRestClient()
        # Clear current session
        cy.session.delete()
        # Create a network from edge_df
        self.edge_df.head()

        source = self.edge_df.columns[1]
        target = self.edge_df.columns[0]

        edge_graph = cy.network.create_from_dataframe(self.edge_df, source_col=source, target_col=target,
                                                      interaction_col=source, name='Edges '
                                                                                   'graph')
        # Add styles to the network
        style = cy.style.create('default')

        style_settings = {
            'NODE_FILL_COLOR': '#de1738',
            'BODE_BORDER_WIDTH': 0,
            'NODE_LABEL_COLOR': '#000000',
            'NETWORK_BACKGROUND_PAINT': '#ffffff'
        }

        style.update_defaults(style_settings)

        # Discrete mappings for specific nodes
        key_value_pair = {
            'N0': '#111E6C'
        }
        style.create_discrete_mapping(column=source, col_type='string', vp='NODE_FILL_COLOR', mappings=key_value_pair)

        cy.layout.apply(network=edge_graph)
        cy.layout.fit(network=edge_graph)
        Image(edge_graph.get_png(height=400))

        return cytoscape_successful

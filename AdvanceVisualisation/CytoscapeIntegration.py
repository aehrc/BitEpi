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

    # Method to add the specific row values of the node table columns
    def update_columns(self, edge_graph, col_name):
        node_table = edge_graph.get_node_table()
        for index, row in self.node_df.iterrows():
            cell_value = self.node_df.at[index, col_name]
            if node_table.get_value(index, col_name) != cell_value:
                print("bloop")
                node_table.set_value(index, col_name, cell_value)

        print(node_table)

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

        # Create edge network from DataFrame
        edge_graph = cy.network.create_from_dataframe(self.edge_df, source_col=source, target_col=target,
                                                      interaction_col=source, name='New network!')

        # Merge node_df in cytoscape
        edge_graph.update_node_table(df=self.node_df, network_key_col='name', data_key_col='name')
        # Add attributes to Node table
        edge_graph.create_node_column(name="Order")
        print(edge_graph.get_node_table().head())
        # Add row values of the column
        self.update_columns(edge_graph, "Order")

        cy.layout.apply(network=edge_graph)
        # Add styles to the network
        my_style = cy.style.create('my_style')

        new_styles = {
            'NODE_FILL_COLOR': 'green',
            'NODE_SIZE': 30,
            'NODE_BORDER_WIDTH': 0,
            'NODE_TRANSPARENCY': 120,
            'NODE_LABEL_COLOR': 'black',

            'EDGE_WIDTH': 3,
            'EDGE_STROKE_UNSELECTED_PAINT': '#d3d3d3',
            'EDGE_LINE_TYPE': 'SOLID',
            'EDGE_TRANSPARENCY': 120,

            'NETWORK_BACKGROUND_PAINT': 'white'
        }

        my_style.update_defaults(new_styles)

        # Discrete mappings for specific regions
        key_value_pair = {
            'N0': 'red',
            'N3': 'red',
        }

        my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_FILL_COLOR',
                                         mappings=key_value_pair)

        cy.style.apply(my_style, edge_graph)

        cy.layout.fit(network=edge_graph)
        Image(edge_graph.get_png(height=400))

        return cytoscape_successful

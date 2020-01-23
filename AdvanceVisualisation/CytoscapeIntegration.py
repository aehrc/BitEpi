import os

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
import json


class CytoscapeIntegration:
    def __init__(self, node_df, edge_df):
        self.node_df = node_df
        self.edge_df = edge_df
        self.json_file_name = 'json_file.json'
        self.json_file_path = os.path.join('OutputData/', self.json_file_name)

    # Method to convert dataframe to a json file
    def dataframe_to_json(self):

        # Create new node and edge dictionaries
        node_dic = self.node_df.to_dict('records')
        edge_dic = self.edge_df.to_dict('records')
        complete_node_list = []
        complete_edge_list = []

        # Add node_dic rows as new individual dictionaries to complete_node_list
        for i in range(len(node_dic)):
            temp_node_dic = {'data': node_dic[i], 'selected': False}
            complete_node_list.append(temp_node_dic)

        for i in range(len(edge_dic)):
            temp_edge_dic = {'data': edge_dic[i], 'selected': False}
            complete_edge_list.append(temp_edge_dic)

        # Create dictionary containing both the node and edge data
        elements_dic = {'nodes': complete_node_list, 'edges': complete_edge_list}
        name_dic = {'name': 'Node_Edge_Network'}

        col_node_list = []
        col_edge_list = []

        # Get column names and column types
        for i in range(len(self.node_df.columns)):
            temp_col_node_dic = {'columnName': self.node_df.columns[i], 'type': 'String'}
            col_node_list.append(temp_col_node_dic)

        for i in range(len(self.edge_df.columns)):
            temp_col_edge_dic = {'columnName': self.edge_df.columns[i], 'type': 'String'}
            col_edge_list.append(temp_col_edge_dic)

        # Dictionary containing all the column details
        col_types_dic = {'node': col_node_list, 'edge': col_edge_list}

        # Complete dictionary of data to be converted to json file
        full_dict = {'data': name_dic, 'elements': elements_dic}

        with open(self.json_file_path, 'w') as outfile:
            json.dump(full_dict, outfile, sort_keys=True, indent=4)

    def cytoscape_successful(self):
        cytoscape_successful = True

        # Create client
        cy = CyRestClient()
        # Clear current session
        cy.session.delete()
        # Create a network from edge_df
        self.edge_df.head()

        # Convert dataframe to json file and save file
        self.dataframe_to_json()

        # Create edge network from json file
        edge_graph = cy.network.create_from(self.json_file_path)

        cy.layout.apply(network=edge_graph)
        # Add styles to the network
        my_style = cy.style.create('my_style')

        new_styles = {
            'NODE_FILL_COLOR': 'red',
            'NODE_SIZE': 30,
            'NODE_BORDER_WIDTH': 0,
            'NODE_TRANSPARENCY': 120,
            'NODE_LABEL_COLOR': 'black',

            'EDGE_WIDTH': 3,
            'EDGE_STROKE_UNSELECTED_PAINT': '#333333',
            'EDGE_LINE_TYPE': 'SOLID',
            'EDGE_TRANSPARENCY': 120,

            'NETWORK_BACKGROUND_PAINT': 'white'
        }

        my_style.update_defaults(new_styles)

        # Discrete mappings for specific regions
        key_value_pair = {
            'N0': '#1974d2',
            'N3': '#1974d2',
            'N6': '#1974d2',
            'N19': '#1974d2'
        }

        my_style.create_discrete_mapping(column='name', col_type='String', vp='NODE_FILL_COLOR',
                                         mappings=key_value_pair)

        cy.style.apply(my_style, edge_graph)

        cy.layout.fit(network=edge_graph)
        Image(edge_graph.get_png(height=400))

        return cytoscape_successful

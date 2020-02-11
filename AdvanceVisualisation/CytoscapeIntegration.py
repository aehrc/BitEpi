# Class to send data to Cytoscape in the form of a json file
import os
from py2cytoscape.data.cyrest_client import CyRestClient
from IPython.display import Image
import json
import simplejson


class CytoscapeIntegration:
    def __init__(self, node_df, edge_df, core_details, interaction_or_edge):
        self.node_df = node_df
        self.edge_df = edge_df
        self.core_details = core_details
        self.interaction_or_edge = interaction_or_edge
        self.json_file_name = 'json_file.json'
        self.json_file_path = os.path.join('OutputData/', self.json_file_name)

    # Method to convert the DataFrames to a json object and save as a .json file
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

        # Add edge_dic rows as new individual dictionaries to complete_edge_list
        for i in range(len(edge_dic)):
            temp_edge_dic = {'data': edge_dic[i], 'selected': False}
            complete_edge_list.append(temp_edge_dic)

        # Create dictionary containing both the node and edge data
        elements_dic = {'nodes': complete_node_list, 'edges': complete_edge_list}
        name_dic = {'name': 'Node_Edge_Network'}

        # Complete dictionary of data to be converted to json file
        full_dict = {'data': name_dic, 'elements': elements_dic}

        # Write data to json file in an ordered format
        with open(self.json_file_path, 'w') as outfile:
            json.dump(full_dict, outfile, sort_keys=True, indent=4)

    # Method to create the network and add styles
    def cytoscape_successful(self, update):

        cytoscape_successful = True

        # Create client
        cy = CyRestClient()
        # Clear current session
        cy.session.delete()

        # Convert DataFrame to json file and save file
        self.dataframe_to_json()

        # Create network from json file
        node_edge_network = cy.network.create_from(self.json_file_path)

        cy.layout.apply(network=node_edge_network)

        # Add styles to the network
        my_style = cy.style.create('my_style')

        # Discrete mappings for specific regions
        order_colour_key_value_pair = {
            '1': '#c99e10',
            '2': '#9b4f0f',
            '3': '#1e434c',
            '4': '#8d230f'
        }

        edge_order_colour_key_value_pair = {
            '2': '#9b4f0f',
            '3': '#1e434c',
            '4': '#8d230f'
        }

        edge_order_size_key_value_pair = {
            '2': '5.0',
            '3': '3.0',
            '4': '1.0'
        }

        order_size_key_value_pair = {
            '1': '25.0',
            '2': '35.0',
            '3': '40.0',
            '4': '50.0'
        }

        order_shape_key_value_pair = {
            '1': 'Ellipse',
            '2': 'Diamond',
            '3': 'Triangle',
            '4': 'Hexagon'
        }

        # If the GUI is being loaded for the first time
        # Then create network with 'default' styles
        if not update:
            new_styles = {
                'NODE_FILL_COLOR': '#363636',
                'NODE_SIZE': 10,
                'NODE_BORDER_WIDTH': 0,
                'NODE_TRANSPARENCY': 255,
                'NODE_LABEL_COLOR': '#323334',

                'EDGE_WIDTH': 3,
                'EDGE_STROKE_UNSELECTED_PAINT': '#a9a9a9',
                'EDGE_LINE_TYPE': 'SOLID',
                'EDGE_TRANSPARENCY': 120,

                'NETWORK_BACKGROUND_PAINT': 'white'
            }

            my_style.update_defaults(new_styles)

            # Add these styles only if the network type is Interaction
            if self.interaction_or_edge == 1:
                my_style.create_discrete_mapping(column='order', col_type='String', vp='NODE_FILL_COLOR',
                                                 mappings=order_colour_key_value_pair)

                my_style.create_discrete_mapping(column='order', col_type='String', vp='NODE_SIZE',
                                                 mappings=order_size_key_value_pair)

                my_style.create_discrete_mapping(column='order', col_type='String', vp='NODE_SHAPE',
                                                 mappings=order_shape_key_value_pair)

            my_style.create_discrete_mapping(column='order', col_type='String',
                                             vp='EDGE_STROKE_UNSELECTED_PAINT',
                                             mappings=edge_order_colour_key_value_pair)

            my_style.create_discrete_mapping(column='order', col_type='String', vp='EDGE_WIDTH',
                                             mappings=edge_order_size_key_value_pair)

        # TODO write the code for other styles
        # TODO If need to update then read the csv file with styles and tweak it as necessary
        else:
            print('No styles specified')

        cy.style.apply(my_style, node_edge_network)

        cy.layout.fit(network=node_edge_network)
        Image(node_edge_network.get_png(height=400))

        return cytoscape_successful

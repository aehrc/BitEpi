# Controller class to join the View and the Model
from ReadWriteData import ReadWriteData
from CytoscapeIntegration import CytoscapeIntegration
from FormGUI import FormGUI


class Controller:
    def __init__(self):
        self.input_file = ''
        self.annotation_file = ''

    # Call validate_input_file of the class Readwrite to check the format of the file names
    def validate_input_files(self, input_file, annotation_file):
        read_write_data = ReadWriteData(input_file)
        valid = read_write_data.validate_input_file()
        if valid:
            self.input_file = input_file
            self.annotation_file = annotation_file
            print('The input file, {}, has been successfully validated.'
                  .format(self.input_file))
            print('The annotation file, {}, has been successfully validated.'
                  .format(self.annotation_file))
            return True
        else:
            print('Error found in input file format.')
            return False

    # Method to call forth the GUI
    def perform_form_functionality(self):
        controller = Controller()
        form = FormGUI(controller)
        form_details = form.form()

        return form_details

    # Method to call Model functionality from View
    # Takes the arguments: GUI filtering requirements, whether graph should be updated or if it's the first time,
    # If the network type is Interaction or Edge
    def perform_core_functionality(self, core_details, update, interaction_or_edge):
        # Validate input file
        read_write_data = ReadWriteData(core_details.iat[0, 0], core_details.iat[0, 1])
        valid = read_write_data.validate_input_file()
        if valid:
            print('The input file, {}, has been successfully validated.'
                  .format(core_details.iat[0, 0]))

            print('The annotation file, {}, has been successfully loaded.'
                  .format(core_details.iat[0, 1]))

            # Get the node_df, edge_df and if the functions were successful
            read_write_done = read_write_data.get_dataframes(interaction_or_edge)
            if read_write_done[2]:
                print(
                    'The input file, {}, has been successfully loaded '
                    'and the output file has been created successfully.'.format(
                        core_details.iat[0, 0]))
                print('Send data to Cytoscape.')

                # Send the DataFrames in order to create the json file and create the network
                integration = CytoscapeIntegration(read_write_done[0], read_write_done[1], core_details,
                                                   interaction_or_edge)
                # Call function and determine if cytoscape worked
                cytoscape_successful = integration.cytoscape_successful(update)

                if cytoscape_successful:
                    print('Successful creation of network!')
                else:
                    print('Network creation unsuccessful, please make sure that Cytoscape is running in the background')
            else:
                print('Error has occurred in Read and/or Write of the file.')
        else:
            print('Error found in input file format.')

        return valid

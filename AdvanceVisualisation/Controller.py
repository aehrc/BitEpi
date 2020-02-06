from ReadWriteData import ReadWriteData
from CytoscapeIntegration import CytoscapeIntegration
from FormGUI import FormGUI


class Controller:
    def __init__(self):
        self.input_file = ''
        self.annotation_file = ''

    def validate_input_files(self, input_file, annotation_file):
        read_write_data = ReadWriteData(input_file)
        valid = read_write_data.validate_input_file()
        if valid:
            self.input_file = input_file
            self.annotation_file = annotation_file
            print('The input file, {}, has been successfully validated.'
                  .format(self.input_file))
            print('The input file, {}, has been successfully validated.'
                  .format(self.annotation_file))
            return True
        else:
            print('Error found in input file format.')
            return False

    def perform_form_functionality(self):
        controller = Controller()
        form = FormGUI(controller)
        form_details = form.form()

        return form_details

    def perfrom_core_functionality(self, core_details):
        read_write_data = ReadWriteData(core_details.iat[0, 0])
        valid = read_write_data.validate_input_file()
        if valid:
            print('The input file, {}, has been successfully validated.'
                  .format(core_details.iat[0, 0]))
            read_write_done = read_write_data.read_data_from_csv()
            if read_write_done[2]:
                print(
                    'The input file, {}, has been successfully loaded '
                    'and the output file has been created successfully.'.format(
                        core_details.iat[0, 0]))
                print('Send data to Cytoscape.')

                integration = CytoscapeIntegration(read_write_done[0], read_write_done[1], core_details)
                # Call function to determine if cytoscape works
                cytoscape_successful = integration.cytoscape_successful()

                if cytoscape_successful:
                    print('Successful creation of network!')
                else:
                    print('Network creation unsuccessful, please make sure that Cytoscape is running in the background')
            else:
                print('Error has occurred in Read and/or Write of the file.')
        else:
            print('Error found in input file format.')

from ReadWriteData import ReadWriteData
from CytoscapeIntegration import CytoscapeIntegration
from FormGUI import FormGUI


class Controller:
    def __init__(self):
        print("Inside Controller")

    def validate_input_files(self, input_file, annotation_file):
        read_write_data = ReadWriteData(input_file)
        valid = read_write_data.validate_input_file()
        if valid:
            print('The input file, {}, has been successfully validated.'
                  .format(input_file))
            return True
        else:
            print('Error found in input file format.')
            return False

    def create_form(self):
        controller = Controller()
        form = FormGUI(controller)
        form_details = form.form()
        print('got form details: ', form_details)

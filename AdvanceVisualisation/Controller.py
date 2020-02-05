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

    def create_form(self):
        print('inside create form')
        controller = Controller()
        form = FormGUI(controller)
        form_details = form.form()
        print(form_details)

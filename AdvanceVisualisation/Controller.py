from ReadWriteData import ReadWriteData
from CytoscapeIntegration import CytoscapeIntegration
from FormGUI import FormGUI


class Controller:
    def __init__(self):
        print("Inside Controller")

    def validate_input_files(self, input_file, annotation_file):
        print('__controller__', input_file, '__', annotation_file)

    def create_form(self):
        controller = Controller()
        form = FormGUI(controller)
        form.form()

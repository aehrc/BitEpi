import argparse
import pandas as pd

from ReadWriteData import ReadWriteData
from CytoscapeIntegration import CytoscapeIntegration
from FormGUI import FormGUI
from Controller import Controller


def gui_main():
    controller = Controller()
    form_details = controller.perform_form_functionality()
    controller.perform_core_functionality(form_details)


def main():
    # Read in the arguments and check for validity
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '-bfile', required=True,
                        help='Argument in the form of: -i <input_file.csv>')
    args = vars(parser.parse_args())
    input_file = args['i']
    details_df = pd.DataFrame([[input_file, True]], columns=['input_file', 'reset'])
    controller = Controller()
    controller.perform_core_functionality(details_df)


if __name__ == '__main__':
    gui_main()
    # main()

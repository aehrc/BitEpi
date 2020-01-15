import argparse
import os

import pandas as pd


def validate_input_file(input_file):
    # If need to make the sub strings globally accessible
    # prefix = ''
    # input_type = ''
    first_job_index = ''
    last_job_index = ''
    order = ''
    extension = ''
    valid = True
    # Split the file name and verify each substring
    new_input_file = input_file.split(".")

    if len(new_input_file) < 1:
        valid = False

    prefix = new_input_file[0][-1]

    # Alpha or Beta
    input_type = new_input_file[1]
    if (input_type == "Alpha" or input_type == "Beta" or input_type == "best") and valid:
        if input_type != "best":
            order = new_input_file[2]
            if order != "1" and order != "2" and order != "3" and order != "4":
                valid = False

            # JobIndex is only there depending on the length of the list
            if len(new_input_file) == 5 and valid:
                valid = False

            if len(new_input_file) == 6 and valid:
                first_job_index = new_input_file[3]
                last_job_index = new_input_file[4]

        extension = new_input_file[-1]

        if (extension != "csv") and valid:
            valid = False

    # print("prefix: {} ".format(prefix))
    # print("input_type: {} ".format(input_type))
    # print("order: {} ".format(order))
    # print("first_job_index: {} ".format(first_job_index))
    # print("last_job_index: {} ".format(last_job_index))
    # print("extension: {} ".format(extension))
    # print("valid: {} ".format(valid))

    return valid


# Method to read in data from a cs file and store in a DataFrame
def read_data(input_file):
    data_read = False

    file_path = os.path.join('../sampleData/', input_file)
    df = pd.read_csv(file_path)

    print(df)
    return data_read


def main():
    # Read in the arguments and check for validity
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "-bfile", required=True, help="Argument in the form of: -i <input_file.csv>")
    args = vars(parser.parse_args())
    input_file = args["i"]
    valid = validate_input_file(input_file)
    if valid:
        print("The input file, {}, has been successfully validated.".format(input_file))
        data_read = read_data(input_file)
        if data_read:
            print("The input file, {}, has been successfully loaded.".format(input_file))
        else:
            print("Could not load data from the input file.")
    else:
        print("Error found in input file format.")


if __name__ == "__main__":
    main()

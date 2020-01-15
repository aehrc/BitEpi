import argparse
import os

import pandas as pd

input_type = ''
order = ''


# Method to validate the input file
def validate_input_file(input_file):
    # prefix = ''
    global order
    global input_type
    first_job_index = ''
    last_job_index = ''
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


# Method to create the Node DataFrame
def create_node_df(df):
    # A DataFrame with the interaction node (concat of all the names of the gene) and SNPs
    node_df = pd.DataFrame(columns=["Node"])
    int_order = int(order)

    # Loop through each row
    # Get the cell values of each column of each row and concat the values
    # Add the values to the new DataFrame
    for index, row in df.iterrows():
        for i in range(int_order):
            # Add the single SNPs at the specific position
            cell_value = df.iat[index, i + 1]
            # Add a cell value to the new DataFrame
            node_df = node_df.append(pd.DataFrame([cell_value], columns=["Node"]), ignore_index=True)
            if (i + 1) >= 2:
                new_cell_value = ''
                if (i + 1) == 2:
                    SNP_A_cell_value = df.iat[index, i]
                    new_cell_value = str(SNP_A_cell_value) + str(cell_value)

                elif (i + 1) == 3:
                    SNP_B_cell_value = df.iat[index, i]
                    SNP_A_cell_value = df.iat[index, i - 1]
                    new_cell_value = str(SNP_A_cell_value) + str(SNP_B_cell_value) + str(cell_value)

                elif (i + 1) == 4:
                    SNP_C_cell_value = df.iat[index, i]
                    SNP_B_cell_value = df.iat[index, i - 1]
                    SNP_A_cell_value = df.iat[index, i - 2]
                    new_cell_value = str(SNP_A_cell_value) + str(SNP_B_cell_value) + str(SNP_C_cell_value) + str(
                        cell_value)

                    # Add a cell value to the new DataFrame
                node_df = node_df.append(pd.DataFrame([new_cell_value], columns=["Node"]), ignore_index=True)

    # print(node_df)
    return node_df


# Method to create the edge DataFrame
def create_edge_df(df):
    # A DataFrame with all the nodes: individual -> quadruple
    edge_df = pd.DataFrame(columns=["Interaction Node"])

    # Interaction node is a concatenation of the SNPs
    # All the concatenated SNPs are connected to it
    # Later on: Check if the SNP is from Alpha and Beta
    # Check if there are any duplicates i.e same SNPs but different ordering

    return edge_df


# Method to read in data from a cs file and store in a DataFrame
def read_data(input_file):
    data_read = False

    file_path = os.path.join('../sampleData/', input_file)
    df = pd.read_csv(file_path)
    print(df)

    # Get the node_df in order to write to csv
    node_df = create_node_df(df)

    # Get the edge_df in order to write to csv
    edge_df = create_edge_df(df)

    # Each time a csv file is read and the dfs are created
    # Append to the existing csv files

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

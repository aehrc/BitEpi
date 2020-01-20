import pandas as pd
import os


class ReadWriteData:
    input_type = ''
    order = ''

    def __init__(self, input_file):
        self.input_file = input_file

    # Method to validate the input file
    @staticmethod
    def validate_input_file(self):
        # prefix = ''
        global order
        global input_type
        first_job_index = ''
        last_job_index = ''
        extension = ''
        valid = True
        # Split the file name and verify each substring
        new_input_file = ReadWriteData.input_file.split(".")

        if len(new_input_file) < 1:
            valid = False

        prefix = new_input_file[0][-1]

        # Alpha or Beta
        input_type = new_input_file[1]
        if (input_type == "Alpha" or input_type == "Beta"
            or input_type == "best") and valid:
            if input_type != "best":
                order = new_input_file[2]
                if order != "1" and order != "2" \
                        and order != "3" and order != "4":
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
    @staticmethod
    def create_node_df(self, df, int_order):
        # A DataFrame with the interaction node
        # (concat of all the names of the gene) and SNPs
        node_df = pd.DataFrame(columns=["Node"])

        # Loop through each row
        # Get the cell values of each column of each row and concat the values
        # Add the values to the new DataFrame
        for index, row in df.iterrows():
            for i in range(int_order):
                # Get the single SNPs at the specific position
                cell_value = df.iat[index, i + 1]
                # Add a cell value to the new DataFrame
                node_df = node_df.append(pd.DataFrame
                                         ([cell_value], columns=["Node"]),
                                         ignore_index=True)
                if (i + 1) >= 2:
                    new_cell_value = ''
                    if (i + 1) == 2:
                        snp_a_cell_value = df.iat[index, i]
                        new_cell_value = str(snp_a_cell_value) + str(cell_value)

                    elif (i + 1) == 3:
                        snp_b_cell_value = df.iat[index, i]
                        snp_a_cell_value = df.iat[index, i - 1]
                        new_cell_value = str(snp_a_cell_value) \
                                         + str(snp_b_cell_value) + str(cell_value)

                    elif (i + 1) == 4:
                        snp_c_cell_value = df.iat[index, i]
                        snp_b_cell_value = df.iat[index, i - 1]
                        snp_a_cell_value = df.iat[index, i - 2]
                        new_cell_value = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                         + str(snp_c_cell_value) + str(
                            cell_value)

                        # Add a cell value to the new DataFrame
                    node_df = node_df.append \
                        (pd.DataFrame([new_cell_value], columns=["Node"]),
                         ignore_index=True)

        # print(node_df)
        return node_df

    # Method to create the edge DataFrame
    @staticmethod
    def create_edge_df(self, df, int_order):
        # A DataFrame with all the nodes:
        # individual SNP -> interaction node
        edge_df = pd.DataFrame(columns=["Target", "Source"])

        for index, row in df.iterrows():
            for i in range(int_order):
                # Get the cell value
                edge_value = df.iat[index, i + 1]

                # If order is one
                # TODO Check if the SNP is from Alpha and Beta
                # If order is greater than 1
                if int_order >= 2:
                    node = ''
                    # Create the node which is a concatenation of all the SNPs
                    if int_order == 2:
                        snp_a_cell_value = df.at[index, "SNP_A"]
                        snp_b_cell_value = df.at[index, "SNP_B"]
                        node = str(snp_a_cell_value) + str(snp_b_cell_value)

                    elif int_order == 3:
                        snp_a_cell_value = df.at[index, "SNP_A"]
                        snp_b_cell_value = df.at[index, "SNP_B"]
                        snp_c_cell_value = df.at[index, "SNP_C"]
                        node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                               + str(snp_c_cell_value)

                    elif int_order == 4:
                        snp_a_cell_value = df.at[index, "SNP_A"]
                        snp_b_cell_value = df.at[index, "SNP_B"]
                        snp_c_cell_value = df.at[index, "SNP_C"]
                        snp_d_cell_value = df.at[index, "SNP_D"]
                        node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                               + str(snp_c_cell_value) + str(snp_d_cell_value)

                    # Add the new edge-node pair to the DataFrame
                    edge_df = edge_df.append(pd.DataFrame([[edge_value, node]],
                                                          columns=["Target", "Source"]),
                                             ignore_index=True)

        # print(edge_df)
        return edge_df

    # Method to check for duplicate nodes in the existing file and new DataFrame
    @staticmethod
    def check_node_duplicates(self, node_df, existing_df):
        # Create a new node DataFrame with the non-duplicated nodes
        temp_node_df = pd.concat([node_df, existing_df])
        # Remove duplicates and keep only the first occurrence of the node
        new_node_df = temp_node_df.drop_duplicates(keep='first')
        return new_node_df

    # Method to check edge duplicates
    @staticmethod
    def check_edge_duplicates(self, edge_df, existing_df):
        # Create a new edge DataFrame with the non-duplicated nodes
        temp_edge_df = pd.concat([edge_df, existing_df])
        # Remove duplicates and keep only the first occurrence of the node
        new_edge_df = temp_edge_df.drop_duplicates(subset=['Target', 'Source'], keep='first')
        return new_edge_df

    # Method to write data in the correct format to csv files
    @staticmethod
    def write_data_to_csv(self, node_df, edge_df, int_order):
        data_written_to_csv = True
        node_file_name = 'nodes.csv'
        edge_file_name = 'edges.csv'
        node_file_path = os.path.join('OutputData/', node_file_name)
        edge_file_path = os.path.join('OutputData/', edge_file_name)

        # Check if files exist
        # Create a new file if there is no existing one
        # If exists then append to the existing file
        if not os.path.isfile(node_file_path):
            print("No existing node file found. Creating a new file nodes.csv")
            node_df.to_csv(node_file_path, encoding='utf-8', index=False)
        else:
            print("Existing node file found. Appending to node.csv")
            # Read in the existing DataFrame and check for duplicates
            existing_df = pd.read_csv(node_file_path)
            # Get a new DataFrame without any duplicated nodes
            new_node_df = ReadWriteData.check_node_duplicates(node_df, existing_df)
            os.remove(node_file_path)
            new_node_df.to_csv(node_file_path, encoding='utf-8', index=False)

        if int_order != 1:
            if not os.path.isfile(edge_file_path):
                print("No existing edges file found. Creating a new file edges.csv")
                edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
            else:
                print("Existing edges file found. Appending to edges.csv")
                # Read in the existing DataFrame and check for duplicates
                existing_df = pd.read_csv(edge_file_path)
                # Get a new DataFrame without any duplicated edges
                new_edge_df = ReadWriteData.check_edge_duplicates(edge_df, existing_df)
                os.remove(edge_file_path)
                new_edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)

        return data_written_to_csv

    # Method to read in data and write data from and to a csv file
    # TODO rename this function
    @staticmethod
    def read_data_from_csv(self):
        read_write_done = True

        file_path = os.path.join('../sampleData/', ReadWriteData.input_file)
        df = pd.read_csv(file_path)
        int_order = int(order)

        if df.empty:
            read_write_done = False
        else:
            # print(df)
            # Get the node_df in order to write to csv
            node_df = ReadWriteData.create_node_df(df, int_order)
            # Get the edge_df in order to write to csv
            # Ignore order 1 as the SNPs are already added to the nodes
            edge_df = pd.DataFrame()
            if int_order != 1:
                edge_df = ReadWriteData.create_edge_df(df, int_order)

            if node_df.empty or (edge_df.empty and int_order != 1):
                read_write_done = False
                print("Error the newly created DataFrames are empty.")
            else:
                data_written = ReadWriteData.write_data_to_csv(node_df, edge_df, int_order)

                if not data_written:
                    read_write_done = False
                    print("Error could not write data to the csv file/s!")

        return read_write_done

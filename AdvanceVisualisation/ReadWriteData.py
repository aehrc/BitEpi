import os
import pandas as pd


class ReadWriteData:
    input_type = ''
    order = ''

    def __init__(self, input_file):
        self.input_file = input_file
        self.int_id = 0

    # Method to validate the input file
    def validate_input_file(self):
        # prefix = ''
        global order
        global input_type
        first_job_index = ''
        last_job_index = ''
        extension = ''
        valid = True
        # Split the file name and verify each substring
        new_input_file = self.input_file.split('.')

        if len(new_input_file) < 1:
            valid = False

        prefix = new_input_file[0][-1]

        # Alpha or Beta
        input_type = new_input_file[1]
        if (input_type == 'Alpha' or input_type == 'Beta' or input_type == 'best') \
                and valid:
            if input_type != 'best':
                order = new_input_file[2]
                if order != '1' and order != '2' \
                        and order != '3' and order != '4':
                    valid = False

                # JobIndex is only there depending on the length of the list
                if len(new_input_file) == 5 and valid:
                    valid = False

                if len(new_input_file) == 6 and valid:
                    first_job_index = new_input_file[3]
                    last_job_index = new_input_file[4]

            extension = new_input_file[-1]

            if (extension != 'csv') and valid:
                valid = False

        return valid

    # Method to create the node DataFrame
    def create_node_df(self, df, int_order):
        # A DataFrame with the interaction node
        # (concat of all the names of the gene) and SNPs
        node_df = pd.DataFrame(columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist'])

        # Determine the input type
        new_alpha = ''
        new_beta = ''
        if input_type == 'Alpha':
            new_alpha = 'Alpha'
            new_beta = 'None'
        else:
            new_alpha = 'None'
            new_beta = 'Beta'

        # Loop through each row
        # Get the cell values of each column of each row and concat the values
        # Add the values to the new DataFrame
        for index, row in df.iterrows():
            for i in range(int_order):
                # Get the single SNPs at the specific position
                cell_value = df.iat[index, i + 1]
                new_order = '1'

                # Determine if the node comes from order=1 or an interaction node
                reason = ''
                if int_order == 1:
                    reason = 'Important'
                else:
                    reason = 'Presentation'

                # Add a cell value to the new DataFrame
                node_df = node_df.append(pd.DataFrame
                                         ([[cell_value, cell_value, new_order, new_alpha, new_beta, reason]],
                                          columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                         ignore_index=True)
            if int_order >= 2:
                new_cell_value = ''
                new_order = ''
                if int_order == 2:
                    snp_b_cell_value = df.at[index, 'SNP_B']
                    snp_a_cell_value = df.at[index, 'SNP_A']
                    new_order = '2'
                    new_cell_value = str(snp_a_cell_value) + str(snp_b_cell_value)

                elif int_order == 3:
                    snp_c_cell_value = df.at[index, 'SNP_C']
                    snp_b_cell_value = df.at[index, 'SNP_B']
                    snp_a_cell_value = df.at[index, 'SNP_A']
                    new_order = '3'
                    new_cell_value = str(snp_a_cell_value) \
                                     + str(snp_b_cell_value) + str(snp_c_cell_value)

                elif int_order == 4:
                    snp_d_cell_value = df.at[index, 'SNP_D']
                    snp_c_cell_value = df.at[index, 'SNP_C']
                    snp_b_cell_value = df.at[index, 'SNP_B']
                    snp_a_cell_value = df.at[index, 'SNP_A']
                    new_order = '4'
                    new_cell_value = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                     + str(snp_c_cell_value) + str(snp_d_cell_value)

                    # Add a cell value to the new DataFrame
                node_df = node_df.append \
                        (pd.DataFrame(
                        [[new_cell_value, new_cell_value, new_order, new_alpha, new_beta, reason]],
                        columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                        ignore_index=True)

        # print(node_df)
        return node_df

    # Method to create the edge DataFrame
    def create_edge_df(self, df, int_order):
        # A DataFrame with all the nodes:
        # individual SNP -> interaction node
        edge_df = pd.DataFrame(columns=['id', 'source', 'target', 'interaction'])

        # If the order is greater than 1 then there exist interactions between nodes
        if int_order >= 2:
            for index, row in df.iterrows():
                for i in range(int_order):
                    # Get the cell value
                    edge_value = df.iat[index, i + 1]

                    # If order is greater than 1
                    if int_order >= 2:
                        node = ''
                        # Create the node which is a concatenation of all the SNPs
                        if int_order == 2:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value)

                        elif int_order == 3:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            snp_c_cell_value = df.at[index, 'SNP_C']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                   + str(snp_c_cell_value)

                        elif int_order == 4:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            snp_c_cell_value = df.at[index, 'SNP_C']
                            snp_d_cell_value = df.at[index, 'SNP_D']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                   + str(snp_c_cell_value) + str(snp_d_cell_value)

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[edge_value, node,
                                                                edge_value,
                                                                node]],
                                                              columns=['id', 'source', 'target', 'interaction']),
                                                 ignore_index=True)

        # print(edge_df)
        return edge_df

    # Method to check for duplicate nodes in the existing file and new DataFrame
    def check_node_duplicates(self, node_df, existing_df):

        # out2.Alpha.3.csv
        # out2.Beta.2.csv
        # out2.Alpha.1.csv

        # Create a new node DataFrame with the non-duplicated nodes
        new_node_df = pd.concat([node_df, existing_df], axis=0).reset_index(drop=True)
        # Remove duplicates and keep only the node from order=1 or the original node
        for index, row in new_node_df.iterrows():
            # Get all the duplicates of a specific id at an index
            # Filter all the rows which has that id
            # Loop through them and look for the specific conditions and drop others

            temp_node = new_node_df.at[index, 'id']
            temp_df = new_node_df.loc[new_node_df['id'] == temp_node]
            found = False
            if not found:
                temp_node_important = temp_df.loc[temp_df['reason to exist'] == 'Important']
                # Not important
                if temp_node_important.empty:
                    found_both_alpha_and_beta = temp_df.loc[(temp_df['Alpha'] == 'Alpha')
                                                            & (temp_df['Beta'] == 'Beta')]
                    # Both Alpha and Beta together aren't present
                    if found_both_alpha_and_beta.empty:
                        found_only_alpha = temp_df.loc[temp_df['Alpha'] == 'Alpha']
                        found_only_beta = temp_df.loc[temp_df['Beta'] == 'Beta']

                        # There exists rows with either one of Alpha or Beta values
                        if (not found_only_alpha.empty) and (not found_only_beta.empty):
                            # Get the first alpha row and the first beta row
                            beta_value = found_only_beta.at[0, 'Beta']
                            # Add the Beta row to the Alpha row
                            found_only_alpha.iloc[0, found_only_alpha.columns.get_loc('Beta')] = beta_value
                            # Create a new row and delete others
                            new_cell = found_only_alpha.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # If only Alpha values were present
                        elif (not found_only_alpha.empty) and found_only_beta.empty:
                            # Create a new row and delete others
                            new_cell = found_only_alpha.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # If only Beta values were found
                        elif found_only_alpha.empty and (not found_only_beta.empty):
                            # Create a new row and delete others
                            new_cell = found_only_beta.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # Both aren't present
                        else:
                            new_cell = temp_df.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True

                    # Both Alpha and Beta are present for Presentation cell
                    else:
                        new_cell = found_both_alpha_and_beta.head(1)
                        new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                        new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                        found = True

                # Important
                else:
                    found_both_alpha_and_beta = temp_df.loc[(temp_df['Alpha'] == 'Alpha')
                                                            & (temp_df['Beta'] == 'Beta')]
                    # Both Alpha and Beta together aren't present
                    if found_both_alpha_and_beta.empty:
                        found_only_alpha = temp_df.loc[temp_df['Alpha'] == 'Alpha']
                        found_only_beta = temp_df.loc[temp_df['Beta'] == 'Beta']

                        # There exists rows with either one of Alpha or Beta values
                        if (not found_only_alpha.empty) and (not found_only_beta.empty):
                            # Get the first alpha row and the first beta row
                            beta_value = found_only_beta.at[0, 'Beta']
                            # Add the Beta row to the Alpha row
                            found_only_alpha.iloc[0, found_only_alpha.columns.get_loc('Beta')] = beta_value
                            # Create a new row and delete others
                            new_cell = found_only_alpha.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # If only Alpha values were present
                        elif (not found_only_alpha.empty) and found_only_beta.empty:
                            # Create a new row and delete others
                            new_cell = found_only_alpha.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # If only Beta values were found
                        elif found_only_alpha.empty and (not found_only_beta.empty):
                            # Create a new row and delete others
                            new_cell = found_only_beta.head(1)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True
                        # Both aren't present
                        else:
                            new_cell = temp_df.head(0)
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                            found = True

                        # Both Alpha and Beta are present for Presentation cell
                    else:
                        new_cell = found_both_alpha_and_beta.head(1)
                        new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                        new_node_df = pd.concat([new_node_df, new_cell], axis=0)
                        found = True
                print(new_node_df)
                print(len(new_node_df))
                print(index)

                if index == len(new_node_df) - 1:
                    print('End')
                    break
                print('_____________________________________________')
            # Major error coz important and other values aren't found
            else:
                print('Abort mission: Duplicates found!')

        print(new_node_df)
        return new_node_df

    # Method to check edge duplicates
    def check_edge_duplicates(self, edge_df, existing_df):
        # Create a new edge DataFrame with the non-duplicated nodes
        temp_edge_df = pd.concat([edge_df, existing_df], axis=0)
        # Remove duplicates and keep only the first occurrence of the node
        new_edge_df = temp_edge_df.drop_duplicates(subset=['target', 'source', 'id'], keep='first')
        return new_edge_df

    # Method to write data in the correct format to csv files
    def write_data_to_csv(self, node_df, edge_df, int_order):
        data_written_to_csv = True
        node_file_name = 'nodes.csv'
        edge_file_name = 'edges.csv'
        node_file_path = os.path.join('OutputData/', node_file_name)
        edge_file_path = os.path.join('OutputData/', edge_file_name)

        # DataFrames to send to Cytoscape
        correct_node_df = pd.DataFrame
        correct_edge_df = pd.DataFrame

        # Check if files exist
        # Create a new file if there is no existing one
        # If exists then append to the existing file
        if not os.path.isfile(node_file_path):
            print('No existing node file found. Creating a new file nodes.csv')
            node_df.to_csv(node_file_path, encoding='utf-8', index=False)
            correct_node_df = node_df
        else:
            print('Existing node file found. Appending to node.csv')
            # Read in the existing DataFrame and check for duplicates
            existing_df = pd.read_csv(node_file_path)
            # Get a new DataFrame without any duplicated nodes
            new_node_df = self.check_node_duplicates(node_df, existing_df)
            os.remove(node_file_path)
            new_node_df.to_csv(node_file_path, encoding='utf-8', index=False)
            correct_node_df = new_node_df

        if not os.path.isfile(edge_file_path):
            print('No existing edges file found. Creating a new file edges.csv')
            edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
            correct_edge_df = edge_df
        else:
            print('Existing edges file found. Appending to edges.csv')
            # Read in the existing DataFrame and check for duplicates
            existing_df = pd.read_csv(edge_file_path)
            if existing_df.empty:
                print('The existing file is empty. Creating a new file edges.csv')
                edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
                correct_edge_df = edge_df
            else:
                # Get a new DataFrame without any duplicated edges
                new_edge_df = self.check_edge_duplicates(edge_df, existing_df)
                os.remove(edge_file_path)
                new_edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
                correct_edge_df = new_edge_df

        # Send the DataFrames (which were checked for duplication) to Cytoscape
        data_written = [correct_node_df, correct_edge_df, data_written_to_csv]

        return data_written

    # Method to read in data and write data from and to a csv file
    # TODO rename this function
    def read_data_from_csv(self):
        read_write_done = True

        file_path = os.path.join('../sampleData/', self.input_file)
        df = pd.read_csv(file_path)
        int_order = int(order)

        if df.empty:
            read_write_done = False
        else:
            # print(df)
            # Get the node_df in order to write to csv
            node_df = self.create_node_df(df, int_order)
            # Get the edge_df in order to write to csv
            edge_df = self.create_edge_df(df, int_order)

            if node_df.empty or (edge_df.empty and int_order != 1):
                read_write_done = False
                print('Error the newly created DataFrames are empty.')
            else:
                data_written = self.write_data_to_csv(node_df, edge_df, int_order)

                if not data_written[2]:
                    read_write_done = False
                    print('Error could not write data to the csv file/s!')

                else:
                    # Send the DataFrames as an array
                    cytoscape_df = [data_written[0], data_written[1], read_write_done]

        return cytoscape_df

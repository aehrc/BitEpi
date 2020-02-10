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
        global order
        global input_type
        first_job_index = ''
        last_job_index = ''
        extension = ''
        valid = True
        # Split the file name and verify each substring
        new_input_file = self.input_file.split('.')
        if len(new_input_file) <= 1:
            valid = False
        else:
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
    def create_interaction_node_df(self, df, int_order):
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
            # Get the single SNPs at the specific position
            snp_a_cell_value = df.at[index, 'SNP_A']
            new_order = '1'

            # Determine if the node comes from order=1 or an interaction node
            reason = ''
            if int_order == 1:
                reason = 'Important'
            else:
                reason = 'Presentation'

                # Add a cell value to the new DataFrame
            node_df = node_df.append(pd.DataFrame
                                     ([[snp_a_cell_value, snp_a_cell_value, new_order, new_alpha, new_beta, reason]],
                                      columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                     ignore_index=True)
            new_cell_value = ''
            new_order = ''
            if 'SNP_B' in df.columns:
                snp_b_cell_value = df.at[index, 'SNP_B']
                new_order = '2'
                # Add a cell value to the new DataFrame
                node_df = node_df.append \
                        (pd.DataFrame(
                        [[snp_b_cell_value, snp_b_cell_value, new_order, new_alpha, new_beta, reason]],
                        columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                        ignore_index=True)

            if 'SNP_C' in df.columns:
                snp_c_cell_value = df.at[index, 'SNP_C']
                new_order = '3'
                # Add a cell value to the new DataFrame
                node_df = node_df.append \
                        (pd.DataFrame(
                        [[snp_c_cell_value, snp_c_cell_value, new_order, new_alpha, new_beta, reason]],
                        columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                        ignore_index=True)

            if 'SNP_D' in df.columns:
                snp_d_cell_value = df.at[index, 'SNP_D']
                new_order = '4'
                # Add a cell value to the new DataFrame
                node_df = node_df.append \
                        (pd.DataFrame(
                        [[snp_d_cell_value, snp_d_cell_value, new_order, new_alpha, new_beta, reason]],
                        columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                        ignore_index=True)

        # print(node_df)
        return node_df

    # Method to create the edge DataFrame for interaction node mode
    def create_interaction_edge_df(self, df, int_order):
        # A DataFrame with all the nodes:
        # individual SNP -> interaction node
        edge_df = pd.DataFrame(columns=['id', 'source', 'target', 'interaction', 'order'])

        # If the order is greater than 1 then there exist interactions between nodes
        if int_order >= 2:
            for index, row in df.iterrows():
                for i in range(int_order):
                    # Get the cell value
                    edge_value = df.iat[index, i + 1]

                    # If order is greater than 1
                    if int_order >= 2:
                        node = ''
                        new_order = ''
                        # Create the node which is a concatenation of all the SNPs
                        if int_order == 2:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value)
                            new_order = '2'

                        elif int_order == 3:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            snp_c_cell_value = df.at[index, 'SNP_C']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                   + str(snp_c_cell_value)
                            new_order = '3'

                        elif int_order == 4:
                            snp_a_cell_value = df.at[index, 'SNP_A']
                            snp_b_cell_value = df.at[index, 'SNP_B']
                            snp_c_cell_value = df.at[index, 'SNP_C']
                            snp_d_cell_value = df.at[index, 'SNP_D']
                            node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                                   + str(snp_c_cell_value) + str(snp_d_cell_value)
                            new_order = '4'

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[edge_value, node,
                                                                edge_value,
                                                                node, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'order']),
                                                 ignore_index=True)

        # print(edge_df)
        return edge_df

    # Method to create the edge DataFrame
    def create_edge_df(self, df, int_order):
        # A DataFrame with all the nodes:
        # individual SNP -> interaction node
        edge_df = pd.DataFrame(
            columns=['id', 'source', 'target', 'interaction', 'edge label', 'Alpha', 'Beta', 'order'])
        # Determine the input type
        new_alpha = ''
        new_beta = ''
        if input_type == 'Alpha':
            new_alpha = 'Alpha'
            new_beta = 'None'
        else:
            new_alpha = 'None'
            new_beta = 'Beta'

        # If the order is greater than 1 then there exist interactions between nodes
        if int_order >= 2:
            for index, row in df.iterrows():
                # If order is greater than 1
                if int_order >= 2:
                    node = ''
                    new_order = ''
                    # Create the node which is a concatenation of all the SNPs
                    if int_order == 2:
                        snp_a_cell_value = df.at[index, 'SNP_A']
                        snp_b_cell_value = df.at[index, 'SNP_B']
                        node = str(snp_a_cell_value) + str(snp_b_cell_value)
                        new_order = '2'

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_a_cell_value,
                                                                snp_b_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

                        edge_df = edge_df.append(pd.DataFrame([[node, snp_b_cell_value,
                                                                snp_a_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

                    elif int_order == 3:
                        snp_a_cell_value = df.at[index, 'SNP_A']
                        snp_b_cell_value = df.at[index, 'SNP_B']
                        snp_c_cell_value = df.at[index, 'SNP_C']
                        node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                               + str(snp_c_cell_value)
                        new_order = '3'

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_a_cell_value,
                                                                snp_b_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_b_cell_value,
                                                                snp_c_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_a_cell_value,
                                                                snp_c_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

                    elif int_order == 4:
                        snp_a_cell_value = df.at[index, 'SNP_A']
                        snp_b_cell_value = df.at[index, 'SNP_B']
                        snp_c_cell_value = df.at[index, 'SNP_C']
                        snp_d_cell_value = df.at[index, 'SNP_D']
                        node = str(snp_a_cell_value) + str(snp_b_cell_value) \
                               + str(snp_c_cell_value) + str(snp_d_cell_value)
                        new_order = '4'

                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_a_cell_value,
                                                                snp_b_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)
                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_b_cell_value,
                                                                snp_c_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)
                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_c_cell_value,
                                                                snp_d_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)
                        # Add the new edge-node pair to the DataFrame
                        edge_df = edge_df.append(pd.DataFrame([[node, snp_d_cell_value,
                                                                snp_a_cell_value,
                                                                node, node, new_alpha, new_beta, new_order]],
                                                              columns=['id', 'source', 'target', 'interaction',
                                                                       'edge label', 'Alpha', 'Beta', 'order']),
                                                 ignore_index=True)

        # print(edge_df)
        return edge_df

    # Method to get the number of interaction nodes connected to a node
    def create_connection_count_df(self, correct_edge_df):
        connection_count_df = pd.DataFrame(columns=['id', 'count'])
        for index, row in correct_edge_df.iterrows():
            temp_node = correct_edge_df.at[index, 'target']
            # TODO Lengths must match to compare??
            count = str(correct_edge_df.loc[correct_edge_df.target == temp_node, 'target'].count())
            connection_count_df = connection_count_df.append(pd.DataFrame([[temp_node, count]]
                                                                          , columns=['id', 'count'])
                                                             , ignore_index=True)
        connection_count_df = connection_count_df.drop_duplicates(subset=['id'], keep='first')
        return connection_count_df

    # Method to check for duplicate nodes in the existing file and new DataFrame
    def check_node_duplicates(self, node_df, existing_df):
        # Delete the count column of the existing DataFrame as it will be merged later
        if 'count' in existing_df.columns:
            existing_df = existing_df.drop('count', axis=1)
        # Create a new node DataFrame with the non-duplicated nodes
        new_node_df = pd.concat([node_df, existing_df], axis=0).reset_index(drop=True)
        # Remove duplicates and keep only the node from order=1 or the original node
        i = len(new_node_df) - 1
        for index, row in new_node_df.iterrows():
            # Get all the duplicates of a specific id at an index
            # Filter all the rows which has that id
            # Loop through them and look for the specific conditions and drop others
            temp_node = new_node_df.at[i, 'id']
            new_order = ''
            temp_df = new_node_df.loc[new_node_df['id'] == temp_node]
            found = False

            if not found:
                temp_node_important = temp_df.loc[temp_df['reason to exist'] == 'Important'].reset_index(drop=True)
                temp_node_presentation = temp_df.loc[temp_df['reason to exist'] == 'Presentation'].reset_index(
                    drop=True)

                # Not important i.e presentation
                if temp_node_important.empty:
                    found_both_alpha_and_beta = temp_node_presentation.loc[(temp_node_presentation['Alpha'] == 'Alpha')
                                                                           & (temp_node_presentation['Beta'] == 'Beta')] \
                        .reset_index(drop=True)
                    new_reason = 'Presentation'
                    # Both Alpha and Beta together aren't present
                    if found_both_alpha_and_beta.empty:
                        found_only_alpha = temp_df.loc[(temp_df['Alpha'] == 'Alpha')
                                                       & (temp_df['reason to exist'] == 'Presentation')].reset_index(
                            drop=True)
                        found_only_beta = temp_df.loc[(temp_df['Beta'] == 'Beta')
                                                      & (temp_df['reason to exist'] == 'Presentation')].reset_index(
                            drop=True)

                        # There exists rows with either one of Alpha or Beta values
                        if (not found_only_alpha.empty) and (not found_only_beta.empty):
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_alpha.at[0, 'order'])
                            new_alpha = found_only_alpha.at[0, 'Alpha']
                            # Get the beta value
                            new_beta = found_only_beta.at[0, 'Beta']

                            # Append the new obj to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, new_alpha, new_beta, new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True
                        # If only Alpha values were present
                        elif (not found_only_alpha.empty) and found_only_beta.empty:
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_alpha.at[0, 'order'])
                            new_alpha = found_only_alpha.at[0, 'Alpha']

                            # Append the new obj to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, new_alpha, 'None', new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True
                        # If only Beta values were found
                        elif found_only_alpha.empty and (not found_only_beta.empty):
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_beta.at[0, 'order'])
                            new_beta = found_only_beta.at[0, 'Beta']
                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, 'None', new_beta, new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)
                            i = i - len(temp_df)
                            found = True
                        # Both aren't present
                        else:
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(temp_node_presentation.at[0, 'order'])
                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, 'None', 'None', new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True

                    # Both Alpha and Beta are present for Presentation cell
                    else:
                        # Remove duplicates
                        new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                        # Create a new DataFrame object
                        # Convert int order value to str
                        new_order = str(found_both_alpha_and_beta.at[0, 'order'])
                        new_beta = found_both_alpha_and_beta.at[0, 'Beta']
                        new_alpha = found_both_alpha_and_beta.at[0, 'Alpha']
                        # Append the new cell to the DataFrame
                        new_node_df = new_node_df.append \
                                (pd.DataFrame(
                                [[temp_node, temp_node, new_order, new_alpha, new_beta, new_reason]],
                                columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                ignore_index=True).reset_index(drop=True)

                        i = i - len(temp_df)
                        found = True

                # Important
                else:
                    found_both_alpha_and_beta = temp_node_important.loc[(temp_node_important['Alpha'] == 'Alpha')
                                                                        & (temp_node_important['Beta'] == 'Beta')] \
                        .reset_index(drop=True)
                    new_reason = 'Important'

                    # Both Alpha and Beta together aren't present
                    if found_both_alpha_and_beta.empty:
                        found_only_alpha = temp_df.loc[(temp_df['Alpha'] == 'Alpha')
                                                       & (temp_df['reason to exist'] == 'Important')].reset_index(
                            drop=True).reset_index(drop=True)
                        found_only_beta = temp_df.loc[(temp_df['Beta'] == 'Beta')
                                                      & (temp_df['reason to exist'] == 'Important')].reset_index(
                            drop=True).reset_index(drop=True)

                        # There exists rows with either one of Alpha or Beta values
                        if (not found_only_alpha.empty) and (not found_only_beta.empty):
                            # Get the first alpha row and the first beta row
                            new_beta = found_only_beta.at[0, 'Beta']
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_alpha.at[0, 'order'])
                            new_alpha = found_only_alpha.at[0, 'Alpha']

                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, new_alpha, new_beta, new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True
                        # If only Alpha values were present
                        elif (not found_only_alpha.empty) and found_only_beta.empty:
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_alpha.at[0, 'order'])
                            new_alpha = found_only_alpha.at[0, 'Alpha']

                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, new_alpha, 'None', new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True
                        # If only Beta values were found
                        elif found_only_alpha.empty and (not found_only_beta.empty):
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(found_only_beta.at[0, 'order'])
                            new_beta = found_only_beta.at[0, 'Beta']
                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, 'None', new_beta, new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)
                            i = i - len(temp_df)
                            found = True
                        # Both aren't present
                        else:
                            # Remove duplicates
                            new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                            # Create a new DataFrame object
                            # Convert int order value to str
                            new_order = str(temp_node_presentation.at[0, 'order'])
                            # Append the new cell to the DataFrame
                            new_node_df = new_node_df.append \
                                    (pd.DataFrame(
                                    [[temp_node, temp_node, new_order, 'None', 'None', new_reason]],
                                    columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                    ignore_index=True).reset_index(drop=True)

                            i = i - len(temp_df)
                            found = True

                    # Both Alpha and Beta are present for Important cell
                    else:
                        # Remove duplicates
                        new_node_df = new_node_df.loc[new_node_df['id'] != temp_node]
                        # Create a new DataFrame object
                        # Convert int order value to str
                        new_order = str(found_both_alpha_and_beta.at[0, 'order'])
                        new_beta = found_both_alpha_and_beta.at[0, 'Beta']
                        new_alpha = found_both_alpha_and_beta.at[0, 'Alpha']
                        # Append the new cell to the DataFrame
                        new_node_df = new_node_df.append \
                                (pd.DataFrame(
                                [[temp_node, temp_node, new_order, new_alpha, new_beta, new_reason]],
                                columns=['id', 'name', 'order', 'Alpha', 'Beta', 'reason to exist']),
                                ignore_index=True).reset_index(drop=True)

                        i = i - len(temp_df)
                        found = True
                # exit for loop once the correct new_node_df is made
                if index == len(new_node_df) - 1:
                    break
            # Major error coz important and other values aren't found
            else:
                print('Abort mission: Duplicates found!')

        # print(new_node_df)
        # print('_______________________________________________________________________________________')
        return new_node_df

    # Method to check edge duplicates
    def check_edge_duplicates(self, edge_df, existing_df, interaction_or_edge):
        # Create a new edge DataFrame with the non-duplicated nodes
        temp_edge_df = pd.concat([edge_df, existing_df], axis=0)
        # Remove duplicates and keep only the first occurrence of the node
        if interaction_or_edge == 1:
            new_edge_df = temp_edge_df.drop_duplicates(subset=['target', 'source', 'id'], keep='first')
        elif interaction_or_edge == 2:
            new_edge_df = temp_edge_df.drop_duplicates(subset=['target', 'source', 'id', 'edge label'], keep='first') \
                .reset_index(drop=True)

            i = len(new_edge_df) - 1
            for index, row in new_edge_df.iterrows():
                temp_node = new_edge_df.at[i, 'target']
                filtered_df = new_edge_df.loc[((new_edge_df['target'] == temp_node)
                                               & (new_edge_df['source'] == temp_node))] \
                    .reset_index(drop=True)

                if not filtered_df.empty:
                    # Remove duplicates
                    new_edge_df = new_edge_df.loc[((new_edge_df['target'] != temp_node)
                                                   & (new_edge_df['source'] != temp_node))].reset_index(drop=True)
                    i = len(new_edge_df) - 1
                else:
                    i = i - 1

                # exit for loop once the correct new_node_df is made
                if i == 0:
                    break

        # print(new_edge_df)
        return new_edge_df.reset_index(drop=True)

    def get_merged_new_node_df(self, new_node_df, connection_count_df):
        if order != '1':
            new_node_df = new_node_df.merge(connection_count_df, how='left')
            new_node_df['count'].fillna('None', inplace=True)

        return new_node_df

    # Method to write data in the correct format to csv files
    def write_data_to_csv(self, node_df, edge_df, int_order, interaction_or_edge):
        data_written_to_csv = True
        node_file_name = 'nodes.csv'
        edge_file_name = 'edges.csv'
        trans_edge_file_name = 'trans_edges.csv'
        trans_node_file_name = 'trans_nodes.csv'

        if interaction_or_edge == 1:
            node_file_path = os.path.join('OutputData/', node_file_name)
            edge_file_path = os.path.join('OutputData/', edge_file_name)
        elif interaction_or_edge == 2:
            node_file_path = os.path.join('OutputData/', trans_node_file_name)
            edge_file_path = os.path.join('OutputData/', trans_edge_file_name)

        # DataFrames to send to Cytoscape
        correct_node_df = pd.DataFrame
        correct_edge_df = pd.DataFrame

        # Check if files exist
        # Create a new file if there is no existing one
        # If exists then append to the existing file
        if not os.path.isfile(edge_file_path):
            print('No existing edges file found. Creating a new file edges.csv')
            new_edge_df = self.check_edge_duplicates(edge_df, edge_df, interaction_or_edge)
            new_edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
            correct_edge_df = new_edge_df
        else:
            print('Existing edges file found. Appending to edges.csv')
            # Read in the existing DataFrame and check for duplicates
            existing_df = pd.read_csv(edge_file_path)
            if existing_df.empty:
                print('The existing file is empty. Creating a new file edges.csv')
                new_edge_df = self.check_edge_duplicates(edge_df, edge_df, interaction_or_edge)
                new_edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
                correct_edge_df = new_edge_df
            else:
                # Get a new DataFrame without any duplicated edges
                new_edge_df = self.check_edge_duplicates(edge_df, existing_df, interaction_or_edge)
                os.remove(edge_file_path)
                new_edge_df.to_csv(edge_file_path, encoding='utf-8', index=False)
                correct_edge_df = new_edge_df

        connection_count_df = self.create_connection_count_df(correct_edge_df)

        if not os.path.isfile(node_file_path):
            print('No existing node file found. Creating a new file nodes.csv')
            new_node_df = self.check_node_duplicates(node_df, node_df)
            new_node_df = self.get_merged_new_node_df(new_node_df, connection_count_df)
            new_node_df.to_csv(node_file_path, encoding='utf-8', index=False)
            correct_node_df = new_node_df
        else:
            print('Existing node file found. Appending to node.csv')
            # Read in the existing DataFrame and check for duplicates
            existing_df = pd.read_csv(node_file_path)
            # Get a new DataFrame without any duplicated nodes
            new_node_df = self.check_node_duplicates(node_df, existing_df)
            new_node_df = self.get_merged_new_node_df(new_node_df, connection_count_df)
            os.remove(node_file_path)
            new_node_df.to_csv(node_file_path, encoding='utf-8', index=False)
            correct_node_df = new_node_df

        # Send the DataFrames (which were checked for duplication) to Cytoscape
        data_written = [correct_node_df, correct_edge_df, data_written_to_csv]

        return data_written

    # Method to read in data and write data from and to a csv file
    def read_data_from_csv(self, interaction_or_edge):
        read_write_done = True

        file_path = os.path.join('../sampleData/', self.input_file)
        df = pd.read_csv(file_path)
        int_order = int(order)

        if df.empty:
            read_write_done = False
        else:
            # print(df)
            # Get the node_df in order to write to csv
            # Get the edge_df in order to write to csv
            if interaction_or_edge == 1:
                node_df = self.create_interaction_node_df(df, int_order)
                edge_df = self.create_interaction_edge_df(df, int_order)
            elif interaction_or_edge == 2:
                node_df = self.create_node_df(df, int_order)
                edge_df = self.create_edge_df(df, int_order)

            if node_df.empty or (edge_df.empty and int_order != 1):
                read_write_done = False
                print('Error the newly created DataFrames are empty.')
            else:
                data_written = self.write_data_to_csv(node_df, edge_df, int_order, interaction_or_edge)
                if not data_written[2]:
                    read_write_done = False
                    print('Error could not write data to the csv file/s!')

                else:
                    # Send the DataFrames as an array
                    cytoscape_df = [data_written[0], data_written[1], read_write_done]

        return cytoscape_df

# Class containing all the GUI functionality
import Tkinter as tk
import tkMessageBox
import pandas as pd


class FormGUI:

    def __init__(self, controller):
        self.controller = controller
        self.hide_bool = False
        self.show_bool = False
        self.highlight_bool = False
        self.gray_bool = False
        self.reset_bool = False

    # Load file upon clicking submit on the GUI
    def load_files(self, input_file, annotation_file, interaction_or_edge):
        details_df = pd.DataFrame([[input_file, annotation_file, True]],
                                  columns=['input_file', 'annotation_file', 'reset'])
        # Call the validate function in Controller class
        valid = self.controller.perform_core_functionality(details_df, False, interaction_or_edge)
        if valid:
            tkMessageBox.showinfo('Success', 'Check Cytoscape to view your network!')
        else:
            tkMessageBox.showinfo('Error', 'Please input valid files')

    # Method to create the form
    def form(self):
        root = tk.Tk()

        # Main window
        canvas = tk.Canvas(root, bg='#763626', height=640, width=640)
        canvas.pack()

        main_title = tk.Label(root, bg='#763626', text='Hi, please fill out the specifics you would like to see!')
        main_title.place(relx=0.1, rely=0.05, relheight=0.05, relwidth=0.75)

        # Frame to input files and load
        file_frame = tk.Frame(root, bg='#2A3132', bd=5)
        file_frame.place(relx=0.5, rely=0.1, relwidth=0.95, relheight=0.2, anchor='n')

        input_file_title = tk.Label(file_frame, text='Prefix file path: ')
        input_file_title.place(relx=0.02, rely=0.1, relheight=0.2, relwidth=0.25)

        input_file_entry = tk.Entry(file_frame)
        input_file_entry.place(relx=0.3, rely=0.1, relwidth=0.675, relheight=0.2)

        annot_file_title = tk.Label(file_frame, text='Annotation file path: ')
        annot_file_title.place(relx=0.02, rely=0.4, relheight=0.2, relwidth=0.25)

        annot_file_entry = tk.Entry(file_frame)
        annot_file_entry.place(relx=0.3, rely=0.4, relwidth=0.675, relheight=0.2)

        # Radio buttons to specify Interaction or Edge network
        var_interaction_or_edge = tk.IntVar()
        var_interaction_or_edge.set(1)
        node_radio_button = tk.Radiobutton(file_frame, text='Show Interaction Nodes',
                                           variable=var_interaction_or_edge, value=1)
        node_radio_button.place(relx=0.02, rely=0.7, relwidth=0.35, relheight=0.2)

        edge_radio_button = tk.Radiobutton(file_frame, text='Show edges', variable=var_interaction_or_edge, value=2)
        edge_radio_button.place(relx=0.4, rely=0.7, relwidth=0.2, relheight=0.2)

        load_button = tk.Button(file_frame, bg='#90afc5', text="Load files",
                                command=lambda: self.load_files(input_file_entry.get(), annot_file_entry.get(),
                                                                var_interaction_or_edge.get()))
        load_button.place(relx=0.725, rely=0.7, relheight=0.2, relwidth=0.25)

        # Frame to specify styles for the network
        view_frame = tk.Frame(root, bg='#336b87', bd=5)
        view_frame.place(relx=0.258, rely=0.31, relwidth=0.468, relheight=0.575, anchor='n')

        view_frame_title = tk.Label(view_frame, bg='#336b87', text='Specify how you would like to put styles!')
        view_frame_title.place(relx=0.04, rely=0.05, relheight=0.1, relwidth=0.95)

        node_colour_title = tk.Label(view_frame, bg='#90afc5', justify='left', text='Node colour by: ')
        node_colour_title.place(relx=0.04, rely=0.2, relheight=0.1, relwidth=0.45)

        node_color_list = ['Order', 'Type', 'Default']
        node_colour_variable = tk.StringVar(view_frame)
        node_colour_variable.set(node_color_list[2])
        node_colour_options = tk.OptionMenu(view_frame, node_colour_variable, *node_color_list,
                                            command=lambda x: self.node_colour(node_colour_variable.get(),
                                                                               input_file_entry.get(),
                                                                               annot_file_entry.get(),
                                                                               var_interaction_or_edge.get()))
        node_colour_options.place(relx=0.525, rely=0.2, relheight=0.1, relwidth=0.45)

        node_size_title = tk.Label(view_frame, bg='#90afc5', text='Node size by: ')
        node_size_title.place(relx=0.04, rely=0.35, relheight=0.1, relwidth=0.45)

        node_size_list = ['Order', 'Type', 'Default']
        node_size_variable = tk.StringVar(view_frame)
        node_size_variable.set(node_size_list[2])
        node_size_options = tk.OptionMenu(view_frame, node_size_variable, *node_size_list,
                                          command=lambda x: self.node_size(node_size_variable.get(),
                                                                           input_file_entry.get(),
                                                                           annot_file_entry.get(),
                                                                           var_interaction_or_edge.get()))
        node_size_options.place(relx=0.525, rely=0.35, relheight=0.1, relwidth=0.45)

        node_shape_title = tk.Label(view_frame, bg='#90afc5', text='Node shape by: ')
        node_shape_title.place(relx=0.04, rely=0.5, relheight=0.1, relwidth=0.45)

        node_shape_list = ['Order', 'Type', 'Default']
        node_shape_variable = tk.StringVar(view_frame)
        node_shape_variable.set(node_shape_list[2])
        node_shape_options = tk.OptionMenu(view_frame, node_shape_variable, *node_shape_list,
                                           command=lambda x: self.node_shape(node_shape_variable.get(),
                                                                             input_file_entry.get(),
                                                                             annot_file_entry.get(),
                                                                             var_interaction_or_edge.get()))
        node_shape_options.place(relx=0.525, rely=0.5, relheight=0.1, relwidth=0.45)

        edge_colour_title = tk.Label(view_frame, bg='#90afc5', text='Edge colour by: ')
        edge_colour_title.place(relx=0.04, rely=0.65, relheight=0.1, relwidth=0.45)

        edge_colour_list = ['Order', 'Default']
        edge_colour_variable = tk.StringVar(view_frame)
        edge_colour_variable.set(edge_colour_list[1])
        edge_colour_options = tk.OptionMenu(view_frame, edge_colour_variable, *edge_colour_list,
                                            command=lambda x: self.edge_colour(edge_colour_variable.get(),
                                                                               input_file_entry.get(),
                                                                               annot_file_entry.get(),
                                                                               var_interaction_or_edge.get()))
        edge_colour_options.place(relx=0.525, rely=0.65, relheight=0.1, relwidth=0.45)

        edge_thickness_title = tk.Label(view_frame, bg='#90afc5', text='Edge thickness by: ')
        edge_thickness_title.place(relx=0.04, rely=0.8, relheight=0.1, relwidth=0.45)

        edge_thickness_list = ['Order', 'Type', 'Default']
        edge_thickness_variable = tk.StringVar(view_frame)
        edge_thickness_variable.set(edge_thickness_list[2])
        edge_thickness_options = tk.OptionMenu(view_frame, edge_thickness_variable, *edge_thickness_list,
                                               command=lambda x: self.edge_thickness(edge_thickness_variable.get(),
                                                                                     input_file_entry.get(),
                                                                                     annot_file_entry.get(),
                                                                                     var_interaction_or_edge.get()))
        edge_thickness_options.place(relx=0.525, rely=0.8, relheight=0.1, relwidth=0.45)

        # Frame to filter data
        filter_frame = tk.Frame(root, bg='#336b87', bd=5)
        filter_frame.place(relx=0.74, rely=0.31, relwidth=0.468, relheight=0.575, anchor='n')

        filter_frame_title = tk.Label(filter_frame, bg='#336b87', text='Specify how you would like to filter!')
        filter_frame_title.place(relx=0.04, rely=0.05, relheight=0.1, relwidth=0.95)

        filter_entry = tk.Entry(filter_frame, font=24)
        filter_entry.place(relx=0.05, rely=0.2, relwidth=0.925, relheight=0.3)

        hide_button = tk.Button(filter_frame, bg='#90afc5', text="Hide",
                                command=lambda: self.hide(input_file_entry.get(), annot_file_entry.get(),
                                                          var_interaction_or_edge.get()))
        hide_button.place(relx=0.04, rely=0.55, relheight=0.1, relwidth=0.45)

        show_button = tk.Button(filter_frame, bg='#90afc5', text="Show",
                                command=lambda: self.show(input_file_entry.get(), annot_file_entry.get(),
                                                          var_interaction_or_edge.get()))
        show_button.place(relx=0.525, rely=0.55, relheight=0.1, relwidth=0.45)

        hl_button = tk.Button(filter_frame, bg='#90afc5', text="Highlight",
                              command=lambda: self.highlight(input_file_entry.get(), annot_file_entry.get(),
                                                             var_interaction_or_edge.get()))
        hl_button.place(relx=0.04, rely=0.7, relheight=0.1, relwidth=0.45)

        gray_button = tk.Button(filter_frame, bg='#90afc5', text="Gray out",
                                command=lambda: self.grayout(input_file_entry.get(),
                                                             annot_file_entry.get(), var_interaction_or_edge.get()))
        gray_button.place(relx=0.525, rely=0.7, relheight=0.1, relwidth=0.45)

        reset_button = tk.Button(filter_frame, bg='#90afc5', text="Reset",
                                 command=lambda: self.reset(input_file_entry.get(),
                                                            annot_file_entry.get(), var_interaction_or_edge.get()))
        reset_button.place(relx=0.04, rely=0.85, relheight=0.1, relwidth=0.45)

        help_button = tk.Button(filter_frame, bg='#90afc5', text="Help",
                                command=lambda: self.help())
        help_button.place(relx=0.525, rely=0.85, relheight=0.1, relwidth=0.45)

        submit_button = tk.Button(root, text="Submit", command=lambda root=root: self.quit(root))
        submit_button.place(relx=0.3, rely=0.9, relheight=0.075, relwidth=0.35)

        root.mainloop()

        # DataFrame to send to Cytoscape_Integration class
        form_details_df = pd.DataFrame([[input_file_entry.get(), annot_file_entry.get(), node_colour_variable.get(),
                                         node_size_variable.get(), node_shape_variable.get()
                                            , edge_colour_variable.get(), edge_thickness_variable.get(),
                                         filter_entry.get()
                                            , self.hide_bool, self.show_bool, self.highlight_bool, self.gray_bool,
                                         self.reset_bool]]
                                       , columns=['input_file', 'annotation_file', 'node_colour', 'node_size',
                                                  'node_shape'
                , 'edge_colour', 'edge_thickness', 'query', 'hide', 'show', 'highlight', 'gray', 'reset'])

        form_details = [form_details_df, True, var_interaction_or_edge.get()]

        return form_details

    # Helper methods to update Cytoscape upon the press of buttons in the GUI
    def hide(self, input_file, annotation_file, interaction_or_edge):
        self.hide_bool = True
        form_details_df = pd.DataFrame([[input_file, annotation_file, self.hide_bool]],
                                       columns=['input_file', 'annotation_file', 'hide'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def show(self, input_file, annotation_file, interaction_or_edge):
        self.show_bool = True
        form_details_df = pd.DataFrame([[input_file, annotation_file, self.show_bool]],
                                       columns=['input_file', 'annotation_file', 'show'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def highlight(self, input_file, annotation_file, interaction_or_edge):
        self.highlight_bool = True
        form_details_df = pd.DataFrame([[input_file, annotation_file, self.highlight_bool]],
                                       columns=['input_file', 'annotation_file', 'highlight'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def grayout(self, input_file, annotation_file, interaction_or_edge):
        self.gray_bool = True
        form_details_df = pd.DataFrame([[input_file, annotation_file, self.gray_bool]],
                                       columns=['input_file', 'annotation_file', 'gray'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def reset(self, input_file, annotation_file, interaction_or_edge):
        self.reset_bool = True
        form_details_df = pd.DataFrame([[input_file, annotation_file, self.reset_bool]],
                                       columns=['input_file', 'annotation_file', 'reset'])
        self.controller.perform_core_functionality(form_details_df, False, interaction_or_edge)

    def help(self):
        tkMessageBox.showinfo('Info', 'This is how the query should be done: \nAaaaaaseemmennnyaaaa hfkewrhgethg4ui5!!')
        pass

    def quit(self, root):
        root.quit()

    def node_colour(self, node_colour, input_file, annotation_file, interaction_or_edge):
        form_details_df = pd.DataFrame([[input_file, annotation_file, node_colour]],
                                       columns=['input_file', 'annotation_file', 'node_colour'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def node_size(self, node_size, input_file, annotation_file, interaction_or_edge):
        form_details_df = pd.DataFrame([[input_file, annotation_file, node_size]],
                                       columns=['input_file', 'annotation_file', 'node_size'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def node_shape(self, node_shape, input_file, annotation_file, interaction_or_edge):
        form_details_df = pd.DataFrame([[input_file, annotation_file, node_shape]],
                                       columns=['input_file', 'annotation_file', 'node_shape'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def edge_colour(self, edge_colour, input_file, annotation_file, interaction_or_edge):
        form_details_df = pd.DataFrame([[input_file, annotation_file, edge_colour]],
                                       columns=['input_file', 'annotation_file', 'edge_colour'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

    def edge_thickness(self, edge_thickness, input_file, annotation_file, interaction_or_edge):
        form_details_df = pd.DataFrame([[input_file, annotation_file, edge_thickness]],
                                       columns=['input_file', 'annotation_file', 'edge_thickness'])
        self.controller.perform_core_functionality(form_details_df, True, interaction_or_edge)

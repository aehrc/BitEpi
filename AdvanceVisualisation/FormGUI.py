import Tkinter as tk


class FormGUI:

    def __init__(self):
        print("Inside GUI")

    # Method to create the form
    def form(self):
        root = tk.Tk()

        canvas = tk.Canvas(root, height=640, width=640)
        canvas.pack()

        main_title = tk.Label(root, text='Hi, please fill out the specifics you would like to see!')
        main_title.place(relx=0.1, rely=0.05, relheight=0.05, relwidth=0.75)

        input_file_frame = tk.Frame(root, bg='#a9a9a9', bd=5)
        input_file_frame.place(relx=0.5, rely=0.1, relwidth=0.95, relheight=0.1, anchor='n')

        file_title = tk.Label(input_file_frame, text='File path: ')
        file_title.place(relx=0.02, rely=0.35, relheight=0.35, relwidth=0.25)

        input_file_entry = tk.Entry(input_file_frame, font=24)
        input_file_entry.place(relx=0.3, rely=0.35, relwidth=0.675, relheight=0.35)

        view_frame = tk.Frame(root, bg='#808080', bd=5)
        view_frame.place(relx=0.258, rely=0.21, relwidth=0.468, relheight=0.65, anchor='n')

        view_frame_title = tk.Label(view_frame, bg='#808080', text='Specify how you would like to put styles!')
        view_frame_title.place(relx=0.04, rely=0.05, relheight=0.1, relwidth=0.95)

        node_colour_title = tk.Label(view_frame, bg='#696969', justify='left', text='Node colour by: ')
        node_colour_title.place(relx=0.04, rely=0.2, relheight=0.1, relwidth=0.45)

        node_color_list = ['None', 'Id', 'Order', 'Type']
        node_colour_variable = tk.StringVar(view_frame)
        node_colour_variable.set(node_color_list[0])
        node_colour_options = tk.OptionMenu(view_frame, node_colour_variable, *node_color_list)
        node_colour_options.place(relx=0.525, rely=0.2, relheight=0.1, relwidth=0.45)

        node_size_title = tk.Label(view_frame, bg='#696969', text='Node size by: ')
        node_size_title.place(relx=0.04, rely=0.35, relheight=0.1, relwidth=0.45)

        node_size_list = ['None', 'Alpha', 'Beta']
        node_size_variable = tk.StringVar(view_frame)
        node_size_variable.set(node_size_list[0])
        node_size_options = tk.OptionMenu(view_frame, node_size_variable, *node_size_list)
        node_size_options.place(relx=0.525, rely=0.35, relheight=0.1, relwidth=0.45)

        node_shape_title = tk.Label(view_frame, bg='#696969', text='Node shape by: ')
        node_shape_title.place(relx=0.04, rely=0.5, relheight=0.1, relwidth=0.45)

        node_shape_list = ['None', 'Alpha', 'Beta']
        node_shape_variable = tk.StringVar(view_frame)
        node_shape_variable.set(node_shape_list[0])
        node_shape_options = tk.OptionMenu(view_frame, node_shape_variable, *node_shape_list)
        node_shape_options.place(relx=0.525, rely=0.5, relheight=0.1, relwidth=0.45)

        edge_colour_title = tk.Label(view_frame, bg='#696969', text='Edge colour by: ')
        edge_colour_title.place(relx=0.04, rely=0.65, relheight=0.1, relwidth=0.45)

        edge_colour_list = ['None', 'Order 1', 'Order 2', 'Order 3', 'Order 4']
        edge_colour_variable = tk.StringVar(view_frame)
        edge_colour_variable.set(edge_colour_list[0])
        edge_colour_options = tk.OptionMenu(view_frame, edge_colour_variable, *edge_colour_list)
        edge_colour_options.place(relx=0.525, rely=0.65, relheight=0.1, relwidth=0.45)

        edge_thickness_title = tk.Label(view_frame, bg='#696969', text='Edge thickness by: ')
        edge_thickness_title.place(relx=0.04, rely=0.8, relheight=0.1, relwidth=0.45)

        edge_thickness_list = ['None', 'Alpha', 'Beta']
        edge_thickness_variable = tk.StringVar(view_frame)
        edge_thickness_variable.set(edge_thickness_list[0])
        edge_thickness_options = tk.OptionMenu(view_frame, edge_thickness_variable, *edge_thickness_list)
        edge_thickness_options.place(relx=0.525, rely=0.8, relheight=0.1, relwidth=0.45)

        filter_frame = tk.Frame(root, bg='#696969', bd=5)
        filter_frame.place(relx=0.74, rely=0.21, relwidth=0.468, relheight=0.65, anchor='n')

        button = tk.Button(root, text="Submit", font=24)
        button.place(relx=0.35, rely=0.9, relheight=0.075, relwidth=0.35)

        root.mainloop()
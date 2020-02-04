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
        input_file_entry.place(relx=0.3, rely=0.35, relwidth=0.7, relheight=0.35)

        view_frame = tk.Frame(root, bg='#808080', bd=5)
        view_frame.place(relx=0.5, rely=0.21, relwidth=0.95, relheight=0.3, anchor='n')

        filter_frame = tk.Frame(root, bg='#696969', bd=5)
        filter_frame.place(relx=0.5, rely=0.52, relwidth=0.95, relheight=0.35, anchor='n')

        button = tk.Button(root, text="Submit", font=24)
        button.place(relx=0.35, rely=0.9, relheight=0.075, relwidth=0.35)

        root.mainloop()

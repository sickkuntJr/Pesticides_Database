import os
import json
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw

class PesticideApp:
    def __init__(self, root, data):
        self.root = root
        self.root.title("Pesticide Database")
        self.data = data

        self.listbox = tk.Listbox(root)
        self.listbox.pack(side="left", fill="both", expand=True)

        self.details_frame = tk.Frame(root)
        self.details_frame.pack(side="right", fill="both", expand=True)

        self.name_label = tk.Label(self.details_frame, text="Name:")
        self.name_label.pack()

        self.cas_label = tk.Label(self.details_frame, text="CAS Number:")
        self.cas_label.pack()

        self.inchi_key_label = tk.Label(self.details_frame, text="InChIKey:")
        self.inchi_key_label.pack()

        self.structure_label = tk.Label(self.details_frame)
        self.structure_label.pack()

        for pesticide in data:
            self.listbox.insert("end", pesticide['name'])

        self.listbox.bind("<<ListboxSelect>>", self.on_select)

    def on_select(self, event):
        index = self.listbox.curselection()[0]
        pesticide = self.data[index]

        self.name_label.config(text=f"Name: {pesticide['name']}")
        self.cas_label.config(text=f"CAS Number: {pesticide['cas_number']}")
        self.inchi_key_label.config(text=f"InChIKey: {pesticide['inchi_key']}")

        image_folder = os.path.join(os.path.dirname(__file__), 'data', 'molecule_images')
        file_name = os.path.join(image_folder, f"{pesticide['name'].replace(' ', '_')}.png")
        if os.path.exists(file_name):
            structure_image = Image.open(file_name)
            structure_photo = ImageTk.PhotoImage(structure_image)
            self.structure_label.config(image=structure_photo, text='')
            self.structure_label.image = structure_photo
        else:
            self.structure_label.config(image='', text='Structure not available')

# Load data from JSON file
script_dir = os.path.dirname(os.path.realpath(__file__))
json_path = os.path.join(script_dir, 'data', 'pesticides_data.json')
with open(json_path, 'r') as f:
    data = json.load(f)

# Create the application
root = tk.Tk()
app = PesticideApp(root, data)
root.mainloop()

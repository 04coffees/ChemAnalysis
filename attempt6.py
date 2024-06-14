import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from PIL import Image, ImageTk
import io
import requests
import webview

class EsterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Ester Chemical Diagram Viewer")

        self.input_label = ttk.Label(root, text="Enter Esters (one per line):")
        self.input_label.grid(row=0, column=0, padx=10, pady=10, sticky="W")

        self.ester_text = scrolledtext.ScrolledText(root, height=20, width=40)
        self.ester_text.grid(row=1, column=0, padx=10, pady=10)

        self.canvas = tk.Canvas(root)
        self.canvas.grid(row=1, column=1, padx=10, pady=10, sticky="NSEW")

        self.diagram_frame = ttk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.diagram_frame, anchor="nw")

        self.scrollbar = ttk.Scrollbar(root, orient="vertical", command=self.canvas.yview)
        self.scrollbar.grid(row=1, column=2, sticky="NS")
        self.canvas.config(yscrollcommand=self.scrollbar.set)

        self.submit_button = ttk.Button(root, text="Submit", command=self.show_diagrams)
        self.submit_button.grid(row=2, column=0, columnspan=3, pady=10)

        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_columnconfigure(1, weight=1)

    def show_diagrams(self):
        self.diagram_frame.destroy()
        self.diagram_frame = ttk.Frame(self.canvas)
        self.canvas.create_window((0, 0), window=self.diagram_frame, anchor="nw")

        esters = self.ester_text.get("1.0", tk.END).strip().split('\n')
        
        for i, ester in enumerate(esters):
            ester_label = ttk.Label(self.diagram_frame, text=ester)
            ester_label.grid(row=i, column=0, padx=5, pady=5, sticky="W")

            smiles = self.get_smiles_from_name(ester)
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = self.mol_to_image(mol)
                    panel = tk.Label(self.diagram_frame, image=img)
                    panel.image = img
                    panel.grid(row=i, column=1, padx=5, pady=5)

                    category = self.categorize_shape(mol)
                    category_label = ttk.Label(self.diagram_frame, text=category)
                    category_label.grid(row=i, column=2, padx=5, pady=5)
                else:
                    error_label = tk.Label(self.diagram_frame, text="Invalid SMILES")
                    error_label.grid(row=i, column=1, padx=5, pady=5)
            else:
                error_label = tk.Label(self.diagram_frame, text="Could not resolve")
                error_label.grid(row=i, column=1, padx=5, pady=5)

        self.canvas.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

    def get_smiles_from_name(self, name):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
        else:
            return None

    def mol_to_image(self, mol):
        img = Draw.MolToImage(mol)
        bio = io.BytesIO()
        img.save(bio, format='PNG')
        bio.seek(0)
        img = Image.open(bio)
        return ImageTk.PhotoImage(img)

    def categorize_shape(self, mol):
        num_rings = Descriptors.RingCount(mol)
        num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        num_atoms = mol.GetNumAtoms()

        if num_rings == 0 and num_rotatable_bonds <= 3:
            category = "Linear and Flexible"
        elif num_rings == 0 and num_rotatable_bonds > 3:
            category = "Linear and Rigid"
        elif num_rings == 1:
            category = "Monocyclic"
        elif num_rings == 2:
            category = "Bicyclic"
        elif num_rings > 2:
            category = "Polycyclic"
        else:
            category = "Unknown"

        return category

def start_gui():
    root = tk.Tk()
    app = EsterApp(root)
    root.mainloop()

if __name__ == "__main__":
    start_gui()

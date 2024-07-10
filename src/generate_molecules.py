import os
import json
import requests
from rdkit import Chem
from rdkit.Chem import Draw

def get_smiles_from_inchikey(inchikey):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        print(f"Retrieved SMILES for InChIKey {inchikey}")
        return response.text.strip()
    else:
        print(f"Failed to retrieve SMILES for InChIKey {inchikey}")
        return None

def generate_structure_images(data, image_folder):
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)
        print(f"Created directory: {image_folder}")

    for pesticide in data:
        if pesticide['inchi_key']:
            smiles = get_smiles_from_inchikey(pesticide['inchi_key'])
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img_path = os.path.join(image_folder, f"{pesticide['name'].replace(' ', '_')}.png")
                    img = Draw.MolToImage(mol)
                    img.save(img_path)
                    print(f"Saved image for {pesticide['name']} at {img_path}")
                else:
                    print(f"Failed to generate molecule from SMILES for {pesticide['name']}")
            else:
                print(f"No SMILES found for {pesticide['name']}")
        else:
            print(f"No InChIKey found for {pesticide['name']}")

if __name__ == "__main__":
    # Full paths for Windows using raw strings
    json_path = r"C:\Users\Chris\Pesticides_Database\data\pesticides_data.json"
    image_folder = r"C:\Users\Chris\Pesticides_Database\data\molecule_images"
    
    with open(json_path, 'r') as f:
        data = json.load(f)
        print(f"Loaded data for {len(data)} pesticides")

    generate_structure_images(data, image_folder)

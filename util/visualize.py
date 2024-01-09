import base64
import io
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

def verify_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None

def render_smiles(smiles_list):
    mol_list = [
        Chem.MolFromSmiles(smile) for smile in smiles_list if verify_smiles(smile)
    ]

    if not mol_list:
        return None

    img = Draw.MolsToGridImage(mol_list, molsPerRow=4, subImgSize=(200, 100), useSVG=False)

    # Save to a byte array in JPEG format
    img_byte_arr = io.BytesIO()
    img.save(img_byte_arr, format='PNG', quality=100, optimize=True, progressive=True, compression='tiff_lzw')

    encoded_img = base64.b64encode(img_byte_arr.getvalue()).decode('ascii')

    return encoded_img

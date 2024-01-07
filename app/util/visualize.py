import io

import base64
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


def verify_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None


def render_smiles(smiles_list):
    mol_list = [Chem.MolFromSmiles(smile) for smile in smiles_list if verify_smiles(smile)]

    # Adjust subImgSize to fit aspect ratio
    img = Draw.MolsToGridImage(
        mol_list, molsPerRow=4, subImgSize=(200, 200), useSVG=False
    )

    img_byte_arr = io.BytesIO()
    img.save(img_byte_arr, format="PNG")
    encoded_img = base64.b64encode(img_byte_arr.getvalue()).decode('ascii')
    return encoded_img

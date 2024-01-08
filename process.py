import random

from loguru import logger

from util.midi3 import df_to_notes
from util.reaction import molecules_to_bond_energy_df
from util.visualize import render_smiles, verify_smiles


def chemsong(mols):
    # get bond energies from mols
    bond_df = molecules_to_bond_energy_df(mols)

    # visualize steps and display in window
    for step, smiles in mols.items():
        img = render_smiles(smiles)

    # play notes
    notes = df_to_notes(bond_df)
    logger.info(f"Notes mapped: {notes}")
    return notes


def process_reaction(entries):
    mols = {}
    for i, entry in enumerate(entries):
        step = i + 1
        smiles = entry
        if not verify_smiles(smiles):
            logger.error(f"Invalid SMILES notation at step {step}")
            return
        # Split the smiles string into a list and assign to the step
        mols[step] = smiles.split()
    notes = chemsong(mols)
    return notes


# Function to add a random step with 1-3 small chemicals
def random_step_value():
    random_chemicals = []  # Extend list as needed
    random_step = ""
    chemletters = ["C", "N", "O", "F", "S", "I", "c1ccccc1", "c1ccncc1"]
    for i in range(random.randint(1, 3)):
        random_step += "C" + random.choice(chemletters) + " "

    logger.info(f"Random step: {random_step}")

    return random_step

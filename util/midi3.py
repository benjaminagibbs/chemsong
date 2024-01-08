import time
from loguru import logger
import pandas as pd


from static.reference.bond_energies import bond_energies
from static.reference.scale_reference import *


# Data Processing Function
def normalize_data(bond_df: pd.DataFrame) -> pd.DataFrame:
    sorted_bond_energies = sorted(set(bond_energies.values()))

    # normalize energies for midi mapping
    bond_df["Energy Index"] = bond_df[
        "Bond Energies"
    ].apply(  # TODO: change to "make sure key ["energies"] matches the key in the dataframe"
        lambda energies: [sorted_bond_energies.index(e) for e in energies]
    )

    return bond_df


# MIDI Generation Function
def generate_midi(data: list[int]) -> list[dict]:
    midi_data = []
    for value in data:
        note = map_to_scale(value)
        midi_data.append({"note_on": note, "velocity": 64, "time": 0})
        midi_data.append({"note_off": note, "velocity": 64, "time": 500})
    return midi_data


# Map Numbers to C Minor Scale
def map_to_scale(note: int) -> int:
    # Mapping the number to a note in scale selected by dropdown

    scale = a_minor()

    # get dropdown value from where it's written:
    with open("static/reference/dropdown_value.txt", "r") as f:
        scale_selection = f.read()
        f.close()

    if scale_selection == "a_major":
        scale = a_major()
    elif scale_selection == "a_minor":
        scale = a_minor()
    elif scale_selection == "a_minor_h":
        scale = a_minor_h()
    elif scale_selection == "b_major":
        scale = b_major()
    elif scale_selection == "b_minor":
        scale = b_minor()

    return scale[note % len(scale)]



# Main synth function
def df_to_notes(df: pd.DataFrame):
    normalized_dataframe = normalize_data(df)

    for step in normalized_dataframe["Energy Index"]:
        time.sleep(0.25)
        midi_data = generate_midi(step)
    
    return midi_data
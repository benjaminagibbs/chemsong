# Chemsong

![Chemsong Window](img/chemsong_window.png)

## Overview

Map the energies of chemical bonds to musical scales.

## Features

- **Interactive Interface**: exists.
- **Chemical Reaction Visualization**: Utilizes RDKit for  rendering molecular structures based on user input.
- **Audio Mapping of Bond Energies**: Converts chemical bond energies into musical notes for an auditory representation of reactions.

## Upcoming Features

- **Object-Oriented Approach**: Future updates will introduce object-oriented programming (OOP) to manage each reaction step, incorporating methods for playing audio, verifying SMILES notation, and rendering images efficiently.
- **Advanced Audio Capabilities**: Enhancement of the audio system to provide more sophisticated and meaningful auditory representations of chemical reactions.

## Project Structure

- `app/`: Contains the main scripts to run Chemsong.
- `util/`: Includes essential subsystems and utilities for Chemsong.
- `RnD/`: Research and development scripts, including preliminary tests and experiments.

## Getting Started

To start using Chemsong, ensure you have Python (version 3.7 or later) installed, along with RDKit, Pandas, and other required libraries as listed in `requirements.txt`. Try firing up a virtual env and running:

`pip install -r requirements.txt`

## Example cmd output:
The bond energies are transcribed from source: Data from J. E. Huheey, E. A. Keiter, and R. L. Keiter, Inorganic Chemistry, 4th ed. (1993)

Bond energy values are output using loguru in the command line, example:

![Chemsong Command Line Output](img/chemsong_cmd_output.png)

## Contributing

Contributions are highly encouraged, whether it's through enhancing audio mappings, improving visualizations, or introducing new features.

## License

This project is open-sourced under the MIT License.

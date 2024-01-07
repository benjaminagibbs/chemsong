from flask import Flask, render_template, request, redirect, url_for

from loguru import logger

from process import *
from static.reference.scale_reference import scale_list


app = Flask(__name__)

sc = scale_list

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', scale_list=scale_list)

@app.route('/process', methods=['GET', 'POST'])
def process():
    if request.method == 'POST':
        musical_scale = request.form['musical_scale']
        smiles_entries = [request.form[key] for key in request.form if key.startswith('smiles_entry_')]

        logger.info(f"Musical Scale: {musical_scale}")
        with open("app/util/dropdown_value.txt", "w") as f:
            f.write(musical_scale)
            f.close()

        logger.info(f"SMILES Entries: {smiles_entries}")

        process_reaction(smiles_entries)  # Handle the processing as per your logic

    return redirect(url_for('index'))

@app.route('/generate_random_step')
def generate_random_step():
    breakpoint()
    # generate random step code goes here
    return render_template('index.html')

@app.route('/smiles-guide')
def smiles_guide():
    return render_template('smiles_guide.html')

if __name__ == '__main__':
    app.run(debug=True)

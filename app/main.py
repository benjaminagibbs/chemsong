from flask import Flask, render_template, request, redirect, url_for

from process import *
from static.reference.scale_reference import scale_list


app = Flask(__name__)

sc = scale_list

@app.route('/')
def index():
    scale_list = sc  # Your list of scales
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process():
    # Get the input from the form
    smiles = request.form.get('smiles')
    process_reaction(smiles)
    return redirect(url_for('index'))

@app.router('/add_step')
def add_step():
    # add step code goes here
    return redirect(url_for('index'))

@app.router('/generate_random_step')
def add_step():
    # generate random step code goes here
    return redirect(url_for('index'))

@app.router('/add_step')
def play_music():
    # play music code goes here
    return redirect(url_for('index'))


@app.route('/reset')
def reset():
    return redirect(url_for('index')) #TODO: check that this works after index is edited.

@app.route('/smiles-guide')
def smiles_guide():
    return render_template('smiles_guide.html')

if __name__ == '__main__':
    app.run(debug=True)

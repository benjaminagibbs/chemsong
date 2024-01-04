from flask import Flask, render_template, request, redirect, url_for
import markdown2

from process import *
from util.scale_reference import scale_list


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

@app.route('/reset')
def reset():
    return redirect(url_for('index'))

@app.route('/smiles-guide')
def smiles_guide():
    return render_template('smiles_guide.html')

if __name__ == '__main__':
    app.run(debug=True)

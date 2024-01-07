from flask import Flask, render_template, request, redirect, url_for, jsonify

from loguru import logger
import io

from process import *
from static.reference.scale_reference import scale_list


app = Flask(__name__)
log_stream = io.StringIO()

# Configure Loguru to write to log_stream
logger.add(log_stream, format="{time} {level} {message}")

@app.route('/get-logs')
def get_logs():
    # Return the logs as JSON
    return jsonify(logs=log_stream.getvalue())


@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html', scale_list=scale_list)

@app.route('/process', methods=['POST'])
def process():
    if request.method == 'POST':
        musical_scale = request.form['musical_scale']
        smiles_entries = [request.form[key] for key in request.form if key.startswith('smiles_entry_')]

        logger.info(f"Musical Scale: {musical_scale}")
        logger.info(f"SMILES Entries: {smiles_entries}")

        # Process the reactions (assuming this function works correctly)
        process_reaction(smiles_entries)

        # Generate images
        images = [render_smiles([smile]) for smile in smiles_entries]
        return jsonify({'result': 'success', 'data': smiles_entries, 'images': images})
    return jsonify({'result': 'error', 'message': 'Invalid request method'})


@app.route('/generate_random_step')
def generate_random_step():
    random_step = random_step_value()
    return jsonify(random_step)

@app.route('/smiles-guide')
def smiles_guide():
    return render_template('smiles_guide.html')

if __name__ == '__main__':
    app.run(debug=True)

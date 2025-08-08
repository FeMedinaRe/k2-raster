
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import os
import subprocess
import json

app = Flask(__name__)
CORS(app)

UPLOAD_FOLDER = os.path.abspath(os.path.join(os.path.dirname(__file__), '../build/bin'))
ALLOWED_EXTENSIONS = {'bin'}

os.makedirs(UPLOAD_FOLDER, exist_ok=True)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


# Endpoint para subir archivos .bin
@app.route('/api/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'success': False, 'error': 'No file part'}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({'success': False, 'error': 'No selected file'}), 400
    if file and allowed_file(file.filename):
        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)
        return jsonify({'success': True, 'filename': file.filename})
    return jsonify({'success': False, 'error': 'Invalid file type'}), 400

# Endpoint para codificar archivos
@app.route('/api/encode', methods=['POST'])
def encode_file():
    data = request.json
    filename = data.get('filename')
    if not filename or not allowed_file(filename):
        return jsonify({'success': False, 'error': 'Invalid filename'}), 400
    input_path = os.path.join(UPLOAD_FOLDER, filename)
    output_path = os.path.join(UPLOAD_FOLDER, 'raster.k2r')
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    encode_exe = os.path.join(project_root, 'build', 'bin', 'encode_k2r')
    get_values_exe = os.path.join(project_root, 'build', 'bin', 'get_values_window_k2r')
    query_txt = os.path.join(project_root, 'build', 'bin', 'query.txt')
    results_bin_backend = os.path.join(os.path.dirname(__file__), 'results.bin')
    results_bin_build = os.path.join(project_root, 'build', 'bin', 'results.bin')
    frontend_dir = os.path.abspath(os.path.join(project_root, 'frontend'))
    historial_path = os.path.join(os.path.dirname(__file__), 'historial_codificaciones.json')
    if not os.path.isfile(encode_exe):
        return jsonify({'success': False, 'error': 'encode_k2r not found'}), 500
    if not os.path.isfile(get_values_exe):
        return jsonify({'success': False, 'error': 'get_values_window_k2r not found'}), 500
    if not os.path.isfile(input_path):
        return jsonify({'success': False, 'error': 'Input file not found'}), 400
    if not os.path.isfile(query_txt):
        return jsonify({'success': False, 'error': 'query.txt not found'}), 400
    # 1. Codificar
    cmd_encode = [encode_exe, input_path, '640', '640', output_path, '-c', '-t', '10', '-k', '2']
    try:
        result_encode = subprocess.run(cmd_encode, capture_output=True, text=True, timeout=60)
        if result_encode.returncode != 0:
            return jsonify({'success': False, 'error': result_encode.stderr or 'Error in encode_k2r'})
        # 2. Ejecutar get_values_window_k2r
        cmd_get_values = [get_values_exe, output_path, query_txt, '-c']
        result_get_values = subprocess.run(cmd_get_values, capture_output=True, text=True, timeout=60)
        if result_get_values.returncode != 0:
            return jsonify({'success': False, 'error': result_get_values.stderr or 'Error in get_values_window_k2r'})
        # 3. Mover results.bin a frontend/
        import shutil
        if os.path.isfile(results_bin_build):
            src_results = results_bin_build
        elif os.path.isfile(results_bin_backend):
            src_results = results_bin_backend
        else:
            return jsonify({'success': False, 'error': 'results.bin not found after get_values_window_k2r'}), 500
        dest_results = os.path.join(frontend_dir, 'raster.bin')
        shutil.move(src_results, dest_results)

        # 4. Guardar historial de build time
        # Buscar la línea con 'k2-raster build time' en el output
        build_time_line = ''
        for line in (result_encode.stdout + '\n' + result_get_values.stdout).split('\n'):
            if 'k2-raster build time' in line:
                build_time_line = line.strip()
                break
        # Leer historial existente o crear uno nuevo
        try:
            if os.path.isfile(historial_path):
                with open(historial_path, 'r', encoding='utf-8') as f:
                    historial = json.load(f)
            else:
                historial = []
        except Exception:
            historial = []
        # Agregar nueva entrada
        historial.append({
            'fecha': __import__('datetime').datetime.now().isoformat(),
            'archivo': filename,
            'build_time': build_time_line
        })
        # Guardar historial actualizado
        with open(historial_path, 'w', encoding='utf-8') as f:
            json.dump(historial, f, ensure_ascii=False, indent=2)

        return jsonify({
            'success': True,
            'output': result_encode.stdout + '\n' + result_get_values.stdout,
            'output_file': 'raster.k2r',
            'results_bin': 'frontend/raster.bin',
            'build_time': build_time_line
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

# Endpoint para servir el historial como texto plano (JSON)
@app.route('/api/historial_raw', methods=['GET'])
def get_historial_raw():
    historial_path = os.path.join(os.path.dirname(__file__), 'historial_codificaciones.json')
    if os.path.isfile(historial_path):
        with open(historial_path, 'r', encoding='utf-8') as f:
            return f.read(), 200, {'Content-Type': 'application/json'}
    else:
        return '[]', 200, {'Content-Type': 'application/json'}

# Endpoint para consultar el historial de codificaciones
@app.route('/api/historial', methods=['GET'])
def get_historial():
    historial_path = os.path.join(os.path.dirname(__file__), 'historial_codificaciones.json')
    if os.path.isfile(historial_path):
        try:
            with open(historial_path, 'r', encoding='utf-8') as f:
                historial = json.load(f)
        except Exception:
            historial = []
    else:
        historial = []
    return jsonify({'success': True, 'historial': historial})
    if not os.path.isfile(input_path):
        return jsonify({'success': False, 'error': 'Input file not found'}), 400
    if not os.path.isfile(query_txt):
        return jsonify({'success': False, 'error': 'query.txt not found'}), 400
    # 1. Codificar
    cmd_encode = [encode_exe, input_path, '640', '640', output_path, '-c', '-t', '10', '-k', '2']
    try:
        result_encode = subprocess.run(cmd_encode, capture_output=True, text=True, timeout=60)
        if result_encode.returncode != 0:
            return jsonify({'success': False, 'error': result_encode.stderr or 'Error in encode_k2r'})
        # 2. Ejecutar get_values_window_k2r
        cmd_get_values = [get_values_exe, output_path, query_txt, '-c']
        result_get_values = subprocess.run(cmd_get_values, capture_output=True, text=True, timeout=60)
        if result_get_values.returncode != 0:
            return jsonify({'success': False, 'error': result_get_values.stderr or 'Error in get_values_window_k2r'})
        # 3. Mover results.bin a frontend/
        import shutil
        if os.path.isfile(results_bin_build):
            src_results = results_bin_build
        elif os.path.isfile(results_bin_backend):
            src_results = results_bin_backend
        else:
            return jsonify({'success': False, 'error': 'results.bin not found after get_values_window_k2r'}), 500
        dest_results = os.path.join(frontend_dir, 'raster.bin')
        shutil.move(src_results, dest_results)

        # 4. Guardar historial de build time
        # Buscar la línea con 'k2-raster build time' en el output
        build_time_line = ''
        for line in (result_encode.stdout + '\n' + result_get_values.stdout).split('\n'):
            if 'k2-raster build time' in line:
                build_time_line = line.strip()
                break
        # Leer historial existente o crear uno nuevo
        try:
            if os.path.isfile(historial_path):
                with open(historial_path, 'r', encoding='utf-8') as f:
                    historial = json.load(f)
            else:
                historial = []
        except Exception:
            historial = []
        # Agregar nueva entrada
        historial.append({
            'fecha': __import__('datetime').datetime.now().isoformat(),
            'archivo': filename,
            'build_time': build_time_line
        })
        # Guardar historial actualizado
        with open(historial_path, 'w', encoding='utf-8') as f:
            json.dump(historial, f, ensure_ascii=False, indent=2)

        return jsonify({
            'success': True,
            'output': result_encode.stdout + '\n' + result_get_values.stdout,
            'output_file': 'raster.k2r',
            'results_bin': 'frontend/raster.bin',
            'build_time': build_time_line
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)

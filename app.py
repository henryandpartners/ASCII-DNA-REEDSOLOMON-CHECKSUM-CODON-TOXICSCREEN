from flask import Flask, render_template, request, jsonify
import zlib
from reedsolo import RSCodec, ReedSolomonError
import base64
from PIL import Image
import io

app = Flask(__name__)

# --- Constants from Gemini 9A.py ---
RS_NSIZE = 255
RS_KSIZE = 223
RS_NSYM = RS_NSIZE - RS_KSIZE

DNA_MAP = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T'
}

REV_DNA_MAP = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}

# --- ASCII Art Conversion ---
ASCII_CHARS = "@%#*+=-:. "

def resize_image(image, new_width=100):
    width, height = image.size
    ratio = height / width
    new_height = int(new_width * ratio)
    resized_image = image.resize((new_width, new_height))
    return resized_image

def grayify(image):
    return image.convert("L")

def pixels_to_ascii(image):
    pixels = image.getdata()
    characters = "".join([ASCII_CHARS[pixel * len(ASCII_CHARS) // 256] for pixel in pixels])
    return characters

def image_to_ascii(image_file, new_width=100):
    try:
        image = Image.open(image_file)
    except Exception as e:
        return f"Error opening image: {e}"

    image = resize_image(image, new_width)
    image = grayify(image)

    ascii_str = pixels_to_ascii(image)

    img_width = image.width
    ascii_str_len = len(ascii_str)
    ascii_img = ""
    for i in range(0, ascii_str_len, img_width):
        ascii_img += ascii_str[i:i+img_width] + "\n"

    return ascii_img

# --- Helper Functions from Gemini 9A.py (refactored for web app) ---

def compress_data(data: bytes) -> bytes:
    return zlib.compress(data, level=9)

def add_error_correction(data: bytes) -> tuple[bytes, int]:
    rs = RSCodec(RS_NSYM)
    padding_needed = (RS_KSIZE - (len(data) % RS_KSIZE)) % RS_KSIZE
    padded_data = data + b'\x00' * padding_needed
    encoded_blocks = []
    for i in range(0, len(padded_data), RS_KSIZE):
        block = padded_data[i : i + RS_KSIZE]
        encoded_block = rs.encode(block)
        encoded_blocks.append(encoded_block)
    return b''.join(encoded_blocks), padding_needed

def encode_to_dna(data: bytes, scramble_key: int = 123) -> str:
    dna_sequence_parts = []
    scrambled_data = bytearray(len(data))
    for i, byte_val in enumerate(data):
        scrambled_data[i] = byte_val ^ ((scramble_key + i) % 256)
    for byte_val in scrambled_data:
        d3 = (byte_val >> 6) & 0b11
        d2 = (byte_val >> 4) & 0b11
        d1 = (byte_val >> 2) & 0b11
        d0 = byte_val & 0b11
        dna_sequence_parts.append(DNA_MAP[d3])
        dna_sequence_parts.append(DNA_MAP[d2])
        dna_sequence_parts.append(DNA_MAP[d1])
        dna_sequence_parts.append(DNA_MAP[d0])
    return "".join(dna_sequence_parts)

def decode_from_dna(dna_sequence: str, scramble_key: int = 123) -> bytes:
    if len(dna_sequence) % 4 != 0:
        raise ValueError("DNA sequence length must be a multiple of 4.")
    base4_digits = [REV_DNA_MAP.get(base, 0) for base in dna_sequence]
    decoded_bytes = bytearray()
    for i in range(0, len(base4_digits), 4):
        byte_val = (base4_digits[i] << 6) | (base4_digits[i+1] << 4) | (base4_digits[i+2] << 2) | base4_digits[i+3]
        decoded_bytes.append(byte_val)
    descrambled_data = bytearray(len(decoded_bytes))
    for i, byte_val in enumerate(decoded_bytes):
        descrambled_data[i] = byte_val ^ ((scramble_key + i) % 256)
    return bytes(descrambled_data)

def remove_error_correction(data_with_ecc: bytes, padding_needed: int) -> bytes:
    rs = RSCodec(RS_NSYM)
    decoded_blocks = []
    for i in range(0, len(data_with_ecc), RS_NSIZE):
        block = data_with_ecc[i : i + RS_NSIZE]
        if len(block) < RS_NSIZE:
            block += b'\x00' * (RS_NSIZE - len(block))
        try:
            message, _, _ = rs.decode(block)
            decoded_blocks.append(message)
        except ReedSolomonError:
            decoded_blocks.append(b'\x00' * RS_KSIZE)
    combined_decoded_data = b''.join(decoded_blocks)
    if padding_needed > 0:
        return combined_decoded_data[:-padding_needed]
    return combined_decoded_data

def decompress_data(data: bytes) -> bytes:
    return zlib.decompress(data)

# --- Main Interface Functions from Gemini 9A.py (adapted for web app) ---

def encode_ascii_to_dna_app(ascii_data: str) -> tuple[str, int]:
    original_data = ascii_data.encode('utf-8')
    compressed_data = compress_data(original_data)
    rs_encoded_data, padding_needed = add_error_correction(compressed_data)
    dna_sequence = encode_to_dna(rs_encoded_data)
    return dna_sequence, padding_needed

def decode_dna_to_ascii_app(dna_data: str, padding_needed: int = 0) -> str:
    decoded_rs_data = decode_from_dna(dna_data)
    decompressed_data = remove_error_correction(decoded_rs_data, padding_needed)
    original_data = decompress_data(decompressed_data)
    try:
        return original_data.decode('utf-8')
    except UnicodeDecodeError:
        return original_data.hex()

# --- Flask Routes ---

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/encode', methods=['POST'])
def encode():
    data = request.json
    ascii_data = data.get('text', '')
    try:
        dna_sequence, padding_needed = encode_ascii_to_dna_app(ascii_data)
        return jsonify({
            'success': True,
            'dna': dna_sequence,
            'padding': padding_needed,
            'base_pairs': len(dna_sequence)
        })
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/decode', methods=['POST'])
def decode():
    data = request.json
    dna_data = data.get('dna', '')
    padding_needed = int(data.get('padding', 0))
    try:
        ascii_result = decode_dna_to_ascii_app(dna_data, padding_needed)
        return jsonify({'success': True, 'text': ascii_result})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

@app.route('/image-to-ascii', methods=['POST'])
def convert_image_to_ascii():
    if 'image' not in request.files:
        return jsonify({'success': False, 'error': 'No image file provided'})

    file = request.files['image']

    if file.filename == '':
        return jsonify({'success': False, 'error': 'No selected file'})

    try:
        ascii_art = image_to_ascii(file)
        return jsonify({'success': True, 'ascii': ascii_art})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080)

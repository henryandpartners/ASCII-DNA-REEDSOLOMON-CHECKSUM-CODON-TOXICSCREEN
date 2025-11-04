import zlib
from reedsolo import RSCodec, ReedSolomonError
import sys

# --- Constants from original code ---
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

# Conceptual 6-base DNA map (for illustration - needs proper definition based on how 6 bases map to bits)
# This is a placeholder and would need a robust design for actual implementation.
DNA_MAP_6_BASE = {
    0: 'A',
    1: 'C',
    2: 'G',
    3: 'T',
    4: 'X', # Placeholder for additional bases
    5: 'Y'  # Placeholder for additional bases
}

REV_DNA_MAP_6_BASE = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
    'X': 4,
    'Y': 5
}

# --- Original Helper Functions ---

def compress_data(data: bytes) -> bytes:
    print(f"[Compression] Original size: {len(data)} bytes")
    compressed_data = zlib.compress(data, level=9)
    print(f"[Compression] Compressed size: {len(compressed_data)} bytes")
    return compressed_data

def add_error_correction(data: bytes) -> tuple[bytes, int]:
    rs = RSCodec(RS_NSYM)

    padding_needed = (RS_KSIZE - (len(data) % RS_KSIZE)) % RS_KSIZE
    padded_data = data + b'\x00' * padding_needed

    encoded_blocks = []
    for i in range(0, len(padded_data), RS_KSIZE):
        block = padded_data[i : i + RS_KSIZE]
        encoded_block = rs.encode(block)
        encoded_blocks.append(encoded_block)

    rs_corrected_data = b''.join(encoded_blocks)
    print(f"[Reed-Solomon] Encoded size: {len(rs_corrected_data)} bytes")
    return rs_corrected_data, padding_needed

def encode_to_dna(data: bytes, scramble_key: int = 123, dna_map=DNA_MAP) -> str:
    dna_sequence_parts = []

    scrambled_data = bytearray(len(data))
    for i, byte_val in enumerate(data):
        scrambled_data[i] = byte_val ^ ((scramble_key + i) % 256)

    for byte_val in scrambled_data:
        # Each byte (8 bits) is split into four 2-bit chunks
        d3 = (byte_val >> 6) & 0b11
        d2 = (byte_val >> 4) & 0b11
        d1 = (byte_val >> 2) & 0b11
        d0 = byte_val & 0b11

        dna_sequence_parts.append(dna_map[d3])
        dna_sequence_parts.append(dna_map[d2])
        dna_sequence_parts.append(dna_map[d1])
        dna_sequence_parts.append(dna_map[d0])

    dna_sequence = "".join(dna_sequence_parts)
    print(f"[DNA Encoding] Generated DNA sequence with {len(dna_sequence)} base pairs")
    return dna_sequence

def decode_from_dna(dna_sequence: str, scramble_key: int = 123, rev_dna_map=REV_DNA_MAP) -> bytes:
    if len(dna_sequence) % 4 != 0:
        raise ValueError("DNA sequence length must be a multiple of 4 for proper byte decoding.")

    base4_digits = []
    for base in dna_sequence:
        if base not in rev_dna_map:
            print(f"Warning: Invalid DNA base '{base}' found during decoding. Using default 0.")
            base4_digits.append(0)
        else:
            base4_digits.append(rev_dna_map[base])

    decoded_bytes = bytearray()
    for i in range(0, len(base4_digits), 4):
        d3 = base4_digits[i]
        d2 = base4_digits[i+1]
        d1 = base4_digits[i+2]
        d0 = base4_digits[i+3]

        byte_val = (d3 << 6) | (d2 << 4) | (d1 << 2) | d0
        decoded_bytes.append(byte_val)

    descrambled_data = bytearray(len(decoded_bytes))
    for i, byte_val in enumerate(decoded_bytes):
        descrambled_data[i] = byte_val ^ ((scramble_key + i) % 256)

    return bytes(descrambled_data)

def remove_error_correction(data_with_ecc: bytes, padding_needed: int) -> bytes:
    rs = RSCodec(RS_NSYM)

    decoded_blocks = []
    errors_fixed_total = 0

    for i in range(0, len(data_with_ecc), RS_NSIZE):
        block = data_with_ecc[i : i + RS_NSIZE]
        if len(block) != RS_NSIZE:
            print(f"Warning: Incomplete RS block (length {len(block)}, expected {RS_NSIZE})")
            if len(block) < RS_NSIZE:
                 block += b'\x00' * (RS_NSIZE - len(block))

        try:
            message, _, errata_pos = rs.decode(block)
            decoded_blocks.append(message)
            errors_fixed_total += len(errata_pos)
        except ReedSolomonError as e:
            print(f"Error: Reed-Solomon could not correct block starting at byte {i}")
            decoded_blocks.append(b'\x00' * RS_KSIZE)

    combined_decoded_data = b''.join(decoded_blocks)

    if padding_needed > 0:
        original_data = combined_decoded_data[:-padding_needed]
    else:
        original_data = combined_decoded_data

    print(f"[Reed-Solomon] Errors fixed: {errors_fixed_total}")
    return original_data

def decompress_data(data: bytes) -> bytes:
    print(f"[Decompression] Compressed size: {len(data)} bytes")
    decompressed_data = zlib.decompress(data)
    print(f"[Decompression] Decompressed size: {len(decompressed_data)} bytes")
    return decompressed_data

# --- New Interface Functions ---

def encode_ascii_to_dna(ascii_data: str, use_6_base_codon_raw: bool = False, use_6_base_codon_no_stop: bool = False) -> str:
    """
    Encodes ASCII data into a DNA sequence.
    """
    original_data = ascii_data.encode('utf-8')
    print(f"Encoding ASCII data: '{ascii_data[:50]}...' (first 50 chars)" if len(ascii_data) > 50 else f"Encoding ASCII data: '{ascii_data}'")

    compressed_data = compress_data(original_data)
    rs_encoded_data, padding_needed = add_error_correction(compressed_data)

    if use_6_base_codon_raw or use_6_base_codon_no_stop:
        raise NotImplementedError("6-base codon encoding is not implemented in this version. "
                                  "It requires a fundamental change to the bit-to-base mapping.")
    else:
        dna_sequence = encode_to_dna(rs_encoded_data, scramble_key=123, dna_map=DNA_MAP)

    return dna_sequence, padding_needed

def decode_dna_to_ascii(dna_data: str, use_6_base_codon_raw: bool = False, use_6_base_codon_no_stop: bool = False, padding_needed: int = 0) -> str:
    """
    Decodes a DNA sequence back into ASCII data.
    """
    print(f"Decoding DNA sequence: '{dna_data[:50]}...' (first 50 chars)" if len(dna_data) > 50 else f"Decoding DNA sequence: '{dna_data}'")

    if use_6_base_codon_raw or use_6_base_codon_no_stop:
        raise NotImplementedError("6-base codon decoding is not implemented in this version. "
                                  "It requires a fundamental change to the base-to-bit mapping.")
    else:
        decoded_rs_data = decode_from_dna(dna_data, scramble_key=123, rev_dna_map=REV_DNA_MAP)

    decompressed_data = remove_error_correction(decoded_rs_data, padding_needed)
    original_data = decompress_data(decompressed_data)

    try:
        return original_data.decode('utf-8')
    except UnicodeDecodeError:
        return original_data.hex() # Fallback to hex if not valid ASCII/UTF-8

def main_interface():
    while True:
        print("\n--- DNA Data Storage System ---")
        print("1. Encode ASCII to DNA")
        print("2. Decode DNA to ASCII")
        print("3. Exit")
        choice = input("Enter your choice: ")

        if choice == '1':
            print("Enter ASCII data (type 'END' on a new line to finish input):")
            lines = []
            while True:
                line = sys.stdin.readline().strip()
                if line == "END":
                    break
                lines.append(line)
            ascii_input = "\n".join(lines)
            
            try:
                dna_output, padding_needed_for_decode = encode_ascii_to_dna(ascii_input)
                print(f"\n--- DNA Sequence Result (Padding needed for decoding: {padding_needed_for_decode}) ---")
                print(dna_output)
            except NotImplementedError as e:
                print(f"Error: {e}")
            except Exception as e:
                print(f"An unexpected error occurred during encoding: {e}")

        elif choice == '2':
            print("Enter DNA sequence (type 'END' on a new line to finish input):")
            lines = []
            while True:
                line = sys.stdin.readline().strip()
                if line == "END":
                    break
                lines.append(line)
            dna_input = "".join(lines)
            
            try:
                padding_str = input("Enter padding needed during original encoding (0 if unknown/none): ")
                padding_needed = int(padding_str) if padding_str.isdigit() else 0
            except ValueError:
                print("Invalid padding value. Using 0.")
                padding_needed = 0

            try:
                decoded_result = decode_dna_to_ascii(dna_input, padding_needed=padding_needed)
                print("\n--- Decoded ASCII Result ---")
                print(decoded_result)
            except NotImplementedError as e:
                print(f"Error: {e}")
            except Exception as e:
                print(f"An unexpected error occurred during decoding: {e}")

        elif choice == '3':
            print("Exiting.")
            break
        else:
            print("Invalid choice. Please try again.")

if __name__ == "__main__":
    main_interface()

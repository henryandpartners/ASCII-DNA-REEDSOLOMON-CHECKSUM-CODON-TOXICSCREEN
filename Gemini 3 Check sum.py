import itertools
import reedsolomon as rs # This is a placeholder, you'd use a real library

class DNAScrambler:
    """
    A class to encrypt ASCII text (including multi-line) to DNA sequences and
    decrypt DNA sequences back to ASCII. It now uses a checksum-based system
    for integrity checking. It also provides functionality to screen DNA sequences
    for specific DNA codons, and for biological translation, protein motif screening,
    codon optimization, and basic toxic protein screening.
    """

    def __init__(self, data_symbols=16, parity_symbols=4): # Reed-Solomon parameters (still present for ECC, but main focus shifts to checksum for basic integrity)
        """
        Initializes the DNAScrambler with all possible DNA bases, generates
        the 64 unique 3-base codons, and sets up the genetic code.
        """
        # Initialize all required attributes
        self.bases = ['A', 'T', 'C', 'G']
        self.all_codons = self._generate_all_codons()
        # Map integers to codons and codons to integers for 64 codons
        self.int_to_codon_map = {i: codon for i, codon in enumerate(self.all_codons)}
        self.codon_to_int_map = {codon: i for i, codon in enumerate(self.all_codons)}
        # Example genetic code (should be filled with all codons)
        self.genetic_code = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", "CCT": "P", "CCC": "P",
            "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
            "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "TGT": "C", "TGC": "C", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
            "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", "TAA": "*", "TAG": "*", "TGA": "*"
        }
        self.start_codons = {"ATG"}
        self.stop_codons = {"TAA", "TAG", "TGA"}
        self.complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        # Example codon usage (E. coli, partial)
        self.codon_usage = {
            'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'],
            'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
            'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
            'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA']
        }
        # Toxic protein motifs
        self.toxic_protein_motifs = {
            "Lysine Rich Region": ["KKKKK", "RRRRR"],
            "Hydrophobic Patch": ["LLLL", "IIII", "VVVV"],
            "Repeated Proline": ["PPPP"],
        }
        self.data_symbols = data_symbols
        self.parity_symbols = parity_symbols
        self.total_symbols = self.data_symbols + self.parity_symbols

    def _generate_all_codons(self):
        """Helper method to generate all 64 unique 3-base DNA codons."""
        return [''.join(p) for p in itertools.product(self.bases, repeat=3)]

    # --- Run-Length Encoding (RLE) Methods ---
    def _rle_compress(self, text):
        """
        Compresses a string using simple Run-Length Encoding.
        Format: <count><char> (e.g., "AAAA" becomes "4A").
        Counts are simple integers. If count is 1, just the char.
        We need to escape numbers in the original text, to avoid confusion with counts.
        A simple escape is to prefix a backslash \ to numbers in the original text.
        e.g., "123ABC" -> "\1\2\3ABC"
        """
        if not text:
            return ""

        compressed_parts = []
        i = 0
        while i < len(text):
            char = text[i]
            count = 1
            j = i + 1
            while j < len(text) and text[j] == char:
                count += 1
                j += 1
            
            # Escape numbers in the original text if they appear
            # to prevent confusion with RLE counts.
            if char.isdigit():
                compressed_parts.append(f"\\{char}")
            elif count > 1:
                compressed_parts.append(f"{count}{char}")
            else:
                compressed_parts.append(char)
            i = j
        return "".join(compressed_parts)

    def _rle_decompress(self, compressed_text):
        """
        Decompresses a string encoded with simple Run-Length Encoding.
        Handles escaped numbers.
        """
        if not compressed_text:
            return ""

        decompressed_parts = []
        i = 0
        while i < len(compressed_text):
            if compressed_text[i] == '\\':
                # Handle escaped number
                if i + 1 < len(compressed_text) and compressed_text[i+1].isdigit():
                    decompressed_parts.append(compressed_text[i+1])
                    i += 2
                else:
                    # Malformed escape, or literal backslash
                    decompressed_parts.append('\\')
                    i += 1
            elif compressed_text[i].isdigit():
                # Read count
                count_str = ""
                while i < len(compressed_text) and compressed_text[i].isdigit():
                    count_str += compressed_text[i]
                    i += 1
                count = int(count_str)
                
                # The next character is the one to repeat
                if i < len(compressed_text):
                    char_to_repeat = compressed_text[i]
                    decompressed_parts.append(char_to_repeat * count)
                    i += 1
                else:
                    # Malformed RLE (count without a character)
                    raise ValueError(f"Malformed RLE: Count '{count_str}' at end of string.")
            else:
                # Just a single character
                decompressed_parts.append(compressed_text[i])
                i += 1
        return "".join(decompressed_parts)

    # --- Checksum-based Encrypt/Decrypt Methods ---
    def _calculate_checksum(self, codon_integers):
        """Calculates a simple sum-based checksum from a list of codon integer representations."""
        # Simple checksum: sum of integer values of codons, modulo 64 to get a codon index
        return sum(codon_integers) % len(self.all_codons)

    def encrypt_with_checksum(self, text):
        """
        Encrypts ASCII text into a DNA sequence using a direct character-to-codon mapping
        and appends a checksum codon for integrity verification.
        Each character maps to one codon. This limits us to 64 unique characters.
        If more characters are needed, a more complex mapping (e.g., 2 codons per char) is required.
        """
        encrypted_codons = []
        codon_integers = []
        
        # Check if all characters can be mapped to a single codon (0-63 ASCII)
        for char in text:
            ascii_val = ord(char)
            if not (0 <= ascii_val < len(self.all_codons)): # Only ASCII 0-63 can be directly mapped
                raise ValueError(f"Character '{char}' (ASCII value {ascii_val}) cannot be directly mapped to a single 3-base codon for checksum encryption. "
                                 "Only ASCII values 0-63 are supported with this simple mapping.")
            
            codon = self.int_to_codon_map[ascii_val]
            encrypted_codons.append(codon)
            codon_integers.append(ascii_val) # Use ASCII value as the integer representation

        # Calculate checksum
        checksum_int = self._calculate_checksum(codon_integers)
        checksum_codon = self.int_to_codon_map[checksum_int]
        
        encrypted_codons.append(checksum_codon) # Append checksum codon
        
        return "".join(encrypted_codons)

    def decrypt_with_checksum(self, dna_sequence):
        """
        Decrypts a DNA sequence that includes a checksum.
        Verifies integrity using the checksum.
        """
        decrypted_text_parts = []
        warnings = []
        
        if len(dna_sequence) % 3 != 0:
            warnings.append(f"Warning: Input DNA sequence length ({len(dna_sequence)} bases) is not a multiple of 3. "
                            "This may indicate corruption or incomplete data.")
            # Attempt to trim to a multiple of 3 for decoding
            dna_sequence = dna_sequence[:len(dna_sequence) - (len(dna_sequence) % 3)]

        if len(dna_sequence) < 3: # Need at least one codon (data) + one codon (checksum)
            raise ValueError("DNA sequence too short to contain data and checksum.")

        received_codons = []
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3].upper()
            if codon not in self.codon_to_int_map:
                warnings.append(f"Unknown codon '{codon}' found at DNA position {i}. Cannot decrypt this segment.")
                # For checksum validation, we might treat unknown codons as 0 or exclude them
                # For now, let's just append 0 to the integer list for checksum calculation
                received_codons.append(0) 
            else:
                received_codons.append(self.codon_to_int_map[codon])

        # The last codon is the checksum
        received_checksum_int = received_codons[-1]
        data_codon_integers = received_codons[:-1] # All but the last codon

        # Recalculate checksum from the data codons
        calculated_checksum_int = self._calculate_checksum(data_codon_integers)

        if calculated_checksum_int != received_checksum_int:
            warnings.append(f"Checksum mismatch! Expected: {self.int_to_codon_map[calculated_checksum_int]} (int: {calculated_checksum_int}), "
                            f"Received: {self.int_to_codon_map[received_checksum_int]} (int: {received_checksum_int}). Data integrity compromised.")
        else:
            warnings.append("Checksum validated: Data integrity confirmed.")

        # Decrypt the data codons
        for codon_int in data_codon_integers:
            if 0 <= codon_int < 128: # Ensure it's a valid ASCII range
                decrypted_text_parts.append(chr(codon_int))
            else:
                decrypted_text_parts.append('?') # Cannot map back to ASCII char

        return "".join(decrypted_text_parts), warnings

    # --- ECC Methods (kept as-is for now, separate from checksum) ---
    def encrypt_with_ecc(self, text):
        """
        Encrypts ASCII text into a DNA sequence with Reed-Solomon error correction.
        """
        byte_data = [ord(char) for char in text]
        
        remainder = len(byte_data) % self.data_symbols
        if remainder != 0:
            padding_needed = self.data_symbols - remainder
            byte_data.extend([0] * padding_needed)

        encoded_data_with_parity = []
        for i in range(0, len(byte_data), self.data_symbols):
            block = byte_data[i : i + self.data_symbols]
            
            # Placeholder for RS encoding:
            # You would replace this with actual Reed-Solomon encoding
            # For example: rs_block = self.rs_encoder.encode(bytearray(block))
            rs_block = list(block) + [0] * self.parity_symbols # Dummy parity
            
            encoded_data_with_parity.extend(rs_block)

        encrypted_dna_parts = []
        
        for val in encoded_data_with_parity:
            # Map byte value (0-255) to 3-base codons. 
            # This requires a more robust mapping than simple ASCII char to codon.
            # Here, we'll map each byte to two 3-base codons (6 bits per codon, so 12 bits total, 256 values fit).
            if len(self.all_codons) * len(self.all_codons) < 256:
                raise RuntimeError("Not enough unique codon pairs to map all 256 byte values for ECC.")

            codon1_idx = val // len(self.all_codons)
            codon2_idx = val % len(self.all_codons)
            encrypted_dna_parts.append(self.int_to_codon_map[codon1_idx])
            encrypted_dna_parts.append(self.int_to_codon_map[codon2_idx])
            
        return "".join(encrypted_dna_parts)

    def decrypt_with_ecc(self, dna_sequence):
        """
        Decrypts a DNA sequence (with Reed-Solomon error correction) back into an ASCII text string.
        """
        decrypted_integer_values_from_dna = []
        i = 0
        while i < len(dna_sequence) - 5: # Process in 6-base chunks (two codons per byte)
            codon1 = dna_sequence[i:i+3].upper()
            codon2 = dna_sequence[i+3:i+6].upper()
            
            if codon1 in self.codon_to_int_map and codon2 in self.codon_to_int_map:
                val = self.codon_to_int_map[codon1] * len(self.all_codons) + self.codon_to_int_map[codon2]
                decrypted_integer_values_from_dna.append(val)
            else:
                # If a segment doesn't map, it's an error. Reed-Solomon can correct.
                decrypted_integer_values_from_dna.append(0) # Placeholder for unknown/error
                
            i += 6

        corrected_byte_data = []
        errors_corrected_count = 0
        # Reed-Solomon decoding logic (placeholder)
        for i in range(0, len(decrypted_integer_values_from_dna), self.total_symbols):
            block = decrypted_integer_values_from_dna[i : i + self.total_symbols]
            
            if len(block) < self.total_symbols:
                break 

            # Placeholder for RS decoding:
            # rs_decoded_block, corrections = self.rs_decoder.decode(bytearray(block))
            rs_decoded_block = block[:self.data_symbols] # Dummy decode
            
            # if corrections > 0:
            #    errors_corrected_count += corrections
            
            corrected_byte_data.extend(rs_decoded_block)

        decrypted_text_parts = []
        for val in corrected_byte_data:
            if 0 <= val <= 127: # Only output valid ASCII characters
                decrypted_text_parts.append(chr(val))

        warnings = []
        if errors_corrected_count > 0:
             warnings.append(f"Reed-Solomon corrected {errors_corrected_count} errors during decryption.")
        
        return "".join(decrypted_text_parts), warnings

    # --- NEW: ECC with Checksum Methods ---
    def encrypt_with_ecc_and_checksum(self, text):
        """
        Encrypts ASCII text into a DNA sequence, applying Reed-Solomon ECC
        and then appending a checksum based on the ECC-encoded data.
        """
        # 1. ECC Encode the text
        # This part generates the integer sequence that would be converted to DNA.
        # We need to expose this intermediate step or refactor _encrypt_with_ecc
        # to return the encoded_data_with_parity (integers).
        byte_data = [ord(char) for char in text]
        
        remainder = len(byte_data) % self.data_symbols
        if remainder != 0:
            padding_needed = self.data_symbols - remainder
            byte_data.extend([0] * padding_needed)

        ecc_encoded_integers = []
        for i in range(0, len(byte_data), self.data_symbols):
            block = byte_data[i : i + self.data_symbols]
            # Placeholder for RS encoding:
            rs_block = list(block) + [0] * self.parity_symbols # Dummy parity
            ecc_encoded_integers.extend(rs_block)

        # 2. Convert ECC-encoded integers to DNA codons
        dna_parts = []
        for val in ecc_encoded_integers:
            if len(self.all_codons) * len(self.all_codons) < 256:
                raise RuntimeError("Not enough unique codon pairs to map all 256 byte values for ECC.")
            codon1_idx = val // len(self.all_codons)
            codon2_idx = val % len(self.all_codons)
            dna_parts.append(self.int_to_codon_map[codon1_idx])
            dna_parts.append(self.int_to_codon_map[codon2_idx])
        
        # 3. Calculate checksum on the ECC-encoded integer values (before DNA conversion)
        # This ensures the checksum verifies the integrity of the data that went into DNA.
        checksum_int = self._calculate_checksum(ecc_encoded_integers)
        checksum_codon = self.int_to_codon_map[checksum_int]

        # 4. Append checksum codon to the DNA sequence
        dna_parts.append(checksum_codon)

        return "".join(dna_parts)

    def decrypt_with_ecc_and_checksum(self, dna_sequence):
        """
        Decrypts a DNA sequence that was encrypted with ECC and a checksum.
        Performs ECC decoding and then checksum verification.
        """
        warnings = []
        
        if len(dna_sequence) < 9: # Needs at least 2 codons for 1 byte (ECC) + 1 codon for checksum
            raise ValueError("DNA sequence too short to contain ECC data and checksum.")
        
        # 1. Extract checksum codon (last 3 bases)
        received_checksum_codon = dna_sequence[-3:].upper()
        
        if received_checksum_codon not in self.codon_to_int_map:
            warnings.append(f"Warning: Unknown checksum codon '{received_checksum_codon}'. Checksum cannot be validated.")
            received_checksum_int = -1 # Indicate unknown
        else:
            received_checksum_int = self.codon_to_int_map[received_checksum_codon]

        # 2. Process the rest of the DNA sequence for ECC decryption
        ecc_dna_segment = dna_sequence[:-3] # DNA without the checksum
        
        decrypted_integer_values_from_ecc_dna = []
        i = 0
        while i < len(ecc_dna_segment) - 5: # Process in 6-base chunks (two codons per byte)
            codon1 = ecc_dna_segment[i:i+3].upper()
            codon2 = ecc_dna_segment[i+3:i+6].upper()
            
            if codon1 in self.codon_to_int_map and codon2 in self.codon_to_int_map:
                val = self.codon_to_int_map[codon1] * len(self.all_codons) + self.codon_to_int_map[codon2]
                decrypted_integer_values_from_ecc_dna.append(val)
            else:
                warnings.append(f"Warning: Unknown codon pair '{codon1}{codon2}' at DNA position {i} (ECC segment).")
                decrypted_integer_values_from_ecc_dna.append(0) # Placeholder for unknown/error
                
            i += 6

        # 3. Perform Reed-Solomon decoding on the integer values
        corrected_byte_data = []
        ecc_errors_corrected_count = 0
        # Reed-Solomon decoding logic (placeholder)
        for i in range(0, len(decrypted_integer_values_from_ecc_dna), self.total_symbols):
            block = decrypted_integer_values_from_ecc_dna[i : i + self.total_symbols]
            
            if len(block) < self.total_symbols:
                break 

            # Placeholder for RS decoding:
            # rs_decoded_block, corrections = self.rs_decoder.decode(bytearray(block))
            rs_decoded_block = block[:self.data_symbols] # Dummy decode
            
            # if corrections > 0:
            #    ecc_errors_corrected_count += corrections
            
            corrected_byte_data.extend(rs_decoded_block)
        
        if ecc_errors_corrected_count > 0:
             warnings.append(f"Reed-Solomon corrected {ecc_errors_corrected_count} errors during decryption.")

        # 4. Recalculate checksum from the *decoded ECC data* (corrected_byte_data)
        # Note: This checksum is calculated on the *byte values* that ECC produced,
        # not directly on the DNA codons. This verifies the integrity of the data stream itself.
        calculated_checksum_int = self._calculate_checksum(corrected_byte_data)

        # 5. Compare checksums
        if received_checksum_int != -1: # Only compare if received checksum was valid
            if calculated_checksum_int != received_checksum_int:
                warnings.append(f"Checksum mismatch! Expected: {self.int_to_codon_map[calculated_checksum_int]} (int: {calculated_checksum_int}), "
                                f"Received: {self.int_to_codon_map[received_checksum_int]} (int: {received_checksum_int}). Data integrity compromised AFTER ECC.")
            else:
                warnings.append("Checksum validated: Data integrity confirmed AFTER ECC.")
        
        # 6. Convert corrected byte data back to ASCII text
        decrypted_text_parts = []
        for val in corrected_byte_data:
            if 0 <= val <= 127:
                decrypted_text_parts.append(chr(val))

        return "".join(decrypted_text_parts), warnings

    def screen_codons_in_dna(self, dna_sequence, codon_type_map):
        """Screens a given DNA sequence for occurrences of specified 3-base DNA codon types."""
        if len(dna_sequence) % 3 != 0:
            raise ValueError("DNA sequence length must be a multiple of 3 for codon screening.")
            
        processed_codon_type_map = {}
        for type_name, codons_list in codon_type_map.items():
            valid_codons = []
            for codon in codons_list:
                codon_upper = codon.strip().upper()
                if len(codon_upper) == 3 and all(base in self.bases for base in codon_upper):
                    valid_codons.append(codon_upper)
                else:
                    print(f"Warning: Invalid codon '{codon}' for type '{type_name}' in screen_codons_in_dna. Skipping.")
            if valid_codons:
                processed_codon_type_map[type_name] = valid_codons
            else:
                print(f"No valid codons found for type '{type_name}'. This type will be skipped in screening.")

        found_codons_by_type = {type_name: [] for type_name in processed_codon_type_map.keys()}
        
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3].upper()
            for type_name, codons_list in processed_codon_type_map.items():
                if codon in codons_list:
                    found_codons_by_type[type_name].append((codon, i))
        
        return found_codons_by_type

    # --- New Biological Functions ---

    def _reverse_complement(self, dna_seq):
        """
        Generates the reverse complement of a DNA sequence.
        """
        dna_seq = dna_seq.upper()
        
        complement = "".join([self.complement_map.get(base, '') for base in dna_seq])
        return complement[::-1]

    def translate_dna_to_protein(self, dna_segment):
        """
        Translates a DNA segment into an amino acid sequence based on the
        standard genetic code.
        """
        if len(dna_segment) % 3 != 0:
            raise ValueError("DNA segment for translation must be a multiple of 3 bases.")
        
        amino_acids = []
        for i in range(0, len(dna_segment), 3):
            codon = dna_segment[i:i+3].upper()
            amino_acid = self.genetic_code.get(codon, 'X')
            amino_acids.append(amino_acid)
        
        return "".join(amino_acids)

    def find_orfs_and_proteins(self, dna_sequence, min_protein_length=10):
        """
        Finds Open Reading Frames (ORFs) and translates them into protein sequences
        from a given DNA sequence, considering all 6 reading frames.
        """
        dna_sequence = dna_sequence.upper()
        orfs_found = []

        strands = {
            "Forward": dna_sequence,
            "Reverse Complement": self._reverse_complement(dna_sequence)
        }

        for strand_name, current_strand_dna in strands.items():
            for frame in range(3):
                current_orf_dna = ""
                in_orf = False
                
                for i in range(frame, len(current_strand_dna) - 2, 3):
                    codon = current_strand_dna[i : i+3]

                    if not all(base in self.bases for base in codon):
                        if in_orf:
                            in_orf = False
                            translated_protein = self.translate_dna_to_protein(current_orf_dna)
                            if len(translated_protein.replace('*', '')) >= min_protein_length:
                                orfs_found.append({
                                    'strand': strand_name,
                                    'frame_start_index': frame,
                                    'start_pos_in_strand': i - len(current_orf_dna),
                                    'dna_sequence': current_orf_dna,
                                    'protein_sequence': translated_protein
                                })
                            current_orf_dna = ""
                        continue

                    amino_acid = self.genetic_code.get(codon, 'X')

                    if amino_acid == 'M' and codon in self.start_codons and not in_orf:
                        in_orf = True
                        current_orf_dna = codon
                    elif in_orf:
                        current_orf_dna += codon
                        if amino_acid == '*':
                            translated_protein = self.translate_dna_to_protein(current_orf_dna)
                            if len(translated_protein.replace('*', '')) >= min_protein_length:
                                orfs_found.append({
                                    'strand': strand_name,
                                    'frame_start_index': frame,
                                    'start_pos_in_strand': i - len(current_orf_dna) + 3,
                                    'dna_sequence': current_orf_dna,
                                    'protein_sequence': translated_protein
                                })
                            in_orf = False
                            current_orf_dna = ""
                
                if in_orf and current_orf_dna:
                    translated_protein = self.translate_dna_to_protein(current_orf_dna)
                    if len(translated_protein.replace('*', '')) >= min_protein_length:
                        orfs_found.append({
                            'strand': strand_name,
                            'frame_start_index': frame,
                            'start_pos_in_strand': (len(current_strand_dna) - len(current_orf_dna)),
                            'dna_sequence': current_orf_dna,
                            'protein_sequence': translated_protein
                        })
        return orfs_found

    def screen_proteins_for_motifs(self, protein_sequences_data, motif_map):
        """
        Screens a list of translated protein sequences (from ORFs) for user-defined
        amino acid motifs.
        """
        found_motifs_by_type = {type_name: [] for type_name in motif_map.keys()}
        
        processed_motif_map = {}
        for type_name, motifs_list in motif_map.items():
            valid_motifs = [m.strip().upper() for m in motifs_list if m.strip()]
            if valid_motifs:
                processed_motif_map[type_name] = valid_motifs
            else:
                print(f"Warning: No valid motifs found for type '{type_name}'. This type will be skipped.")

        for prot_idx, orf_data in enumerate(protein_sequences_data):
            protein_seq = orf_data['protein_sequence']
            
            for type_name, motifs_to_find in processed_motif_map.items():
                for motif in motifs_to_find:
                    start = 0
                    while True:
                        pos = protein_seq.find(motif, start)
                        if pos == -1:
                            break
                        found_motifs_by_type[type_name].append((motif, prot_idx, pos, protein_seq)) # Added protein_seq to output
                        start = pos + len(motif)
                        
        return found_motifs_by_type

    def optimize_codons(self, dna_sequence, target_organism_codon_usage=None):
        """
        Optimizes a DNA sequence for codon usage bias based on a target organism's
        preferred codons.
        
        Args:
            dna_sequence (str): The DNA sequence to optimize. Must be a multiple of 3.
            target_organism_codon_usage (dict): A dictionary mapping amino acids (1-letter code)
                                                to a list of preferred codons for the target organism.
                                                If None, uses the internal self.codon_usage (E. coli example).
        Returns:
            str: The codon-optimized DNA sequence.
            list: A list of warnings/notes regarding optimization.
        """
        if len(dna_sequence) % 3 != 0:
            raise ValueError("DNA sequence for codon optimization must be a multiple of 3 bases.")
        
        warnings = []
        optimized_dna_parts = []
        
        # Use provided codon usage or default
        current_codon_usage = target_organism_codon_usage if target_organism_codon_usage is not None else self.codon_usage

        for i in range(0, len(dna_sequence), 3):
            original_codon = dna_sequence[i:i+3].upper()
            amino_acid = self.genetic_code.get(original_codon, 'X') # 'X' for unknown/invalid codons

            if amino_acid == 'X':
                optimized_dna_parts.append(original_codon) # Keep original if unknown
                warnings.append(f"Warning: Unknown codon '{original_codon}' at position {i}. Skipping optimization for this codon.")
                continue

            # Get preferred codons for this amino acid
            preferred_codons = current_codon_usage.get(amino_acid)

            if preferred_codons:
                # If the original codon is already one of the preferred, keep it
                if original_codon in preferred_codons:
                    optimized_dna_parts.append(original_codon)
                else:
                    # Choose the first preferred codon as the optimized one
                    optimized_codon = preferred_codons[0]
                    optimized_dna_parts.append(optimized_codon)
                    warnings.append(f"Optimized '{original_codon}' (for {amino_acid}) to '{optimized_codon}' at position {i}.")
            else:
                optimized_dna_parts.append(original_codon) # Keep original if no preferred codon info
                warnings.append(f"No preferred codons found for amino acid '{amino_acid}'. Keeping original codon '{original_codon}' at position {i}.")
                
        return "".join(optimized_dna_parts), warnings

    def screen_for_toxic_proteins(self, protein_sequences_data, custom_toxic_motifs=None):
        """
        Screens a list of translated protein sequences for known or custom toxic protein motifs.
        
        Args:
            protein_sequences_data (list): List of dictionaries, each containing 'protein_sequence'.
            custom_toxic_motifs (dict, optional): A dictionary of custom motifs to screen for.
                                                  Key: motif type name, Value: list of motifs (strings).
                                                  If None, uses the internal self.toxic_protein_motifs.
                                                  
        Returns:
            dict: A dictionary of detected toxic motifs, structured by type.
                  Format: {'Type Name': [(motif, protein_index, start_position_in_protein, protein_sequence), ...]}
        """
        toxic_motifs_to_screen = custom_toxic_motifs if custom_toxic_motifs is not None else self.toxic_protein_motifs
        
        found_toxic_motifs = {type_name: [] for type_name in toxic_motifs_to_screen.keys()}

        for prot_idx, orf_data in enumerate(protein_sequences_data):
            protein_seq = orf_data['protein_sequence']
            
            for type_name, motifs_list in toxic_motifs_to_screen.items():
                for motif in motifs_list:
                    start = 0
                    while True:
                        pos = protein_seq.find(motif.upper(), start) # Ensure motif is uppercase for search
                        if pos == -1:
                            break
                        found_toxic_motifs[type_name].append((motif, prot_idx, pos, protein_seq))
                        start = pos + len(motif)
                        
        return found_toxic_motifs

# --- User Interface (CLI) ---
def main():
    """
    Main function to run the command-line interface for the DNAScrambler.
    Allows users to encrypt, decrypt, screen DNA codons, and perform biological analysis.
    """
    scrambler = DNAScrambler(data_symbols=16, parity_symbols=4)

    def get_multi_line_input(prompt_start, end_marker='END', is_dna_input=False):
        print(f"\n{prompt_start} (type '{end_marker}' on a new line to finish input):")
        lines = []
        while True:
            line = input()
            if line.upper() == end_marker:
                break
            lines.append(line)
        
        raw_input_string = "\n".join(lines)
        
        if is_dna_input:
            cleaned_dna = "".join(filter(str.isalpha, raw_input_string)).upper()
            return cleaned_dna
        else:
            return raw_input_string

    while True:
        print("\n--- DNA Scrambler & Bio-Analyzer ---")
        print("1. Encrypt ASCII Art to DNA (With Checksum)")
        print("2. Decrypt DNA sequence to ASCII Art (With Checksum)")
        print("3. Screen DNA for specific Codon Types (DNA Level Analysis)")
        print("4. Biological Analysis (Find ORFs & Screen Proteins)")
        print("5. Encrypt ASCII Art to DNA (With Reed-Solomon ECC)")
        print("6. Decrypt DNA sequence to ASCII Art (With Reed-Solomon ECC)")
        print("7. Codon Optimization")
        print("8. Toxic Protein Screening")
        print("9. Encrypt ASCII Art to DNA (ECC + Checksum)") # New Option 1
        print("10. Decrypt DNA sequence to ASCII Art (ECC + Checksum)") # New Option 2
        print("11. Exit") # Updated Exit option

        choice = input("Enter your choice (1-11): ") # Updated range

        if choice == '1':
            ascii_art_text = get_multi_line_input("Enter your ASCII Art")
            if not ascii_art_text.strip():
                print("No ASCII art entered. Encryption skipped.")
                continue
            
            print(f"\nOriginal ASCII Art (Character Count: {len(ascii_art_text)}):\n{ascii_art_text}")
            
            # --- RLE Compression ---
            compressed_text = scrambler._rle_compress(ascii_art_text)
            print(f"\nCompressed ASCII Art (RLE Character Count: {len(compressed_text)}):\n{compressed_text}")
            print(f"Compression Ratio: {len(ascii_art_text) / len(compressed_text):.2f}x")
            # --- End RLE Compression ---

            try:
                # Encrypt the compressed text with checksum
                encrypted_dna = scrambler.encrypt_with_checksum(compressed_text) 
                print(f"\nEncrypted DNA Sequence (Base Count: {len(encrypted_dna)}) [WITH CHECKSUM]:\n{encrypted_dna}")
            except ValueError as e:
                print(f"Encryption Error: {e}")
        
        elif choice == '2':
            dna_input = get_multi_line_input("Enter DNA sequence to decrypt (with checksum)", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for decryption. Skipping.")
                continue
            print(f"DNA Sequence for Decryption (Input Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            try:
                # Decrypt to RLE format with checksum verification
                decrypted_rle_text, warnings = scrambler.decrypt_with_checksum(dna_input)
                print(f"\nDecrypted RLE Text (Character Count: {len(decrypted_rle_text)}) [WITH CHECKSUM]:\n{decrypted_rle_text}")
                
                # --- RLE Decompression ---
                try:
                    decompressed_text = scrambler._rle_decompress(decrypted_rle_text)
                    print(f"\nDecompressed ASCII Art (Character Count: {len(decompressed_text)}):\n{decompressed_text}")
                except ValueError as e:
                    print(f"RLE Decompression Error: {e}. Displaying raw decrypted RLE text.")
                    decompressed_text = decrypted_rle_text
                # --- End RLE Decompression ---

                if warnings:
                    print("\n--- Decryption Warnings/Notes ---")
                    for warning in warnings:
                        print(f"- {warning}")
            except ValueError as e:
                print(f"Decryption Error: {e}")

        elif choice == '3':
            dna_input = get_multi_line_input("Enter DNA sequence to screen", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for screening. Skipping.")
                continue
            print(f"DNA Sequence for Screening (Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            print("\n--- Choose Codon Types for Screening ---")
            print("1. Standard Stop Codons (TAA, TAG, TGA)")
            print("2. Standard Start Codon (ATG)")
            print("3. Custom Codon Types")
            screen_choice = input("Enter your screening choice (1-3): ")
            codon_types_to_screen = {}
            if screen_choice == '1':
                codon_types_to_screen['Stop Codons'] = ["TAA", "TAG", "TGA"]
            elif screen_choice == '2':
                codon_types_to_screen['Start Codon'] = ["ATG"]
            elif screen_choice == '3':
                print("\n--- Define Custom Codon Types ---")
                print("Enter a type name (e.g., 'Cysteine Codons').")
                print("Then, enter the 3-base codons for that type, separated by commas (e.g., UGC,UGU).")
                print("Type 'done' as the type name when you are finished.")
                while True:
                    type_name = input("\nEnter codon type name (or 'done' to finish): ").strip()
                    if type_name.lower() == 'done':
                        break
                    codon_list_str = input(f"Enter codons for '{type_name}' (e.g., UGC,UGU): ").strip()
                    codons = [c.strip().upper() for c in codon_list_str.split(',') if c.strip()]
                    valid_custom_codons = []
                    for c in codons:
                        if len(c) == 3 and all(base in scrambler.bases for base in c):
                            valid_custom_codons.append(c)
                        else:
                            print(f"Warning: Invalid codon '{c}' provided. Skipping.")
                    if valid_custom_codons:
                        codon_types_to_screen[type_name] = valid_custom_codons
                    else:
                        print(f"No valid codons provided for type '{type_name}'. This type will be skipped.")
            else:
                print("Invalid screening choice. Please try again.")
                continue
            if not codon_types_to_screen:
                print("No codon types selected for screening.")
                continue
            try:
                screen_results = scrambler.screen_codons_in_dna(dna_input, codon_types_to_screen)
                print("\n--- Codon Screening Results ---")
                found_any_result = False
                for type_name, findings in screen_results.items():
                    if findings:
                        found_any_result = True
                        print(f"\n{type_name} Found:")
                        for codon, pos in findings:
                            print(f"  - Codon: '{codon}' at DNA position: {pos}")
                if not found_any_result:
                    print("No specified codon types found in the DNA sequence.")
            except ValueError as e:
                print(f"Screening Error: {e}")

        elif choice == '4':
            dna_input = get_multi_line_input("Enter DNA sequence for biological analysis", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for biological analysis. Skipping.")
                continue
            print(f"DNA Sequence for Biological Analysis (Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            try:
                min_len_str = input("Enter minimum protein length to report (e.g., 10 for 10 amino acids): ")
                min_protein_length = int(min_len_str)
                if min_protein_length < 1:
                    print("Minimum length must be at least 1. Setting to 10.")
                    min_protein_length = 10
            except ValueError:
                print("Invalid length. Setting minimum protein length to 10.")
                min_protein_length = 10
            print("\n--- Performing Biological Analysis ---")
            print("Finding Open Reading Frames (ORFs) and translating to protein sequences...")
            try:
                orfs_data = scrambler.find_orfs_and_proteins(dna_input, min_protein_length)
                if not orfs_data:
                    print("No ORFs/potential proteins found meeting the minimum length criteria.")
                else:
                    print(f"\nFound {len(orfs_data)} potential protein sequences:")
                    for i, orf in enumerate(orfs_data):
                        print(f"\n--- Protein {i+1} ---")
                        print(f"  Strand: {orf['strand']}")
                        print(f"  Frame Start Index: {orf['frame_start_index']}")
                        print(f"  Start Position in Strand: {orf['start_pos_in_strand']}")
                        print(f"  DNA Sequence ({len(orf['dna_sequence'])} bases): {orf['dna_sequence']}")
                        print(f"  Protein Sequence ({len(orf['protein_sequence'])} AAs): {orf['protein_sequence']}")
                if orfs_data:
                    perform_motif_screening = input("\nDo you want to screen these proteins for custom amino acid motifs? (yes/no): ").lower()
                    if perform_motif_screening == 'yes':
                        print("\n--- Define Custom Protein Motifs ---")
                        print("Enter a motif type name (e.g., 'Zinc Finger').")
                        print("Then, enter amino acid motifs (e.g., CXXC, HHHH). Separate with commas.")
                        print("Amino acids are 1-letter codes (e.g., M, V, S, *, X).")
                        print("Type 'done' as the type name when you are finished.")
                        protein_motifs_to_screen = {}
                        while True:
                            type_name = input("\nEnter motif type name (or 'done' to finish): ").strip()
                            if type_name.lower() == 'done':
                                break
                            motif_list_str = input(f"Enter motifs for '{type_name}' (e.g., CXXC, HHHH): ").strip()
                            motifs = [m.strip().upper() for m in motif_list_str.split(',') if m.strip()]
                            valid_motifs = []
                            for m in motifs:
                                # Fix: Check if all chars in motif are valid AA codes or 'X'
                                if all(char in scrambler.genetic_code.values() or char == 'X' for char in m):
                                    valid_motifs.append(m)
                                else:
                                    print(f"Warning: Invalid motif '{m}' for type '{type_name}'. Contains unrecognized amino acid codes. Skipping.")
                            if valid_motifs:
                                protein_motifs_to_screen[type_name] = valid_motifs
                            else:
                                print(f"No valid motifs provided for type '{type_name}'. This type will be skipped.")
                        if protein_motifs_to_screen:
                            print("\n--- Protein Motif Screening Results ---")
                            motif_screen_results = scrambler.screen_proteins_for_motifs(orfs_data, protein_motifs_to_screen)
                            found_any_motif = False
                            for type_name, findings in motif_screen_results.items():
                                if findings:
                                    found_any_motif = True
                                    print(f"\n{type_name} Motifs Found:")
                                    for motif, prot_idx, pos, prot_seq in findings: # Added prot_seq
                                        print(f"  - Motif: '{motif}' in Protein {prot_idx+1} (AA pos: {pos})")
                                        print(f"    Full Protein Sequence: {prot_seq}") # Display full protein
                            if not found_any_motif:
                                print("No specified protein motifs found in the translated sequences.")
                        else:
                            print("No protein motifs selected for screening.")
                else:
                    print("No proteins found to screen for motifs.")
            except ValueError as e:
                print(f"Biological Analysis Error: {e}")

        elif choice == '5':
            ascii_art_text = get_multi_line_input("Enter ASCII Art to encrypt with ECC")
            if not ascii_art_text.strip():
                print("No ASCII art entered. Encryption skipped.")
                continue

            print(f"\nOriginal ASCII Art (Character Count: {len(ascii_art_text)}):\n{ascii_art_text}")
            
            # --- RLE Compression ---
            compressed_text = scrambler._rle_compress(ascii_art_text)
            print(f"\nCompressed ASCII Art (RLE Character Count: {len(compressed_text)}):\n{compressed_text}")
            print(f"Compression Ratio: {len(ascii_art_text) / len(compressed_text):.2f}x")
            # --- End RLE Compression ---

            try:
                # Encrypt the compressed text
                encrypted_dna = scrambler.encrypt_with_ecc(compressed_text)
                print(f"\nEncrypted DNA Sequence (Base Count: {len(encrypted_dna)}) [WITH ECC]:\n{encrypted_dna}")
            except ValueError as e:
                print(f"ECC Encryption Error: {e}")
            except RuntimeError as e: # Catch the specific error for insufficient codon pairs
                print(f"ECC Encryption Configuration Error: {e}")


        elif choice == '6':
            dna_input = get_multi_line_input("Enter DNA sequence to decrypt with ECC", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for decryption. Skipping.")
                continue
            print(f"DNA Sequence for Decryption (Input Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            try:
                # Decrypt to RLE format
                decrypted_rle_text, warnings = scrambler.decrypt_with_ecc(dna_input)
                print(f"\nDecrypted RLE Text (Character Count: {len(decrypted_rle_text)}) [WITH ECC]:\n{decrypted_rle_text}")
                
                # --- RLE Decompression ---
                try:
                    decompressed_text = scrambler._rle_decompress(decrypted_rle_text)
                    print(f"\nDecompressed ASCII Art (Character Count: {len(decompressed_text)}):\n{decompressed_text}")
                except ValueError as e:
                    print(f"RLE Decompression Error: {e}. Displaying raw decrypted RLE text.")
                    decompressed_text = decrypted_rle_text
                # --- End RLE Decompression ---

                if warnings:
                    print("\n--- Decryption Warnings/Notes ---")
                    for warning in warnings:
                        print(f"- {warning}")
            except ValueError as e:
                print(f"ECC Decryption Error: {e}")
                
        elif choice == '7':
            dna_input = get_multi_line_input("Enter DNA sequence to optimize codons", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for optimization. Skipping.")
                continue
            
            print(f"Original DNA Sequence for Optimization (Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)

            print("\n--- Codon Optimization Options ---")
            print("1. Use default (E. coli example) codon usage table")
            print("2. Provide custom codon usage table (Amino Acid: Codons)")
            opt_choice = input("Enter your choice (1-2): ")

            custom_codon_usage = None
            if opt_choice == '2':
                custom_codon_usage = {}
                print("\n--- Enter Custom Codon Usage ---")
                print("For each amino acid (1-letter code, e.g., 'A' for Alanine), enter its preferred codons, separated by commas.")
                print("Example: A: GCT,GCC,GCG. Type 'done' for amino acid when finished with a particular amino acid, or 'exit' to stop.")
                while True:
                    aa_input = input("Enter Amino Acid (e.g., 'A', or 'done' to finish for this AA, 'exit' to stop entering tables): ").strip().upper()
                    if aa_input == 'EXIT':
                        break
                    if aa_input == 'DONE':
                        continue # Allow "done" to just go to next AA
                    
                    if len(aa_input) != 1 or not aa_input.isalpha():
                        print("Invalid amino acid code. Please enter a single letter.")
                        continue
                    
                    codons_str = input(f"Enter preferred codons for '{aa_input}' (e.g., TTT,TTC): ").strip()
                    codons = [c.strip().upper() for c in codons_str.split(',') if c.strip()]
                    
                    valid_aa_codons = []
                    for c in codons:
                        if len(c) == 3 and all(base in scrambler.bases for base in c) and scrambler.genetic_code.get(c) == aa_input:
                            valid_aa_codons.append(c)
                        else:
                            print(f"Warning: Invalid codon '{c}' for amino acid '{aa_input}' or codon does not map to this amino acid. Skipping.")
                    
                    if valid_aa_codons:
                        custom_codon_usage[aa_input] = valid_aa_codons
                    else:
                        print(f"No valid codons provided for amino acid '{aa_input}'. Skipping this amino acid.")
            
            try:
                optimized_dna, opt_warnings = scrambler.optimize_codons(dna_input, custom_codon_usage)
                print(f"\nOptimized DNA Sequence (Base Count: {len(optimized_dna)}):\n{optimized_dna}")
                if opt_warnings:
                    print("\n--- Optimization Warnings/Notes ---")
                    for warning in opt_warnings:
                        print(f"- {warning}")
            except ValueError as e:
                print(f"Codon Optimization Error: {e}")

        elif choice == '8':
            dna_input = get_multi_line_input("Enter DNA sequence for toxic protein screening", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for screening. Skipping.")
                continue
            print(f"DNA Sequence for Toxic Protein Screening (Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            
            try:
                min_len_str = input("Enter minimum protein length to consider for ORFs (e.g., 10 amino acids): ")
                min_protein_length = int(min_len_str)
                if min_protein_length < 1:
                    print("Minimum length must be at least 1. Setting to 10.")
                    min_protein_length = 10
            except ValueError:
                print("Invalid length. Setting minimum protein length to 10.")
                min_protein_length = 10

            print("\n--- Finding ORFs and translating proteins for toxic screening ---")
            orfs_data = scrambler.find_orfs_and_proteins(dna_input, min_protein_length)
            
            if not orfs_data:
                print("No ORFs/potential proteins found meeting the minimum length criteria. Cannot screen for toxic proteins.")
                continue

            print(f"\nFound {len(orfs_data)} potential protein sequences to screen.")
            
            screen_toxic_choice = input("\nDo you want to use custom toxic motifs or default motifs? (custom/default): ").lower()
            toxic_motifs_config = None
            if screen_toxic_choice == 'custom':
                toxic_motifs_config = {}
                print("\n--- Define Custom Toxic Protein Motifs ---")
                print("Enter a motif type name (e.g., 'Transmembrane Domain').")
                print("Then, enter amino acid motifs (e.g., AAAA, LLLV). Separate with commas.")
                print("Amino acids are 1-letter codes.")
                print("Type 'done' as the type name when you are finished.")
                while True:
                    type_name = input("\nEnter motif type name (or 'done' to finish): ").strip()
                    if type_name.lower() == 'done':
                        break
                    motif_list_str = input(f"Enter motifs for '{type_name}' (e.g., AAAA, LLLV): ").strip()
                    motifs = [m.strip().upper() for m in motif_list_str.split(',') if m.strip()]
                    if motifs:
                        toxic_motifs_config[type_name] = motifs
                    else:
                        print(f"No motifs provided for type '{type_name}'. Skipping this type.")
            elif screen_toxic_choice != 'default':
                print("Invalid choice. Using default toxic motifs.")

            try:
                toxic_screen_results = scrambler.screen_for_toxic_proteins(orfs_data, toxic_motifs_config)
                print("\n--- Toxic Protein Screening Results ---")
                found_any_toxic_motif = False
                for type_name, findings in toxic_screen_results.items():
                    if findings:
                        found_any_toxic_motif = True
                        print(f"\n{type_name} Found:")
                        for motif, prot_idx, pos, prot_seq in findings:
                            print(f"  - Motif: '{motif}' in Protein {prot_idx+1} (AA pos: {pos})")
                            print(f"    Full Protein Sequence: {prot_seq}")
                if not found_any_toxic_motif:
                    print("No specified toxic protein motifs found in the translated sequences.")
            except ValueError as e:
                print(f"Toxic Protein Screening Error: {e}")

        elif choice == '9': # NEW Encrypt with ECC + Checksum
            ascii_art_text = get_multi_line_input("Enter ASCII Art to encrypt with ECC and Checksum")
            if not ascii_art_text.strip():
                print("No ASCII art entered. Encryption skipped.")
                continue

            print(f"\nOriginal ASCII Art (Character Count: {len(ascii_art_text)}):\n{ascii_art_text}")
            
            # --- RLE Compression ---
            compressed_text = scrambler._rle_compress(ascii_art_text)
            print(f"\nCompressed ASCII Art (RLE Character Count: {len(compressed_text)}):\n{compressed_text}")
            print(f"Compression Ratio: {len(ascii_art_text) / len(compressed_text):.2f}x")
            # --- End RLE Compression ---

            try:
                encrypted_dna = scrambler.encrypt_with_ecc_and_checksum(compressed_text)
                print(f"\nEncrypted DNA Sequence (Base Count: {len(encrypted_dna)}) [ECC + CHECKSUM]:\n{encrypted_dna}")
            except ValueError as e:
                print(f"Encryption Error (ECC + Checksum): {e}")
            except RuntimeError as e:
                print(f"Encryption Configuration Error (ECC + Checksum): {e}")

        elif choice == '10': # NEW Decrypt with ECC + Checksum
            dna_input = get_multi_line_input("Enter DNA sequence to decrypt (ECC + Checksum)", is_dna_input=True)
            if not dna_input:
                print("No DNA sequence entered for decryption. Skipping.")
                continue
            print(f"DNA Sequence for Decryption (Input Base Count: {len(dna_input)} bases cleaned from input):")
            if len(dna_input) > 200:
                print(f"{dna_input[:100]}...{dna_input[-100:]}")
            else:
                print(dna_input)
            try:
                decrypted_rle_text, warnings = scrambler.decrypt_with_ecc_and_checksum(dna_input)
                print(f"\nDecrypted RLE Text (Character Count: {len(decrypted_rle_text)}) [ECC + CHECKSUM]:\n{decrypted_rle_text}")
                
                # --- RLE Decompression ---
                try:
                    decompressed_text = scrambler._rle_decompress(decrypted_rle_text)
                    print(f"\nDecompressed ASCII Art (Character Count: {len(decompressed_text)}):\n{decompressed_text}")
                except ValueError as e:
                    print(f"RLE Decompression Error: {e}. Displaying raw decrypted RLE text.")
                    decompressed_text = decrypted_rle_text
                # --- End RLE Decompression ---

                if warnings:
                    print("\n--- Decryption Warnings/Notes (ECC + Checksum) ---")
                    for warning in warnings:
                        print(f"- {warning}")
            except ValueError as e:
                print(f"Decryption Error (ECC + Checksum): {e}")

        elif choice == '11': # Updated Exit option
            print("Exiting the program. Goodbye!")
            break
            
        else:
            print("Invalid choice. Please enter a number between 1 and 11.") # Updated range

if __name__ == "__main__":
    main()

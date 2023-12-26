import tkinter as tk
from tkinter import messagebox

class ReedSolomonApp:
    def __init__(self, master):
        self.master = master
        master.title("Reed-Solomon Encoder/Decoder")

        self.message_label = tk.Label(master, text="Enter Message:")
        self.message_label.pack()

        self.message_entry = tk.Entry(master, width=30)
        self.message_entry.pack()

        self.nsym_label = tk.Label(master, text="Number of Error Correction Symbols:")
        self.nsym_label.pack()

        self.nsym_entry = tk.Entry(master, width=5)
        self.nsym_entry.pack()

        self.encode_button = tk.Button(master, text="Encode", command=self.encode_message)
        self.encode_button.pack()

        self.error_label = tk.Label(master, text="Introduce Error at Position:")
        self.error_label.pack()

        self.error_entry = tk.Entry(master, width=5)
        self.error_entry.pack()

        self.decode_button = tk.Button(master, text="Decode", command=self.decode_message)
        self.decode_button.pack()

        self.original_label = tk.Label(master, text="Original Message:")
        self.original_label.pack()

        self.original_message_text = ""
        self.codeword = bytearray()
        self.decoded_message_text = ""

    def encode_message(self):
        try:
            self.original_message_text = self.message_entry.get()
            self.nsym = int(self.nsym_entry.get())
            self.codeword = encode_rs(self.original_message_text.encode('utf-8'), self.nsym)
            messagebox.showinfo("Encoded Message", self.codeword)
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter valid values.")

    def decode_message(self):
        try:
            codeword_with_errors = bytearray(self.codeword)
            error_position = int(self.error_entry.get())
            codeword_with_errors[error_position] ^= 0x01
            decoded_message = decode_rs(codeword_with_errors, self.nsym)
            self.decoded_message_text = decoded_message.decode('utf-8')
            self.original_label.config(text="Original Message: " + self.original_message_text)
        except ValueError as e:
            messagebox.showerror("Error", str(e))

def encode_rs(message, nsym):
    codeword = bytearray(message)
    codeword.extend([0] * nsym)

    for i in range(len(message)):
        factor = codeword[i]
        if factor != 0:
            for j in range(nsym):
                codeword[i + j] ^= GF_MUL(factor, GF_EXP[j])

    return codeword

def decode_rs(codeword, nsym):
    syndromes = calculate_syndromes(codeword, nsym)
    error_locator, _ = euclidean_algorithm(syndromes, nsym)
    error_positions = find_error_positions(error_locator, nsym)
    
    if len(error_positions) > nsym:
        raise ValueError("Too many errors to correct")

    error_values = find_errors_chien(codeword, nsym, error_positions)
    correct_errors_forney(codeword, nsym, error_locator, error_positions, error_values)

    return codeword[:len(codeword) - nsym]

def calculate_syndromes(codeword, nsym):
    syndromes = [0] * nsym
    for i in range(nsym):
        for j in range(len(codeword)):
            syndromes[i] ^= GF_MUL(codeword[j], GF_EXP[i * j % 255])
    return syndromes

def euclidean_algorithm(syndromes, nsym):
    a = [1] + [0] * (nsym + 1)
    b = [0] + [1] + [0] * nsym
    for i in range(nsym):
        temp_b = b[-1]
        for j in range(nsym + 1):
            b[j] = b[j] - temp_b * syndromes[i] * a[j]
        a.pop(0)
        a.append(0)
    return b, syndromes

def find_error_positions(error_locator, nsym):
    error_positions = []
    for i in range(1, nsym + 1):
        root = GF_EXP[(255 - GF_LOG[i]) % 255]
        if evaluate_polynomial(error_locator, root) == 0:
            error_positions.append(nsym - i)
    return error_positions

def find_errors_chien(codeword, nsym, error_positions):
    error_values = []
    for i in range(len(error_positions)):
        root = GF_EXP[(255 - GF_LOG[error_positions[i]]) % 255]
        derivative = evaluate_derivative(codeword, root)
        error_value = GF_DIV(evaluate_polynomial(codeword, root), derivative)
        error_values.append(error_value)
    return error_values

def correct_errors_forney(codeword, nsym, error_locator, error_positions, error_values):
    for i in range(len(error_positions)):
        root = GF_EXP[(255 - GF_LOG[error_positions[i]]) % 255]
        error_position = nsym - error_positions[i] - 1
        for j in range(len(error_locator) - 1):
            codeword[j] ^= GF_MUL(error_locator[j + 1], error_values[i]) * GF_INV(root)

GF_EXP = [0] * 512
GF_LOG = [0] * 256

def initialize_galois_field():
    prim_poly = 0x11d
    x = 1
    for i in range(255):
        GF_EXP[i] = x
        GF_LOG[x] = i
        x <<= 1
        if x & 0x100:
            x ^= prim_poly

initialize_galois_field()

def GF_MUL(a, b):
    if a == 0 or b == 0:
        return 0
    return GF_EXP[(GF_LOG[a] + GF_LOG[b]) % 255]

def GF_INV(a):
    return GF_EXP[255 - GF_LOG[a]]

def GF_DIV(a, b):
    if b == 0:
        raise ZeroDivisionError("Division by zero")
    if a == 0:
        return 0
    return GF_MUL(a, GF_INV(b))

def evaluate_polynomial(coeffs, x):
    result = 0
    for i in range(len(coeffs)):
        result = GF_MUL(result, x) ^ coeffs[i]
    return result

def evaluate_derivative(coeffs, x):
    result = 0
    for i in range(1, len(coeffs)):
        result = GF_MUL(result, x) ^ GF_MUL(coeffs[i], i)
    return result

if __name__ == "__main__":
    root = tk.Tk()
    app = ReedSolomonApp(root)
    root.mainloop()
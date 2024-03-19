import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np
from scipy.interpolate import LinearNDInterpolator

# Global variables to store loaded data
E_11, E_33, nu_13, voids, theta, pores_voids = [], [], [], [], [], []

# Global variables for interpolators
Voids_interpolator = None
Theta_t_interpolator = None
Pores_Voids_Ratio_interpolator = None

# Load data from CSV file function
def load_data(file_path):
    global E_11, E_33, nu_13, voids, theta, pores_voids, Voids_interpolator, Theta_t_interpolator, Pores_Voids_Ratio_interpolator
    results = np.genfromtxt(file_path, delimiter=",")[1:,:]
    E_11 = np.divide(results[:, 4],70000)
    E_33 =  np.divide(results[:, 5],70000)
    nu_13 = results[:, 6]
    voids = results[:, 0]
    theta = np.multiply(results[:, 1], 2)
    pores_voids = results[:, 2]

    # Build interpolators
    Voids_interpolator = LinearNDInterpolator((E_11, E_33, nu_13), voids)
    Theta_t_interpolator = LinearNDInterpolator((E_11, E_33, nu_13), theta)
    Pores_Voids_Ratio_interpolator = LinearNDInterpolator((E_11, E_33, nu_13), pores_voids)

# Function to calculate and display output
def calculate_output():
    if not E_11.size or not E_33.size or not nu_13.size or not voids.size or not theta.size or not pores_voids.size:
        # Data not loaded, prompt the user to load a file
        result_label.config(text="Please load a CSV file first.")
        return

    E_11_val = float(entry_E_11.get())/float(entry_bulk_modulus.get())
    E_33_val = float(entry_E_33.get())/float(entry_bulk_modulus.get())
    nu_13_val = float(entry_nu_13.get())

    voids_output = Voids_interpolator(E_11_val, E_33_val, nu_13_val) * 100
    theta_t_output = Theta_t_interpolator(E_11_val, E_33_val, nu_13_val)
    pores_voids_output = Pores_Voids_Ratio_interpolator(E_11_val, E_33_val, nu_13_val)

    result_label.config(text=f"Voids: {voids_output:.2f}%\nθₜ: {theta_t_output:.2f}°\nPores/Voids Ratio: {pores_voids_output:.2f}")

# Function to browse and select a CSV file
def browse_file():
    file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
    entry_output_file.delete(0, tk.END)
    entry_output_file.insert(0, file_path)

    # Load data from the selected CSV file
    load_data(file_path)

# GUI setup
root = tk.Tk()
root.title("Defect Descriptors Calculator")

# Labels and entry widgets for input
label_output_file = ttk.Label(root, text="Select CSV File:")
label_output_file.grid(row=0, column=0, padx=5, pady=5)
entry_output_file = ttk.Entry(root)
entry_output_file.grid(row=0, column=1, columnspan=2, padx=5, pady=5)
browse_button = ttk.Button(root, text="Browse", command=browse_file)
browse_button.grid(row=0, column=3, padx=5, pady=5)

# New input for Bulk Material Modulus
label_bulk_modulus = ttk.Label(root, text="Bulk Material Modulus (MPa):")
label_bulk_modulus.grid(row=1, column=0, padx=5, pady=5)
entry_bulk_modulus = ttk.Entry(root)
entry_bulk_modulus.grid(row=1, column=1, padx=5, pady=5)

# Labels and entry widgets for other inputs
label_E_11 = ttk.Label(root, text="E₁₁ (MPa):")
label_E_11.grid(row=2, column=0, padx=5, pady=5)
entry_E_11 = ttk.Entry(root)
entry_E_11.grid(row=2, column=1, padx=5, pady=5)

label_E_33 = ttk.Label(root, text="E₃₃ (MPa):")
label_E_33.grid(row=3, column=0, padx=5, pady=5)
entry_E_33 = ttk.Entry(root)
entry_E_33.grid(row=3, column=1, padx=5, pady=5)

label_nu_13 = ttk.Label(root, text="ν₁₃:")
label_nu_13.grid(row=4, column=0, padx=5, pady=5)
entry_nu_13 = ttk.Entry(root)
entry_nu_13.grid(row=4, column=1, padx=5, pady=5)

# Button to calculate and display output
calculate_button = ttk.Button(root, text="Calculate", command=calculate_output)
calculate_button.grid(row=5, column=0, columnspan=4, pady=10)

# Label to display output
result_label = ttk.Label(root, text="")
result_label.grid(row=6, column=0, columnspan=4, pady=10)

root.mainloop()


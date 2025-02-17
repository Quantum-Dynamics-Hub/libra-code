# *********************************************************************************
# * Copyright (C) 2024 Qingxin Zhang, Alexey V. Akimov
# *
# * This file is distributed under the terms of the GNU General Public License
# * as published by the Free Software Foundation, either version 3 of
# * the License, or (at your option) any later version.
# * See the file LICENSE in the root directory of this distribution
# * or <http://www.gnu.org/licenses/>.
# *
# *********************************************************************************/
#
#
#
"""

.. module:: plotting
   :platform: Unix, Windows
   :synopsis: This module implements function for plotting the UV-Vis spectrums of CP2K TD-DFT calculations
       Example:
.. moduleauthor:: Qingxin Zhang, Alexey V. Akimov

"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os


def parse_spectrum_data_from_log(file_path):
    energy_levels, intensities = [], []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(' TDDFPT|') and len(line.split()) >= 7:
                parts = line.split()
                try:
                    energy_levels.append(float(parts[2]))
                    intensities.append(float(parts[-1]))
                except ValueError:
                    continue
    return np.array(energy_levels), np.array(intensities)


def gaussian_broadening(energy_levels, intensities, fwhm, num_points=1000):
    if len(energy_levels) == 0 or len(intensities) == 0:
        raise ValueError("Energy levels or intensities array is empty. Check input data.")
    x = np.linspace(np.min(energy_levels) - 5 * fwhm, np.max(energy_levels) + 5 * fwhm, num_points)
    y = np.sum([intensity * np.exp(-((x - energy) ** 2) / (2 * (fwhm / 2.35482) ** 2))
                for energy, intensity in zip(energy_levels, intensities)], axis=0)
    return x, y


def generate_and_save_spectra(log_file_pattern, output_folder, fwhm=0.1, num_points=1000):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    log_files = glob.glob(log_file_pattern)
    if not log_files:
        raise FileNotFoundError(f"No files matched the pattern: {log_file_pattern}")

    for log_file in log_files:
        print(f"Processing: {log_file}")
        energy_levels, intensities = parse_spectrum_data_from_log(log_file)
        x, y = gaussian_broadening(energy_levels, intensities, fwhm, num_points)
        base_name = os.path.basename(log_file).replace('.log', '_spectrum.txt')
        output_file_path = os.path.join(output_folder, base_name)
        np.savetxt(output_file_path, np.column_stack((x, y)), header="Energy (eV)\tIntensity (a.u.)")
        print(f"Spectrum data saved to: {output_file_path}")


def process_spectra(log_file_pattern, output_folder, fwhm=0.1, num_points=1000):
    generate_and_save_spectra(log_file_pattern, output_folder, fwhm, num_points)
    data_dir = output_folder
    if not os.path.exists(data_dir):
        print(f"Catalog {data_dir} not exist, please check path!")
        return

    file_list = [f for f in os.listdir(data_dir) if f.startswith("step_") and f.endswith("_spectrum.txt")]
    if not file_list:
        print(f"No matching file was found, please check the file name or path!")
        return

    file_list = [os.path.join(data_dir, f) for f in file_list]
    x_values_list = []
    y_values_list = []

    for file_name in file_list:
        try:
            print(f"Processing file: {file_name}")
            data = np.loadtxt(file_name)
            x_values_list.append(data[:, 0])
            y_values_list.append(data[:, 1])
        except (FileNotFoundError, ValueError) as e:
            print(f"Error processing {file_name}: {e}")

    if x_values_list and y_values_list:
        x_values_avg = np.mean(x_values_list, axis=0)
        y_values_avg = np.mean(y_values_list, axis=0)
        average_data = np.column_stack((x_values_avg, y_values_avg))
        np.savetxt('average_spectrum_data.txt', average_data, fmt='%.6f', header='X_avg Y_avg', comments='')
        plt.figure(figsize=(10, 6))
        plt.plot(x_values_avg, y_values_avg, label='Average Spectrum Curve', color='blue')
        plt.xlabel('Energy (eV)', fontsize=14)
        plt.ylabel('Intensity (a.u.)', fontsize=14)
        plt.title('Average Spectrum Curve', fontsize=16)
        plt.legend(fontsize=12)
        plt.tight_layout()
        plt.savefig('average_spectrum.png', dpi=300)
        plt.show()
    else:
        print("No valid data file found")

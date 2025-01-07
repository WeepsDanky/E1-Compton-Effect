#!/usr/bin/env python3
"""
Mini-Study 2: Detector Refinement
- We include a more realistic detector model: Gaussian energy smearing,
  possibly a small shift due to the PMT window, etc.
- We still assume single-scatter at fixed angles.
"""

import numpy as np
import matplotlib.pyplot as plt
from compton_helpers import (E0, sample_theta_klein_nishina,
                             compton_scattered_energy,
                             gaussian_energy_smearing,
                             plot_spectrum)

def mini_study_2_detector(num_photons=100000, angles_deg=[30, 90, 135], 
                          FWHM_at_662=0.07, PMT_offset=2.0):
    """
    :param FWHM_at_662: e.g. 0.07 => 7% FWHM at 662 keV. Adjust to match your detector.
    :param PMT_offset: keV offset to simulate, e.g. energy lost in a PMT window.
    """
    angles_rad = np.radians(angles_deg)
    results = {}

    for angle_deg, angle_rad in zip(angles_deg, angles_rad):
        # Ideal Compton scattered energy (single scatter)
        E_scattered = compton_scattered_energy(E0, angle_rad)
        
        # Again, let's create a small distribution around that energy 
        # to mimic partial physics, but now we apply a realistic resolution 
        # function (Gaussian).
        tiny_spread = 0.01 * E_scattered
        energies = np.random.uniform(E_scattered - tiny_spread, 
                                     E_scattered + tiny_spread, 
                                     num_photons)
        
        # PMT window offset: let's assume each photon effectively loses some keV 
        # before detection, or the calibration shifts the measured peak.
        energies_measured = energies - PMT_offset
        energies_measured[energies_measured < 0] = 0.0
        
        # Apply Gaussian smearing for detector resolution
        energies_smeared = gaussian_energy_smearing(energies_measured, 
                                                    FWHM_at_662=FWHM_at_662)
        
        results[angle_deg] = energies_smeared

    # Plot comparison
    plt.figure(figsize=(8,6))
    for angle_deg in angles_deg:
        energies = results[angle_deg]
        plot_spectrum(energies,
                      bins=80,
                      label=f"{angle_deg} deg - Refined Detector")
    
    plt.xlabel("Energy (keV)")
    plt.ylabel("Normalized Counts")
    plt.title("Mini-Study 2: Detector Refinement")
    plt.legend()
    plt.show()
    
    return results

if __name__ == "__main__":
    mini_study_2_detector()

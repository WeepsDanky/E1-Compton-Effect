#!/usr/bin/env python3
"""
Mini-Study 1: Baseline Simulation
- We simulate single-scattering events from Cs-137 (662 keV) at given angles.
- We ignore detector resolution or any detailed geometry.
- Output: Compare the raw scattered energy distribution to experimental data.
"""

import numpy as np
import matplotlib.pyplot as plt
from compton_helpers import (E0, sample_theta_klein_nishina,
                             compton_scattered_energy,
                             plot_spectrum)

def mini_study_1_baseline(num_photons=100000, angles_deg=[30, 90, 135]):
    """
    Baseline simulation: no detector smearing, no geometry correction.
    """
    # Convert angles to radians
    angles_rad = np.radians(angles_deg)
    
    # We'll store results in a dictionary if we want to compare later
    results = {}
    
    # We'll do a simple approach: directly use the known angle to compute scattered energy
    # (i.e., we do not randomize angle from the Klein-Nishina distribution;
    #  we are specifically simulating a scattering event at a fixed angle).
    #
    # Alternatively, you could random-sample from the full distribution of angles 
    # if your experiment measures a broad geometry. But let's do "fixed angle" to 
    # mimic a single-scatter experiment at each angle.

    for angle_deg, angle_rad in zip(angles_deg, angles_rad):
        # We simulate num_photons all scattering at 'angle_rad'
        E_scattered = compton_scattered_energy(E0, angle_rad)
        
        # In a real experiment, you get a distribution of scattered energies 
        # from a peak. But if we assume purely single-scatter with no resolution,
        # it's basically a delta function. So let's artificially broaden it by 
        # sampling a small spread to mimic partial physics or multiple scattering:
        
        # For baseline, let's just do a minuscule uniform spread to produce a tiny distribution.
        tiny_spread = 0.01 * E_scattered
        energies = np.random.uniform(E_scattered - tiny_spread, 
                                     E_scattered + tiny_spread, 
                                     num_photons)
        
        # Save
        results[angle_deg] = energies
        
    # Plot each angle's spectrum
    plt.figure(figsize=(8,6))
    for angle_deg in angles_deg:
        energies = results[angle_deg]
        plot_spectrum(energies,
                      bins=50,
                      label=f"{angle_deg} deg - Baseline")
    
    plt.xlabel("Energy (keV)")
    plt.ylabel("Normalized Counts")
    plt.title("Mini-Study 1: Baseline Simulation")
    plt.legend()
    plt.show()
    
    return results

if __name__ == "__main__":
    mini_study_1_baseline()

#!/usr/bin/env python3
"""
Mini-Study 3 (Extended): Geometry Refinement with Lead Pass-Through
-------------------------------------------------------------------
Now includes a simplified model for internal backscatter in the NaI(Tl) detector.

Key changes:
- Part F: Each photon has a probability `p_backscatter` of doing one Compton 
  scatter in the iodine. A fraction of those might "escape", leaving a 
  backscatter peak in the deposited energy.
"""

import numpy as np
import matplotlib.pyplot as plt

# Same helpers as before:
from compton_helpers import (
    E0,
    sample_theta_klein_nishina,      # for random scattering angles
    compton_scattered_energy,        # Compton formula
    gaussian_energy_smearing,        # Detector resolution
    plot_spectrum                    # Quick histogram plotting
)

def simulate_iodine_backscatter(energies_in, p_backscatter=0.1, p_escape=0.5):
    """
    Simplified model of internal backscatter in the NaI detector.

    :param energies_in: 1D array of photon energies (keV) arriving in the detector.
    :param p_backscatter: Probability that a photon undergoes ONE Compton scatter 
                          in the iodine crystal.
    :param p_escape: Probability that the scattered photon escapes the crystal, 
                     leaving only partial energy deposit.

    Returns a 1D array of "detected" energies, i.e., the energy deposit after 
    possible Compton scatter + escape.
    """

    # Copy the input energies so we can modify them
    energies_out = np.array(energies_in, copy=True)
    n_photons = len(energies_out)

    # Determine which photons will attempt a backscatter
    backscatter_flags = np.random.rand(n_photons) < p_backscatter

    # For those that do backscatter, pick a scattering angle from the KN distribution:
    n_bs = np.count_nonzero(backscatter_flags)
    if n_bs > 0:
        # Random angles for those that backscatter
        thetas = sample_theta_klein_nishina(E0, size=n_bs)  
        # ^ We can refine by using energies_out as the incident energy, but 
        #   for simplicity let's assume all near 662 keV or the "current" 
        #   photon energy. A more advanced approach would do it per-photon.

        # Compute the scattered photon energy
        # If we want to be more accurate, do it individually: 
        # scattered_E = [compton_scattered_energy(E_in, theta) for each photon]
        # But let's do a single average E_in for simplicity.
        E_in_bs = energies_out[backscatter_flags]  
        scattered_E = compton_scattered_energy(E_in_bs, thetas)

        # Now, decide if the scattered photon escapes or not:
        # If it escapes, the deposit is E_in_bs - scattered_E
        # If it is fully absorbed, deposit is basically E_in_bs
        # For a real approach, you'd do multiple layers or partial absorption.
        escapes = np.random.rand(n_bs) < p_escape

        # Energy deposit for escaping photons
        deposit_escape = E_in_bs - scattered_E
        deposit_escape[deposit_escape < 0] = 0  # no negative

        # For photons that do NOT escape, deposit is ~ all E_in_bs
        deposit_full = E_in_bs  # simplified

        # Combine the two
        deposit_total = np.where(escapes, deposit_escape, deposit_full)

        # Update energies_out
        energies_out[backscatter_flags] = deposit_total
    
    return energies_out


def mini_study_3_geometry_lead(num_photons=100000,
                               angles_deg=[30, 90, 135],
                               FWHM_at_662=0.07,
                               PMT_offset=2.0,
                               angle_uncertainty_deg=1.0,
                               fraction_background=0.1,
                               lead_thickness_cm=1.0,
                               mu_lead=1.2,
                               p_backscatter=0.1,
                               p_escape=0.5):
    """
    Extended geometry refinement function, adding lead pass-through
    AND a naive model of iodine backscatter in the NaI detector.

    :param p_backscatter: Probability that a photon Compton-scatters in the crystal.
    :param p_escape: Probability that the scattered photon escapes the crystal 
                     (thus leaving a partial deposit).
    All other params as before.
    """

    angles_rad = np.radians(angles_deg)
    results = {}
    
    # Probability that a 662 keV photon passes through 'lead_thickness_cm' of lead
    transmission_prob = np.exp(-mu_lead * lead_thickness_cm)

    for angle_deg, angle_rad in zip(angles_deg, angles_rad):
        # --------------------------------------------------------------
        # Part A: Main single-scatter from the target, with angle misalignment
        # --------------------------------------------------------------
        angle_dev = np.random.normal(0, np.radians(angle_uncertainty_deg), num_photons)
        actual_angles = angle_rad + angle_dev
        
        E_scattered_main = compton_scattered_energy(E0, actual_angles)
        
        # --------------------------------------------------------------
        # Part B: Background scattering from environment
        # --------------------------------------------------------------
        n_bg = int(fraction_background * num_photons)
        thetas_bg = sample_theta_klein_nishina(E0, size=n_bg)
        E_scattered_bg = compton_scattered_energy(E0, thetas_bg)
        
        E_combined = np.concatenate((E_scattered_main, E_scattered_bg), axis=0)
        
        # --------------------------------------------------------------
        # Part C: Photons passing through lead (unscattered survivors)
        # --------------------------------------------------------------
        fraction_lead_pool = 0.2
        n_lead_attempt = int(fraction_lead_pool * num_photons)
        n_lead_pass = int(transmission_prob * n_lead_attempt)
        
        thetas_lead = sample_theta_klein_nishina(E0, size=n_lead_pass)
        E_scattered_lead = compton_scattered_energy(E0, thetas_lead)
        
        E_final = np.concatenate((E_combined, E_scattered_lead), axis=0)
        
        # --------------------------------------------------------------
        # Part D: Additional multi-scatter spread (some random broadening)
        # --------------------------------------------------------------
        tiny_spread = 0.01 * E_final
        energies_before_detector = E_final + np.random.uniform(-tiny_spread, tiny_spread)
        
        # --------------------------------------------------------------
        # Part E: Simplified PMT offset (shift) 
        #         (representing detection threshold, window losses, etc.)
        # --------------------------------------------------------------
        energies_measured = energies_before_detector - PMT_offset
        energies_measured[energies_measured < 0] = 0.0
        
        # --------------------------------------------------------------
        # Part F: *NEW* Iodine backscatter inside the detector
        #         (Simplified single-scatter + possible escape)
        # --------------------------------------------------------------
        energies_after_iodine = simulate_iodine_backscatter(
            energies_measured,
            p_backscatter=p_backscatter,
            p_escape=p_escape
        )
        
        # --------------------------------------------------------------
        # Part G: Apply detector resolution (Gaussian smearing)
        # --------------------------------------------------------------
        energies_smeared = gaussian_energy_smearing(energies_after_iodine, FWHM_at_662)
        
        # Save results for each angle
        results[angle_deg] = energies_smeared
    
    # --------------------------------------------------------------
    # Plot all angles
    # --------------------------------------------------------------
    plt.figure(figsize=(9,6))
    for angle_deg in angles_deg:
        energies = results[angle_deg]
        plot_spectrum(energies,
                      bins=80,
                      label=f"{angle_deg}° + Iodine Backscatter")
    
    plt.xlabel("Energy (keV)")
    plt.ylabel("Normalized Counts")
    plt.title("Mini-Study 3 (Extended): Geometry + Lead + Iodine Backscatter")
    plt.legend()
    plt.show()
    
    return results


if __name__ == "__main__":
    # Example usage:
    mini_study_3_geometry_lead(
        num_photons=100000,
        angles_deg=[30],
        FWHM_at_662=0.07,
        PMT_offset=2.0,
        angle_uncertainty_deg=1.0,
        fraction_background=1.0,
        lead_thickness_cm=10.0,
        mu_lead=1.2,
        p_backscatter=0.2,  # try 20% chance of single scatter in the crystal
        p_escape=0.7        # if scattered, 70% chance it escapes
    )

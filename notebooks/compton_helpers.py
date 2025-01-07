
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# 1. Global Constants and Helper Functions
# =============================================================================

# Physical Constants
E0 = 662.0              # keV (Initial photon energy for Cs-137)
m_e = 511.0             # keV (electron rest mass)
r0 = 2.8179403262e-15    # m (classical electron radius) - for completeness
c = 3e8                  # m/s (speed of light) - might be useful if needed

# -----------------------------------------------------------------------------
def compton_scattered_energy(E_in, theta):
    """
    Returns the scattered photon energy (keV) for a photon with initial 
    energy E_in (keV) at scattering angle theta (radians).
    Formula: E' = E0 / [1 + (E0/m_e)*(1 - cos(theta))]
    """
    return E_in / (1.0 + (E_in / m_e) * (1.0 - np.cos(theta)))

# -----------------------------------------------------------------------------
def klein_nishina_pdf(theta, E_in):
    """
    Returns the (unnormalized) Klein-Nishina differential cross section
    for scattering angle theta (radians) and incident photon energy E_in (keV).
    
    dσ/dΩ = r0^2 * (E'/E0)^2 * (E'/E0 + E0/E' - sin^2(theta))
    We omit multiplicative constants for simplicity since we only need
    a *relative* PDF to sample angles.
    """
    E_out = compton_scattered_energy(E_in, theta)
    ratio = E_out / E_in
    return (ratio ** 2) * (ratio + 1/ratio - (np.sin(theta))**2)

# -----------------------------------------------------------------------------
def sample_theta_klein_nishina(E_in, size=1):
    """
    Randomly sample scattering angles from the Klein-Nishina distribution
    for an incident photon energy E_in. Uses an acceptance-rejection method.
    """
    # We'll sample theta in [0, pi].
    # Because Compton scattering is forward-biased for high energies,
    # but let's keep it general.

    # 1. Create an array of random thetas in [0, pi].
    # 2. Accept/reject using the Klein-Nishina PDF as weighting.
    
    thetas = []
    
    # We’ll do a loop until we gather 'size' accepted angles.
    # For large 'size', consider more efficient sampling or vectorization.
    # Here, we do a simple approach for clarity.
    
    max_pdf = None
    
    # Pre-compute a small grid to find the approximate maximum for acceptance-rejection
    test_th = np.linspace(0, np.pi, 500)
    test_pdf = klein_nishina_pdf(test_th, E_in)
    max_pdf = np.max(test_pdf)
    
    count_accepted = 0
    
    while count_accepted < size:
        # Propose a random angle in [0, pi]
        th_proposal = np.random.rand() * np.pi
        
        # Evaluate PDF at that angle
        pdf_val = klein_nishina_pdf(th_proposal, E_in)
        
        # Accept/reject
        if np.random.rand() * max_pdf <= pdf_val:
            thetas.append(th_proposal)
            count_accepted += 1
    
    return np.array(thetas)

# -----------------------------------------------------------------------------
def gaussian_energy_smearing(energy_array, FWHM_at_662=0.07):
    """
    Applies a Gaussian smearing to an array of energies to simulate
    detector resolution.
    
    :param energy_array: Array of photon energies (keV).
    :param FWHM_at_662: Fractional FWHM at 662 keV, e.g. 0.07 means 7% FWHM.
                       Typically FWHM scales ~ sqrt(E), so we can assume
                       sigma(E) = (FWHM_at_662 * E / 2.355) * sqrt(E / 662).
    :return: Smeared array of energies.
    """
    # For a detector, a common approach is: 
    # sigma(E) ∝ sqrt(E). We assume the fraction at 662 keV is known.
    smeared = []
    for E in energy_array:
        # scale the sigma
        sigma_E = (FWHM_at_662 * 662.0 / 2.355) * np.sqrt(E / 662.0)  
        # shift and randomize
        E_smeared = np.random.normal(E, sigma_E)
        if E_smeared < 0:
            E_smeared = 0.0  # no negative energies
        smeared.append(E_smeared)
    return np.array(smeared)

# -----------------------------------------------------------------------------
def plot_spectrum(energies, bins, label, color=None, alpha=0.5):
    """
    Helper to plot a histogram (spectrum) of energies.
    """
    plt.hist(energies, bins=bins, histtype='step', density=True, 
             label=label, color=color, alpha=alpha)

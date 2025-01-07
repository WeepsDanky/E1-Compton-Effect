import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Helper Functions for Plotting Spectra
def plot_spectrum(energy, bins=80, label=None, linestyle='-', color=None):
    """
    Quick histogram plotting for energy data.
    """
    counts, bin_edges = np.histogram(energy, bins=bins, range=(0, 2000))
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    plt.step(bin_centers, counts, label=label, linestyle=linestyle, color=color)
    return counts, bin_centers


# Figure 1: Experimental Setup Schematic
def figure_1_schematic():
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.text(0.1, 0.8, "Cs-137 Source", bbox=dict(facecolor='red', alpha=0.5), fontsize=12)
    ax.text(0.5, 0.8, "Aluminum Rod (Scatterer)", bbox=dict(facecolor='blue', alpha=0.5), fontsize=12)
    ax.text(0.9, 0.8, "Lead Shielding", bbox=dict(facecolor='gray', alpha=0.5), fontsize=12)
    ax.text(0.9, 0.4, "NaI(Tl) Detector", bbox=dict(facecolor='green', alpha=0.5), fontsize=12)
    ax.text(0.9, 0.3, "PMT + MCA", bbox=dict(facecolor='yellow', alpha=0.5), fontsize=12)
    ax.arrow(0.1, 0.8, 0.3, 0, head_width=0.05, head_length=0.05, fc='black', ec='black')
    ax.arrow(0.5, 0.8, 0.3, -0.4, head_width=0.05, head_length=0.05, fc='black', ec='black')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    plt.title("Experimental Setup Schematic")
    plt.show()


# Figure 2: Comparison of Baseline Simulated Spectrum and Measured Data
def figure_2_baseline_vs_measured(experiment_file, simulation_data):
    exp_data = pd.read_csv(experiment_file)
    simulated_energy = simulation_data

    plt.figure(figsize=(10, 6))
    plot_spectrum(exp_data['Energy E / keV'], label='Measured Data at 45°', color='blue')
    plot_spectrum(simulated_energy, label='Baseline Simulation at 45°', color='orange')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend()
    plt.title("Comparison of Baseline Simulated Spectrum and Measured Data")
    plt.show()


# Figure 3: Side-by-Side Spectra Showing Old vs. New Simulation vs. Experiment
def figure_3_old_vs_new_simulation(experiment_file, baseline_sim, refined_sim):
    exp_data = pd.read_csv(experiment_file)
    baseline_energy = baseline_sim
    refined_energy = refined_sim

    plt.figure(figsize=(10, 6))
    plot_spectrum(exp_data['Energy E / keV'], label='Measured Data', color='blue')
    plot_spectrum(baseline_energy, label='Baseline Simulation', color='orange')
    plot_spectrum(refined_energy, label='Refined Simulation', color='green')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend()
    plt.title("Side-by-Side Spectra: Old vs. New Simulation vs. Experiment")
    plt.show()


# Figure 4: Final Comparison of Refined Simulation vs. Experiment
def figure_4_final_comparison(angle_files, refined_sim):
    plt.figure(figsize=(12, 8))
    for angle, file in angle_files.items():
        exp_data = pd.read_csv(file)
        plot_spectrum(exp_data['Energy E / keV'], label=f'Measured Data {angle}°')
    plot_spectrum(refined_sim, label='Refined Simulation', linestyle='--', color='black')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend()
    plt.title("Final Comparison of Refined Simulation vs. Experiment")
    plt.show()


# Figure 5: Error Analysis Breakdown for Each Mini-Study
def figure_5_error_breakdown():
    studies = ["Baseline", "Detector Refinement", "Geometry Refinement"]
    stat_errors = [10, 7, 4]  # Example values
    sys_errors = [5, 3, 2]

    x = np.arange(len(studies))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width/2, stat_errors, width, label='Statistical Errors')
    ax.bar(x + width/2, sys_errors, width, label='Systematic Errors')

    ax.set_xlabel('Mini-Studies')
    ax.set_ylabel('Errors (%)')
    ax.set_title('Error Analysis Breakdown for Each Mini-Study')
    ax.set_xticks(x)
    ax.set_xticklabels(studies)
    ax.legend()
    plt.show()


# Figure 6: Impact of Energy-Dependent Detector Resolution
def figure_6_resolution_impact(baseline_sim, refined_sim):
    plt.figure(figsize=(10, 6))
    plot_spectrum(baseline_sim, label='Baseline (Uniform Resolution)', color='orange')
    plot_spectrum(refined_sim, label='Refined (Energy-Dependent Resolution)', color='green')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend()
    plt.title("Impact of Energy-Dependent Detector Resolution")
    plt.show()


# Figure 7: Geometry Refinements and Their Impact
def figure_7_geometry_impact(sim_no_geometry, sim_with_geometry):
    plt.figure(figsize=(10, 6))
    plot_spectrum(sim_no_geometry, label='Simulation Without Geometry Refinement', color='red')
    plot_spectrum(sim_with_geometry, label='Simulation With Geometry Refinement', color='green')
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend()
    plt.title("Geometry Refinements and Their Impact")
    plt.show()


# Main Execution: Generate All Figures
if __name__ == "__main__":
    # Example file paths and synthetic data for testing
    example_experiment_file = "data/experiment/1128/cs137-sspht022-2048bin-100s-2gain-45degree-no-shielding.csv"
    example_simulation_baseline = np.random.normal(500, 100, 10000)
    example_simulation_refined = np.random.normal(480, 50, 10000)

    # File paths for multiple angles
    angle_files = {
        30: "data/experiment/1129/cs137-sspht022-2048bin-600s-2gain-30degree-no-al.csv",
        60: "data/experiment/1129/cs137-sspht022-2048bin-600s-2gain-60degree-no-al.csv",
        90: "data/experiment/1129/cs137-sspht022-2048bin-600s-2gain-90degree-no-al.csv",
        120: "data/experiment/1129/cs137-sspht022-2048bin-600s-2gain-120degree-no-al.csv",
        135: "data/experiment/1129/cs137-sspht022-2048bin-600s-2gain-135degree-no-al.csv"
    }

    # Generate each figure
    figure_1_schematic()
    figure_2_baseline_vs_measured(example_experiment_file, example_simulation_baseline)
    figure_3_old_vs_new_simulation(example_experiment_file, example_simulation_baseline, example_simulation_refined)
    figure_4_final_comparison(angle_files, example_simulation_refined)
    figure_5_error_breakdown()
    figure_6_resolution_impact(example_simulation_baseline, example_simulation_refined)
    figure_7_geometry_impact(example_simulation_baseline, example_simulation_refined)

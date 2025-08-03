import numpy as np
import matplotlib.pyplot as plt

def parse_fortran_output(filename):
    """Parse Fortran output file into time, current, and voltage arrays."""
    times, currents, voltages = [], [], []
    
    with open(filename, 'r') as file:
        for line in file:
            if 'Time(s)=' in line:
                parts = line.strip().split()
                time = float(parts[1])
                current = float(parts[3])
                voltage = float(parts[5])
                times.append(time)
                currents.append(current)
                voltages.append(voltage)
    
    return np.array(times), np.array(currents), np.array(voltages)

def plot_rlc_response(times, currents, voltages):
    """Plot current and voltage vs. time."""
    plt.figure(figsize=(12, 6))
    
    # Plot current (I)
    plt.subplot(2, 1, 1)
    plt.plot(times, currents, 'b-', linewidth=1.5, label='Current (I)')
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.title('RLC Circuit Response')
    plt.grid(True)
    plt.legend()
    
    # Plot capacitor voltage (V_C)
    plt.subplot(2, 1, 2)
    plt.plot(times, voltages, 'r-', linewidth=1.5, label='Voltage (V_C)')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.grid(True)
    plt.legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # --- Parse Fortran output ---
    output_file = "rlc_output.dat"  # Replace with your output file
    times, currents, voltages = parse_fortran_output(output_file)
    
    # --- Plot results ---
    plot_rlc_response(times, currents, voltages)
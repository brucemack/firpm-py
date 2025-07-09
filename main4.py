import numpy as np
import math 
from scipy.signal import lfilter, firwin
import matplotlib.pyplot as plt

def fir_freq_response(coefficients, sampling_rate):
    """
    Calculates the frequency response of an FIR filter.

    Args:
        coefficients: A list or NumPy array of FIR filter coefficients.
        sampling_rate: The sampling rate of the signal.

    Returns:
        A tuple containing:
            - frequencies: A NumPy array of frequencies (Hz).
            - frequency_response: A NumPy array of complex frequency response values.
    """
    # Calculate the frequency response using the DFT (FFT)
    frequency_response = np.fft.fft(coefficients)
    
    # Number of coefficients
    num_coeffs = len(coefficients)
    
    # Generate the corresponding frequency values
    frequencies = np.fft.fftfreq(num_coeffs, 1 / sampling_rate)
    
    return frequencies, frequency_response


# Coefficients for a half-band filter
def half_band_lpf(N: int):
    # K is the number of 1's in the rectangle, which is 
    # N/2 by definition of the half-band filter.
    K = float(N) / 2.0
    result = []
    for n in range(-int(N/2), int(N/2) + 1):
        if n == 0:
            h = (1.0 / float(N)) * K
        else:
            h = (1.0 / float(N)) * math.sin(np.pi * n * K / float(N)) / math.sin(np.pi * n / float(N))
        # Hamming window function
        w = 0.54 + 0.46 * math.cos(2.0 * np.pi * n / N)
        b = h * w
        if abs(b) < 1e-16:
            b = 0
        result.append(b)
    return result

fs = 16000
h = half_band_lpf(19)
print(h)

# Calculate the frequency response
frequencies, frequency_response = fir_freq_response(h, fs)
db_data = 20 * np.log10(np.abs(frequency_response[:len(frequencies)//2]))

fig, ax = plt.subplots()
ax.plot(frequencies[:len(frequencies)//2], db_data)
ax.set_xscale('log') 
#ax.set_yscale('log') 
ax.set_title('Magnitude Response')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Magnitude (dB)')
ax.grid(True)
plt.show()

# Simulate a signal being filtered
t = np.linspace(0, int(fs / 100), int(fs / 100), endpoint=False)

ft0 = 134
a0 = 1
omega0 = 2.0 * 3.1415926 * ft0 / fs
signal0 = a0 * np.sin(omega0 * t)

# Apply the LP filter
signal1 = lfilter(h, [1.0], signal0)

fig, ax = plt.subplots()
ax.plot(t, signal1)
ax.set_title('Filtered Signal')
ax.set_xlabel('Time')
ax.set_ylabel('Magnitude')
ax.grid(True)
plt.show()

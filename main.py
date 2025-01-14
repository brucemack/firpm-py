import math
import firpm 
import numpy as np
import matplotlib.pyplot as plt

print("Figure 7 - N=24 low-pass filter")
h, dev = firpm.design(24, 1, 2, [ 0.00, 0.08, 0.16, 0.5 ], [ 1.0, 0.0 ], [ 1.0, 1.0 ])
print(h)
print("Dev dB", 20 * math.log10(dev))

print("Figure 9 - N=32 bandpass filter")
h, dev = firpm.design(32, 1, 3, [ 0.00, 0.10, 0.20, 0.35, 0.425, 0.5 ], 
            [ 0.0, 1.0, 0.0 ], [ 10.0, 1.0, 10.0 ])
print(h)
print("Dev", dev)
print("Dev dB", 20 * math.log10(dev))

print("Figure 11 - N=50 bandpass filter with unequal weighting in the stopbands")
h, dev = firpm.design(50, 1, 3, [ 0.00, 0.15, 0.20, 0.30, 0.35, 0.5 ], 
            [ 0.0, 1.0, 0.0 ], [ 10.0, 1.0, 100.0 ])
print(h)
print("Dev", dev)
print("Dev dB", 20 * math.log10(dev))

print("Figure 13 - N=31 bandstop filter")
h, dev = firpm.design(31, 1, 3, [ 0.00, 0.1, 0.15, 0.35, 0.42, 0.5 ], 
            [ 1.0, 0.0, 1.0 ], [ 1.0, 50.0, 1.0 ])
print(h)
print("Dev", dev)
print("Dev dB", 20 * math.log10(dev))
# Do an IFFT of this one
x = np.fft.ifft(np.array(h))
# Plot "left half"
plt.plot(np.abs(x[0 : round(len(x) / 2)]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.title('Parks-McClellan Filter')
plt.grid(True)
plt.show()

print("Figure 19 - N=20 Hilbert transformer")
h, dev = firpm.design(20, 3, 1, [ 0.05, 0.5 ], [ 1.0 ], [ 1.0 ])
print(h)
print("Dev dB", 20 * math.log10(dev))
# Do an IFFT of this one
y = np.fft.ifft(np.array(h))
points = round(len(y) / 2)
# Make an x-axis
x = []
for i in range(0, points):
    x.append(0.5 * (i / points))
# Plot "left half"
plt.plot(x, np.abs(y[0 : points]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude')
plt.title('Parks-McClellan Filter Design - Hilbert Transformer')
plt.grid(True)
plt.show()

# Plot "left half"
plt.plot(x, np.angle(y[0 : points]))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase')
plt.title('Parks-McClellan Filter Design - Hilbert Transformer')
plt.grid(True)
plt.show()

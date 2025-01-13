import math
import firpm 

print("Figure 7 - N=24 low-pass filter")
h, dev = firpm.design(24, 1, 2, [ 0.00, 0.08, 0.16, 0.5 ], [ 1.0, 0.0 ], [ 1.0, 1.0 ])
h.dump()
print("Dev dB", 20 * math.log10(dev))

print("Figure 9 - N=32 bandpass filter")
h, dev = firpm.design(32, 1, 3, [ 0.00, 0.10, 0.20, 0.35, 0.425, 0.5 ], 
            [ 0.0, 1.0, 0.0 ], [ 10.0, 1.0, 10.0 ])
h.dump()
print("Dev", dev)
print("Dev dB", 20 * math.log10(dev))

print("Figure 11 - N=50 bandpass filter with unequal weighting in the stopbands")
h, dev = firpm.design(50, 1, 3, [ 0.00, 0.15, 0.20, 0.30, 0.35, 0.5 ], 
            [ 0.0, 1.0, 0.0 ], [ 10.0, 1.0, 100.0 ])
h.dump()
print("Dev", dev)
print("Dev dB", 20 * math.log10(dev))

print("Figure 19 - N=20 Hilbert transformer")
h, dev = firpm.design(20, 3, 1, [ 0.05, 0.5 ], [ 1.0 ], [ 1.0 ])
h.dump()
print("Dev dB", 20 * math.log10(dev))

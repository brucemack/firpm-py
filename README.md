Overview
========

A literal translation of the famous Parks-McClellan algorithm for FIR filter
design from the original FORTRAN to pure Python. This was done for learning 
purposes and to provide a basis for further implementation work in Python 
and other languages.

Validation
==========

The results match the examples provided in the original paper with a high
degree of accuracy. Please see [validate.py](validate.py) for unit tests that demonstrate
matching results covering 
all of the examples in the original paper _A Computer Program for Designing Optimum FIR Linear Phase Digital Filters_. No output 
means that the tests ran successfully.

Usage Notes
===========

The design function in the firpm module computes the coefficients
for the realization of the optimal filter. Here we describe the 
parameters in order:

*nfilt* is the length of the filter, also known as the number of taps.

*jtype* is the type of filter.  Use 1 for a multiple passband/stopband filter, 2 for a differentiator, 3 for a Hilbert transform filter.

*nbands* is the number of bands in the description of the optimal filter.
A "band" is essentially a range of frequencies for which a gain 
constraint is specified. Bands must be specified in order and they 
need not be contiguous. Bands should not be overlapping. A maximum of
10 bands is allowed.

*edges* is an array that contains the pairs of frequencies that 
define the start and end 
frequency of each band. The format is [ start0, end0, start1, end1, ... ].
Frequencies are normalized to the sampling rate of the system. So 
the minimum frequency is 0 and the maximum frequency is 0.5.

*gains* is an array that contains the desired response per band.  One 
entry per band. Each entry is a number from 0 to 1.

*weights* is an array that contains the weight per band.  One entry 
per band. These coefficients are used to allow the designer to 
specify the relative magnitude of the error in each band.  

*lgrid* is the grid density used for estimation, which defaults to 16.

References
==========

[A PDF copy of the original paper](https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/062_computer%20program.pdf) (not sure if this is a legal copy since it's an IEEE published 
paper).

[James McClellan's original Master's thesis from Rice Uiversity is here](https://repository.rice.edu/server/api/core/bitstreams/a924e584-8512-4852-9801-c602985dc0da/content)

[A helpful reference to the FORTRAN code](https://michaelgellis.tripod.com/dsp/pgm21.html)


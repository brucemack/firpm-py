Overview
========

A literal translation of the famous Parks-McClellan algorithm for FIR filter
design from the original FORTRAN to pure Python. This was done for learning 
purposes and to provide a basis for further implementation work in Python 
and other languages.

Validation
==========

The results match the examples provided in the original paper with a high
degree of accuracy. Please see [validate.py](blob/main/validate.py) for unit tests that demonstrate
matching results covering 
all of the examples in the original paper _A Computer Program for Designing Optimum FIR Linear Phase Digital Filters_.

References
==========

[A PDF copy of the original paper](https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/062_computer%20program.pdf) (not sure if this is a legal copy since it's an IEEE published 
paper).

[James McClellan's original Master's thesis from Rice Uiversity is here](https://repository.rice.edu/server/api/core/bitstreams/a924e584-8512-4852-9801-c602985dc0da/content)

[A helpful reference to the FORTRAN code](https://michaelgellis.tripod.com/dsp/pgm21.html)


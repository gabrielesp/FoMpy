##FoMPy Instructions

This project extracts several useful figure of merits (FoM) from a semiconductor's IV curve:

- [x] Threshold voltage (V<sub>T</sub>)
- [x] OFF-current (I<sub>OFF</sub>)
- [x] ON-current (I<sub>ON</sub>)
- [x] Sub-threshold slope (SS)
- [x] Drain induced barrier lowering (DIBL)
- [x] Static power (P<sub>static</sub>)

## QUICK START

To run the **FoMPy** do the following:

$ python3 example.py

So far the only external librariy used other than the ones in ./src 
is numpy. If you do not have numpy installed in your system, and this is a 
unix-based execute the following command:

$ pip3 install numpy

## Code Examples

Inside example.py several code examples have been implemented to check if the code runs properly.

## References

[1] A. Ortiz-Conde, F.J. García Sánchez, J.J. Liou, A. Cerdeira, M. Estrada, Y. Yue, *A review of recent MOSFET threshold voltage extraction methods*,
In Microelectronics Reliability, Volume 42, Issues 4–5, 2002, Pages 583-596, ISSN 0026-2714,
[https://doi.org/10.1016/S0026-2714(02)00027-6](https://doi.org/10.1016/S0026-2714(02)00027-6).
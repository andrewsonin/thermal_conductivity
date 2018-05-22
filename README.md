# thermal_conductivity

Andrew Sonin, 5113. Variant 4.

Two-layer six-point difference scheme for solving the heat equation.

u'_t = u''_{xx} + e^t * sin(x), x \in [-\pi/2, \pi/2], t \in [0, T]

u(t=0) = sin(x), u(x=-\pi/2) = -cosh(t) = -u(x=\pi/2)

exact solution: u = cosh(t) * sin(x)

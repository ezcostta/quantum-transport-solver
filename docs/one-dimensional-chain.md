# One-Dimensional Tight-Binding Scattering Solver

## 1. Physical model

The first implemented model is a one-dimensional tight-binding chain connected to two identical semi-infinite leads.

The Hamiltonian convention is

[
H
=

## \sum_j \epsilon_j c_j^\dagger c_j

t\sum*j
\left(
c*{j+1}^\dagger c_j

- c*j^\dagger c*{j+1}
  \right).
  ]

The corresponding stationary Schrödinger equation is

[
E\psi_j
=======

## \epsilon_j\psi_j

## t\psi\_{j-1}

t\psi\_{j+1}.
]

Rearranging,

[
(E-\epsilon_j)\psi_j

- t\psi\_{j-1}
- # t\psi\_{j+1}

0.

]

This sign convention is important. It was also the source of the first bug found by the physics tests.

---

## 2. Lead dispersion

In the left and right leads, the onsite energy is taken as

[
\epsilon=0.
]

For a plane wave

[
\psi_j = e^{ikj},
]

the tight-binding equation gives

[
E = -2t\cos k.
]

More generally, if the lead onsite energy is (\epsilon_0),

[
E = \epsilon_0 - 2t\cos k.
]

Therefore,

[
\cos k = \frac{\epsilon_0-E}{2t}.
]

The code implements this in `OneDimensionalLead.wave_number`.

Only energies satisfying

[
\left|\frac{\epsilon_0-E}{2t}\right| \leq 1
]

belong to the propagating band. Energies outside this interval raise a `ValueError`.

---

## 3. Scattering ansatz

The finite scattering region contains (N) sites:

[
j=0,1,\dots,N-1.
]

The left lead occupies

[
j\leq -1,
]

and the right lead occupies

[
j\geq N.
]

For an incoming wave from the left, the wavefunction is written as

[
\psi_j
======

e^{ikj}

- r e^{-ikj},
  \qquad j\leq -1,
  ]

where (r) is the reflection amplitude.

In the right lead,

[
\psi_j
======

\tau e^{ikj},
\qquad j\geq N,
]

where (\tau) is the transmission amplitude.

The unknowns of the problem are therefore

[
\psi_0,\psi_1,\dots,\psi_{N-1},r,\tau.
]

The solver constructs a linear system

[
A x = b,
]

where

[
x =
(\psi_0,\psi_1,\dots,\psi_{N-1},r,\tau)^T.
]

---

## 4. Equations inside the scattering region

For each internal site (j), the equation is

[
(E-\epsilon_j)\psi_j

- t\psi\_{j-1}
- # t\psi\_{j+1}

0.

]

For (1\leq j\leq N-2), both neighbors are inside the scattering region.

For (j=0), the left neighbor is in the left lead:

[
\psi\_{-1}
=========

e^{-ik}

- r e^{ik}.
  ]

For (j=N-1), the right neighbor is in the right lead:

[
\psi_N
======

\tau e^{ikN}.
]

These boundary substitutions are what couple the internal scattering-region amplitudes to the reflection and transmission amplitudes.

---

## 5. Left boundary equation

The solver also imposes the tight-binding equation at the last site of the left lead, (j=-1):

[
E\psi\_{-1}
==========

## -t\psi\_{-2}

t\psi_0.
]

Using

[
\psi\_{-1}
=========

e^{-ik}

- r e^{ik},
  ]

and

[
\psi\_{-2}
=========

e^{-2ik}

- r e^{2ik},
  ]

one obtains an equation involving (\psi_0) and (r).

This equation is the first row of the matrix system.

---

## 6. Right boundary equation

Similarly, the solver imposes the tight-binding equation at the first site of the right lead, (j=N):

[
E\psi_N
=======

## -t\psi\_{N-1}

t\psi\_{N+1}.
]

Using

[
\psi_N
======

\tau e^{ikN},
]

and

[
\psi\_{N+1}
==========

\tau e^{ik(N+1)},
]

one obtains an equation involving (\psi\_{N-1}) and (\tau).

This equation is the last row of the matrix system.

---

## 7. Structure of the linear system

The unknown vector is ordered as

[
x =
(\psi_0,\psi_1,\dots,\psi_{N-1},r,\tau)^T.
]

Therefore:

- columns (0,\dots,N-1) correspond to the scattering-region wavefunction;
- column (N) corresponds to (r);
- column (N+1) corresponds to (\tau).

The matrix has size

[
(N+2)\times(N+2),
]

because there are:

- (N) equations for the scattering region;
- one equation at the left boundary;
- one equation at the right boundary.

The system is solved with

[
x = A^{-1}b,
]

implemented numerically using

```python
np.linalg.solve(A, b)
```

rather than explicitly computing the inverse.

---

## 8. Reflection and transmission probabilities

For identical left and right leads, the group velocities are equal.

Therefore,

[
R = |r|^2,
]

and

[
T = |\tau|^2.
]

The probability conservation check is

[
R+T=1.
]

This is implemented in `observables.py`.

For non-identical leads, the transmission probability will later need a velocity factor:

[
T =
\frac{v_R}{v_L}
|\tau|^2.
]

This will become important when we generalize to matrix-valued leads, graphene, and Bogoliubov-de Gennes systems.

---

## 9. Physics tests

The current test suite checks four basic physical properties.

### Probability conservation

For a Hermitian scattering region connected to identical leads,

[
R+T=1.
]

### Uniform chain

If the scattering region is identical to the leads,

[
\epsilon_j=0,
]

there should be no reflection:

[
R=0,
\qquad
T=1.
]

This test caught a sign error in the boundary equations.

### Barrier

For a nonzero onsite potential, the system should scatter:

[
0<T<1.
]

### Energy outside the band

If

[
|E|>2t,
]

there are no propagating states in the lead, so the solver raises an exception.

---

## 10. Why this model matters

Although this is only a scalar one-dimensional chain, it already contains the complete logic of a scattering calculation:

1. define lead modes;
2. write an incoming plus reflected wave on the left;
3. write a transmitted wave on the right;
4. impose boundary matching;
5. assemble a linear system;
6. solve for scattering amplitudes;
7. compute observables;
8. verify conservation laws.

The future graphene and BdG solvers will follow the same structure, but with matrix-valued wavefunctions and block Hamiltonians.

In that sense, this model is not a toy that will be discarded. It is the minimal working prototype of the full quantum-transport architecture.

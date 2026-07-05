# One-Dimensional Tight-Binding Scattering Solver

## 1. Introduction

The first implemented model in **QTransport** is a one-dimensional tight-binding chain connected to two identical semi-infinite leads.

Although simple, this model already contains the complete workflow of a quantum scattering calculation:

1. Define the lead eigenmodes.
2. Construct the scattering wavefunction.
3. Apply the boundary conditions.
4. Assemble the linear system.
5. Solve for the scattering amplitudes.
6. Compute the reflection and transmission probabilities.
7. Verify the physical conservation laws.

The future graphene, bilayer graphene, and Bogoliubov-de Gennes (BdG) implementations will follow exactly the same workflow, differing only in the internal structure of the wavefunction and Hamiltonian matrices.

---

# 2. Physical Model

The Hamiltonian of the one-dimensional chain is

$$
H
=

\sum_j
\epsilon_j
c_j^\dagger c_j
---------------

t
\sum_j
\left(
c_{j+1}^\dagger c_j
+
c_j^\dagger c_{j+1}
\right),
$$

where

- $t$ is the hopping energy;
- $\epsilon_j$ is the onsite potential.

The stationary Schrödinger equation becomes

$$
E\psi_j
=======

## \epsilon_j\psi_j

## t\psi_{j-1}

t\psi_{j+1}.
$$

Rearranging,

$$
(E-\epsilon_j)\psi_j
+
t\psi_{j-1}
+
t\psi_{j+1}
===========

0.


$$

This equation is the one implemented by the solver.

---

# 3. Geometry

The system consists of

```text
Left Lead        Scattering Region          Right Lead

←∞ ... -2 -1 | 0 1 2 ... N−1 | N N+1 ... +∞ →
```

The scattering region contains

$$
N
$$

sites,

$$
j=0,\ldots,N-1.
$$

The leads are semi-infinite.

---

# 4. Lead Dispersion

Inside the leads,

$$
\epsilon_j=\epsilon_0.
$$

Assuming a plane wave

$$
\psi_j=e^{ikj},
$$

the dispersion relation is

$$
E
=

## \epsilon_0

2t\cos k.
$$

Therefore,

$$
\cos k
======

\frac{\epsilon_0-E}{2t}.
$$

The method

```python
OneDimensionalLead.wave_number()
```

computes

$$
k
=

\arccos
\left(
\frac{\epsilon_0-E}{2t}
\right).
$$

Only energies satisfying

$$
\left|
\frac{\epsilon_0-E}{2t}
\right|
\le1
$$

correspond to propagating states.

Outside this interval the solver raises

```python
ValueError
```

because the lead has no propagating modes.

---

# 5. Scattering Ansatz

An incoming electron comes from the left lead.

The wavefunction is written as

## Left lead

$$
\psi_j
======

e^{ikj}
+
r
e^{-ikj},
\qquad
j\le -1,
$$

where

- the first term is the incoming wave;
- $r$ is the reflection amplitude.

---

## Right lead

Only transmitted waves are allowed,

$$
\psi_j
======

\tau
e^{ikj},
\qquad
j\ge N,
$$

where

$$
\tau
$$

is the transmission amplitude.

---

# 6. Unknown Variables

The solver treats the following quantities as unknowns:

$$
\psi_0,
\psi_1,
\ldots,
\psi_{N-1},
r,
\tau.
$$

Therefore the unknown vector is

$$
x=
\begin{pmatrix}
\psi_0\
\psi_1\
\vdots\
\psi_{N-1}\
r\
\tau
\end{pmatrix}.
$$

The matrix system has

$$
N+2
$$

unknowns.

---

# 7. Equations Inside the Scattering Region

For every site inside the scattering region,

$$
(E-\epsilon_j)\psi_j
+
t\psi_{j-1}
+
t\psi_{j+1}
===========

0.


$$

For interior sites

$$
1\le j\le N-2,
$$

both neighbors belong to the scattering region.

For the edge sites,

- the left neighbor belongs to the left lead;
- the right neighbor belongs to the right lead.

These substitutions introduce the unknown reflection and transmission amplitudes into the system.

---

# 8. Left Boundary Equation

The Schrödinger equation is also imposed at

$$
j=-1.
$$

Using

$$
\psi_{-1}
=========

e^{-ik}
+
r
e^{ik},
$$

and

$$
\psi_{-2}
=========

e^{-2ik}
+
r
e^{2ik},
$$

the equation

$$
E\psi_{-1}
==========

## -t\psi_{-2}

t\psi_0
$$

becomes an equation involving only

- $\psi_0$;
- $r$.

This is the first row of the linear system.

---

# 9. Right Boundary Equation

Similarly,

$$
j=N
$$

satisfies

$$
E\psi_N
=======

## -t\psi_{N-1}

t\psi_{N+1}.
$$

Using

$$
\psi_N
======

\tau
e^{ikN},
$$

and

$$
\psi_{N+1}
==========

\tau
e^{ik(N+1)},
$$

one obtains an equation involving

- $\psi_{N-1}$;
- $\tau$.

This becomes the last row of the matrix system.

---

# 10. Matrix Structure

The solver assembles

$$
A,x=b.
$$

The unknown ordering is

|   Column | Unknown      |
| -------: | ------------ |
|        0 | $\psi_0$     |
|        1 | $\psi_1$     |
| $\vdots$ | $\vdots$     |
|    $N-1$ | $\psi_{N-1}$ |
|      $N$ | $r$          |
|    $N+1$ | $\tau$       |

The matrix size is

$$
(N+2)\times(N+2).
$$

The solution is obtained with

```python
x = np.linalg.solve(A, b)
```

instead of explicitly computing

$$
A^{-1}.
$$

This is numerically more stable and efficient.

---

# 11. Reflection and Transmission

For identical left and right leads,

$$
R
=

|r|^2,
$$

and

$$
T
=

|\tau|^2.
$$

The solver checks

$$
R+T=1,
$$

which expresses probability conservation.

For different leads this will later become

$$
T
=

\frac{v_R}{v_L}
|\tau|^2,
$$

where

- $v_L$ is the group velocity in the left lead;
- $v_R$ is the group velocity in the right lead.

This generalization will be essential for graphene and BdG systems.

---

# 12. Physics Tests

The current implementation contains four physics tests.

## Probability Conservation

$$
R+T=1.
$$

This verifies that the scattering matrix is unitary.

---

## Uniform Chain

If

$$
\epsilon_j=0,
$$

the scattering region is identical to the leads.

Therefore,

$$
R=0,
\qquad
T=1.
$$

This test detected the first implementation bug (an incorrect sign in the boundary equations).

---

## Potential Barrier

For

$$
\epsilon_j\ne0,
$$

the transmission must satisfy

$$
0<T<1.
$$

---

## Energy Outside the Band

For

$$
|E-\epsilon_0|>2t,
$$

the leads possess no propagating states.

The solver therefore raises

```python
ValueError
```

instead of attempting an unphysical calculation.

---

# 13. Code Organization

The implementation is currently divided into three modules.

## `leads.py`

Defines the semi-infinite lead.

Responsibilities:

- compute the dispersion relation;
- compute the wave number;
- compute the group velocity.

---

## `solver.py`

Builds the complete linear system

$$
A,x=b
$$

and solves it.

Responsibilities:

- assemble the boundary equations;
- assemble the scattering-region equations;
- solve the linear system;
- return the scattering amplitudes.

---

## `observables.py`

Computes physical quantities from the amplitudes.

Responsibilities:

- reflection probability;
- transmission probability;
- probability conservation.

---

# 14. Why This Model Matters

Although the current implementation describes only a scalar one-dimensional chain, it already contains the complete logical structure of the future transport framework.

The same algorithmic steps will be used for

- graphene;
- bilayer graphene;
- ferromagnetic systems;
- superconducting systems;
- Bogoliubov-de Gennes Hamiltonians.

The only difference will be that scalar quantities become vectors and matrices.

For this reason, the one-dimensional chain should not be viewed as a temporary toy model. It is the minimal working prototype upon which the entire **QTransport** framework will be built.

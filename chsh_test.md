---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.15.2
kernelspec:
  display_name: qiskit
  language: python
  name: qiskit
---

---
The CHSH test
---

Bell's theorem asserts that quantum mechanics is incompatible with local
hidden variable theories. Certain constraints on the outcome of
measurements implied by hidden variable theories give rise to so-called
Bell inequalities. One such inequality is the CHSH inequality, named
after the authors Clauser, Horne, Shimony, and Holt, is used to
experimentally prove Bell's theorem.

Suppose Alice and Bob are each given a particle and they can perform on
their respective particle one of two measurements: $A_{0}$ or $A_{1}$
for Alice, and $B_{0}$ or $B_{1}$ for Bob. Moreover, Alice and Bob are
sufficiently far apart so that their actions cannot influence on another
(locality hypothesis). These measurements correspond to certain (hidden)
properties of the particles, $a_{0}$ for $A_{0}$, $a_{1}$ for $A_{1}$
and so forth, so that measurements of these properties always result in
$\pm1$. Now consider the combination

$$
C=a_{0}b_{0}+a_{0}b_{1}+a_{1}b_{0}-a_{1}b_{1}=\left(a_{0}+a_{1}\right)b_{0}+\left(a_{0}-a_{1}\right)b_{1}\,.
$$(combination)

Note that if $a_{0}$ and $a_{1}$ have equal values, for instance,
$a_{0}=a_{1}=1$, then the second term on the right is zero, and the
combination $C$ is $\pm2$. The same happens if $a_{0}$ and $a_{1}$ have
opposite values, because then the first term on the right is zero, and
the combination $C$ is again $\pm2$.

Now each time Alice and Bob receive their particle, they perform an
experiment, which is one of the measurements we described. After
repeating this experiment many times, the average value of the
combination {eq}`combination`  will be less than $2$, and it is the sum of
the average values of each term,

$$
\left\langle C\right\rangle =\left\langle a_{0}b_{0}\right\rangle +\left\langle a_{0}b_{1}\right\rangle +\left\langle a_{1}b_{0}\right\rangle -\left\langle a_{1}b_{1}\right\rangle \leq2\label{eq:CHSH}
$$

This is the CHSH inequality.

In the post "CHSH game", we have shown how quantum mechanics can violate
the CHSH inequality. Let us briefly state the main arguments. This time
Alice and Bob are given each a qubit which belong to an entangled pair

$$
\psi=\frac{1}{\sqrt{2}}\left(\left|00\right\rangle +\left|11\right\rangle \right)
$$(bell-state)

And they can perform on their qubit one of two operations, as before.
These operations correspond to unitary hermitian operations acting on a
qubit. The analog of the combination {eq}`combination` is the two-qubit operator

$$
T=A_{0}\otimes B_{0}+A_{0}\otimes B_{1}+A_{1}\otimes B_{0}-A_{1}\otimes B_{1}
$$

We showed that the average value of this operator with respect to the
Bell state {eq}`bell-state`  is

$$
\left\langle T\right\rangle _{\psi}\leq2\sqrt{2}\,,
$$ 

which can exceed the classical upper bound of $2$, deduced from the local hidden
variables hypothesis.

+++

Now we will test the CHSH inequality on a real quantum computer, using the Estimator primitive. 

```{code-cell} ipython3
from qiskit import QuantumCircuit
from qiskit.extensions import UnitaryGate
from numpy import pi, cos, sin, sqrt
from qiskit.quantum_info.operators import Operator
# Runtime imports
from qiskit_ibm_runtime import QiskitRuntimeService, Estimator, Session, Options
```

```{code-cell} ipython3
#First we write down the hermitian unitary matrices corresponding to the operations Alice and Bob can perform

A0_matrix = [ [cos(pi/8), sin(pi/8)],[sin(pi/8),-cos(pi/8)]]
A1_matrix = [ [cos(pi/8), -sin(pi/8)],[-sin(pi/8),-cos(pi/8)]]
B0_matrix = [ [1, 0],[0,-1]]
B1_matrix = [ [cos(pi/4), sin(pi/4)],[sin(pi/4),-cos(pi/4)]]
```

```{code-cell} ipython3
#now we convert the matrices to the Operator class and take the tensor products

A0 = Operator(A0_matrix)
A1 = Operator(A1_matrix)
B0 = Operator(B0_matrix)
B1 = Operator(B1_matrix)

A0B0 = A0.tensor(B0)
A0B1 = A0.tensor(B1)
A1B0 = A1.tensor(B0)
A1B1 = A1.tensor(B1)

ops = [A0B0,A0B1,A1B0,A1B1]

```

```{code-cell} ipython3
#Next we implement the circuit corresponding the the Bell state

bellqc = QuantumCircuit(2)
bellqc.h(0)
bellqc.cx(0,1)
#for each operator in ops, we apply the estimator, so we need the same number of Bell circuits:
circuits= [bellqc]*4
```

```{code-cell} ipython3
# Now we choose our backend and start our session
service = QiskitRuntimeService(channel='ibm_quantum', instance="ibm-q-research-2/fed-uni-ufabc-1/main", token="b991d7c5d1d6efcf2ec598451d82108070c330d72f3b033846e4aaecff48ae1c65575b7463772cf57e1fc8591310f30e299ed962d96f4e6de74e33d44f4db836")

# Set options, which can be overwritten at job level.
#options = Options(optimization_level=1,resilience_level=1)
options = Options()
options.resilience_level=1
options.optimization_level=1

# Select the system with the fewest number of jobs in the queue
#backend = service.least_busy(simulator=False, operational=True)
backend = service.backend("ibm_lagos")
# Initialize your session
session = Session(backend=backend)
#backend.name
```

```{code-cell} ipython3
#Now we create our Estimator instance and run our circuits with the observables in ops
estimator = Estimator(session=session,options=options)

#calculate [<A0B0>,<A0B1>,<A1B0>,<A1B1>]

job = estimator.run(circuits, ops, shots=int(1e4))
print(f">>> Job ID: {job.job_id()}")
print(f">>> Job Status: {job.status()}")  
```

```{code-cell} ipython3
#now we retrieve our results and calculate the average of the CHSH operator
result = service.job("cjmbo7cvcjlre5d4oarg").result().values
#CHSH average
CHSH = result[0] + result[1] + result[2] - result[3]
```

```{code-cell} ipython3
#Compare with the classical upperbound 2

print("The experiment violates the CHSH inequality by", CHSH-2)
```

```{code-cell} ipython3
#let us compare this result with that of a simulator
from qiskit.primitives import Estimator
estimator = Estimator()
job = estimator.run(circuits, ops, shots=int(1e4))
result = job.result().values
#CHSH average
CHSH2 = result[0] + result[1] + result[2] - result[3]
print("The experiment violates the CHSH inequality by", CHSH2-2)
```

```{code-cell} ipython3

```

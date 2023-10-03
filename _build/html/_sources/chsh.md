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
The CHSH game
---

The CHSH is a nonlocal game where two players (Alice and Bob) cooperate
to win. They are allowed to communicate before the game starts, but once
the game starts, they are not allowed to communicate anymore. The game
itself is defined by a simple rule: Alice and Bob are each provided a
classical bit, $0$ or $1$, and they must respond with another classical
bit. Let's say Alice is provided the classical bit $x\in\{0,1\}$ and her
answer is $a\in\{0,1\}$. Likewise, Bob is provided the classical bit
$y\in\{0,1\}$ and his answer is $b\in\{0,1\}$. They win the CHSH game if
the following rule is satisfied:

$$
x\wedge y=a\oplus b
$$(chsh-rule) 

where $\wedge$ is the logical conjunction and $\oplus$ is the logical disjunction
(equivalently, bitwise addition modulo 2). Therefore one can see for
each possible pair $\left(x,y\right)$ the answer $\left(a,b\right)$
which satisfies the above proposition. If the proposition is true, we
say that Alice and Bob have won the game, and that they have lost
otherwise.

| (x,y)    | Alice and Bob win |
|----------|-------------------|
| (0,0)    | a=b               |
| (1,0)    | a=b               |
| (0,1)    | a=b               |
| (1,1)    | a≠b              |

*Winning scenarios for Alice and Bob*




So, let's say Alice and Bob agree beforehand that whatever bit they
receive, their answer will be $0$, therefore, the same, $a=b$. Thus,
they will win the game $3/4$ of the time. Incidently, the best possible
classical strategy has at most $75\%$ success rate. What is interesting
about this game is that there are quantum strategies with better success
rate. In fact, the best possible quantum strategy has success rate of at
most $\frac{1}{2}+\frac{\sqrt{2}}{2}\approx85\%$. To see this, let us
first explain what we mean by a quantum strategy.

Before the game starts, Alice and Bob prepare the Bell state

$$
\psi=\frac{1}{\sqrt{2}}\left(\left|00\right\rangle +\left|11\right\rangle \right)
$$(bell-state)

and each gets one qubit. For the sake of definiteness, let's say Bob
gets the first one (from the right), while Alice gets the second one
(from the right). Then, they agree to apply certain hermitian unitary
operators (like Pauli matrices) to their qubit depending on the bit they
receive. Let's call these operators $A_{x}$ (for Alice), and $B_{y}$
(for Bob). So, if Bob gets $y=1$, he applies the rotation $B_{1}$ to his
qubit. Their respective answer will be the measured state of their qubit
in the standard computational basis. In terms of these quantities, the
success rate of the CHSH game for a given choice
$S=\left\{ A_{x},B_{y},\psi\right\}$ can be written as 

$$
\begin{aligned}
\omega_{CHSH} & =\frac{1}{8}\sum_{x,y\in\left\{ 0,1\right\} }\left(-1\right)^{x\wedge y}\left\langle A_{x}\otimes B_{y}\right\rangle _{\psi}+\frac{1}{2}\nonumber \\
 & =\frac{1}{8}\left(\left\langle A_{0}\otimes B_{0}\right\rangle _{\psi}+\left\langle A_{0}\otimes B_{1}\right\rangle _{\psi}+\left\langle A_{1}\otimes B_{0}\right\rangle _{\psi}-\left\langle A_{1}\otimes B_{1}\right\rangle _{\psi}\right)+\frac{1}{2}\end{aligned}
$$(success)

In order to define their quantum strategy, let us examine the
expectation values of possible outcomes, in other words, let us examine
the expectation values
$\left\langle A_{x}\otimes B_{y}\right\rangle _{\psi}$ (the subindex
$\psi$ just means that the average value is taken with respect to the
Bell state {eq}`bell-state`. Since $A_{x}$ and $B_{y}$ are hermitian
unitary operators, and therefore square to the identity,
$A_{x}^{2}=B_{x}^{2}=\mathbb{I}$, we have that the square of

$$
T=A_{0}\otimes B_{0}+A_{0}\otimes B_{1}+A_{1}\otimes B_{0}-A_{1}\otimes B_{1}$$
is $$\begin{aligned}
T^{2} & =\left(A_{0}\otimes B_{0}+A_{0}\otimes B_{1}+A_{1}\otimes B_{0}-A_{1}\otimes B_{1}\right)^{2}\\
 & =\left(A_{0}\otimes\left(B_{0}+B_{1}\right)+A_{1}\otimes\left(B_{0}-B_{1}\right)\right)^{2}\\
 & =\mathbb{I}\otimes\left(B_{0}+B_{1}\right)^{2}+\mathbb{I}\otimes\left(B_{0}-B_{1}\right)^{2}+A_{0}A_{1}\otimes\left(B_{0}+B_{1}\right)\left(B_{0}-B_{1}\right)+A_{1}A_{0}\otimes\left(B_{0}-B_{1}\right)\left(B_{0}+B_{1}\right)\\
 & =4\mathbb{I}\otimes\mathbb{I}-\left[A_{0},A_{1}\right]\otimes\left[B_{0},B_{1}\right]\end{aligned}
$$

Now using the following property for bounded linear operators
$\left\Vert A\otimes B\right\Vert \leq\left\Vert A\right\Vert \left\Vert B\right\Vert$
and the fact that unitary operators $A$ satisfy
$\left\Vert A\right\Vert =1$, we have

$$
\left\Vert \left[A_{0},A_{1}\right]\otimes\left[B_{0},B_{1}\right]\right\Vert \leq\left\Vert \left[A_{0},A_{1}\right]\right\Vert \left\Vert \left[B_{0},B_{1}\right]\right\Vert \leq4
$$

where we used the triangle inequality to compute
$\left\Vert A_{0}A_{1}-A_{1}A_{0}\right\Vert \leq\left\Vert A_{0}A_{1}\right\Vert +\left\Vert A_{1}A_{0}\right\Vert \leq2\left\Vert A_{0}\right\Vert \left\Vert A_{1}\right\Vert =2$.
Therefore, we have that

$$
\left\Vert T^{2}\right\Vert \leq8\Rightarrow\left\Vert T\right\Vert \leq2\sqrt{2}
$$

By the Cauchy-Schwarz inequality,

$$
\left|\left\langle \psi,T\psi\right\rangle \right|\leq\left\Vert \psi\right\Vert \left\Vert T\psi\right\Vert \leq\left\Vert T\right\Vert \,,
$$

so 

$$
\left\langle T\right\rangle _{\psi}\leq2\sqrt{2}
$$

Using this last
result in success rate $\omega_{CHSH}$ given in {eq}`success`, we get the bound

$$
\omega_{CHSH}=\frac{1}{8}\left\langle T\right\rangle _{\psi}+\frac{1}{2}\leq\frac{\sqrt{2}+2}{4}\approx85\%\,.
$$(success-rate)

Now let us go back to the problem at hand and figure out what the
operators $A_{x}$ and $B_{y}$ must be. We know they are hermitian
unitary matrices on $\mathbb{C}_{2}$, and a general unitary matrix
acting on a qubit is of the form

$$
U\left(\theta,\phi,\lambda\right)=\left(\begin{array}{cc}
\cos\left(\frac{\theta}{2}\right) & -e^{i\lambda}\sin\frac{\theta}{2}\\
e^{i\phi}\sin\frac{\theta}{2} & e^{i\left(\phi+\lambda\right)}\cos\frac{\theta}{2}
\end{array}\right)\,,\,\,\lambda,\phi\in\left[0,2\pi\right]\,,\,\,\theta\in\left[0,\pi\right]\,.
$$

The condition that $U$ is also hermitian or self-adjoint further imposes
the conditions

$$
e^{i\left(\phi+\lambda\right)}=-1\,,\,\,e^{2i\left(\phi+\lambda\right)}=1\,,
$$

which together imply $\phi+\lambda=k\pi$, for odd $k\in\mathbb{Z}$.
Thus, we are able to eliminate one of the angles, i.e.,
$\lambda=k\pi-\phi$, to obtain

$$
U\left(\theta,\phi\right)=\left(\begin{array}{cc}
\cos\left(\frac{\theta}{2}\right) & e^{-i\phi}\sin\frac{\theta}{2}\\
e^{i\phi}\sin\frac{\theta}{2} & -\cos\frac{\theta}{2}
\end{array}\right)\,,\,\,\phi\in\left[0,2\pi\right]\,,\,\,\theta\in\left[0,\pi\right]\,.$$

Let us use this transformation to compute the probabilites of different
measurement outcomes for Alice and Bob. Let's call the angles Alices
uses $\left(\theta_{a},\phi_{a}\right)$ and those of Bob
$\left(\theta_{b},\phi_{b}\right)$. Then, a direct computation yields

$$
\begin{aligned}
\left[U\left(\theta_{a},\phi_{a}\right)\otimes U\left(\theta_{b},\phi_{b}\right)\right]\psi & =\frac{1}{\sqrt{2}}\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}+e^{-i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)\left|00\right\rangle \\
 & +\frac{1}{\sqrt{2}}e^{i\phi_{b}}\left(\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}-e^{-i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right)\left|01\right\rangle \\
 & +\frac{1}{\sqrt{2}}e^{i\phi_{a}}\left(\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}-\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}e^{-i\left(\phi_{a}+\phi_{b}\right)}\right)\left|10\right\rangle \\
 & +\frac{1}{\sqrt{2}}\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}+e^{i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)\left|11\right\rangle \end{aligned}
$$

The probability that Alice's and Bob's measurements have the same result
is 

$$\begin{aligned}
prob\left(a=b\right) & =\frac{1}{2}\left|\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}+e^{-i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)\right|^{2}+\frac{1}{2}\left|\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}+e^{i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)\right|^{2}\\
 & =\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right)^{2}+2\left(\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)^{2}+\left[e^{i\left(\phi_{a}+\phi_{b}\right)}+e^{-i\left(\phi_{a}+\phi_{b}\right)}\right]\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\\
 & =\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right)^{2}+2\left(\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)^{2}+2\cos\left(\phi_{a}+\phi_{b}\right)\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\end{aligned}$$

For $\cos\left(\phi_{a}+\phi_{b}\right)=1$, this reduces to

$$
prob\left(a=b\right)=\left(\cos\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}+\sin\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)^{2}=\cos^{2}\left(\frac{\theta_{a}-\theta_{b}}{2}\right)\,.
$$(a)


The probability that Alice's and Bob's measurements have a differing
result is 

$$
\begin{aligned}
prob\left(a\neq b\right) & =\frac{1}{2}\left|\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}-e^{-i\left(\phi_{a}+\phi_{b}\right)}\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right|^{2}+\frac{1}{2}\left|\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}-\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}e^{-i\left(\phi_{a}+\phi_{b}\right)}\right|^{2}\\
 & =\left(\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)^{2}+\left(\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right)^{2}-\left[e^{i\left(\phi_{a}+\phi_{b}\right)}+e^{-i\left(\phi_{a}+\phi_{b}\right)}\right]\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\\
 & =\left(\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\right)^{2}+\left(\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\right)^{2}-2\cos\left(\phi_{a}+\phi_{b}\right)\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}\cos\frac{\theta_{a}}{2}\sin\frac{\theta_{b}}{2}\end{aligned}
$$

For $\cos\left(\phi_{a}+\phi_{b}\right)=1$, this reduces to

$$
prob\left(a\neq b\right)=\left(\sin\frac{\theta_{a}}{2}\cos\frac{\theta_{b}}{2}-\sin\frac{\theta_{b}}{2}\cos\frac{\theta_{a}}{2}\right)^{2}=\sin^{2}\left(\frac{\theta_{a}-\theta_{b}}{2}\right)\,
$$(anotb)

For the sake of simplicity, we shall assume $\phi_{a}=\phi_{b}=0$, so
the unitary operators are given by $U(\theta,0,\pi)$

$$
U\left(\theta\right)=\left(\begin{array}{cc}
\cos\left(\frac{\theta}{2}\right) & \sin\frac{\theta}{2}\\
\sin\frac{\theta}{2} & -\cos\frac{\theta}{2}
\end{array}\right)\,,\,\,\theta\in\left[0,\pi\right]\,.
$$

Now that we have obtained the respective probabilities, we can choose
the appropriate operations that Alice and Bob must perform to increase
the likelihood of success. Suppose they are given the pair
$\left(x,y\right)=\left(0,0\right)$. In order to win the game, they must
choose their bits so that $a=b$. The highest chance of success is given
by expression {eq}`success-rate`, so let us equate {eq}`a`  with
{eq}`anotb`:
$\cos^{2}\alpha=\frac{\sqrt{2}+2}{4}\Leftrightarrow\cos2\alpha=\frac{\sqrt{2}}{2}$.
So $\alpha=\frac{\pi}{8}$ or $\alpha=-\frac{\pi}{8}$. Since
$\alpha=\frac{\theta_{a}-\theta_{b}}{2}$, a possible choice is
$\theta_{a}=\frac{\pi}{4}$ and $\theta_{b}=0$. In this case, Alice and
Bob must agree beforehand that if given bit $0$, Alice must apply the
unitary operation $U\left(\frac{\pi}{4}\right)$, while Bob must apply
$U\left(0\right)$. Likewise, for $\left(x,y\right)=\left(0,1\right)$ the
winning strategy is $a=b$, so we must have $\theta_{b}=\frac{\pi}{2}$ so
that $\alpha=-\frac{\pi}{8}$. Hence, Alice still applies
$U\left(\frac{\pi}{4}\right)$ , but Bob applies
$U\left(\frac{\pi}{2}\right)$. For $\left(x,y\right)=\left(1,0\right)$
the winning strategy is again $a=b$, so Bob applies $U\left(0\right)$,
while Alice applies $U\left(-\frac{\pi}{4}\right)$. Finally, for
$\left(x,y\right)=\left(1,1\right)$ Alice and Bob must provided distinct
bits, $a\neq b$. In this case,
$\sin^{2}\alpha=\frac{\sqrt{2}+2}{4}\Leftrightarrow\cos2\alpha=-\frac{\sqrt{2}}{2}$,
so $\alpha=\pm\frac{3\pi}{4}$. Alice applies
$U\left(-\frac{\pi}{4}\right)$ as in the previous case, and Bob applies
$U\left(\frac{\pi}{2}\right)$. In summary, the proposed strategy is 


| (x,y)    | Alice and Bob win | Alice's Gate            | Bob's Gate              |
|----------|-------------------|-------------------------|-------------------------|
| (0,0)    | a=b               | U(π/4)                  | U(0)                    |
| (0,1)    | a=b               | U(π/4)                  | U(π/2)                  |
| (1,0)    | a=b               | U(-π/4)                 | U(0)                    |
| (1,1)    | a≠b              | U(-π/4)                 | U(π/2)                  |

*Winning strategy for Alice and Bob*


       

```{code-cell} ipython3
#the following has been adapted from https://learning.quantum-computing.ibm.com/course/basics-of-quantum-information/entanglement-in-action#the-chsh-game

#definition of the game

from numpy.random import randint


def chsh_game(strategy):
    """Plays the CHSH game
    Args:
        strategy (callable): A function that takes two bits (as `int`s) and
            returns two bits (also as `int`s). The strategy must follow the
            rules of the CHSH game.
    Returns:
        int: 1 for a win, 0 for a loss.
    """
    # Referee chooses x and y randomly
    x, y = randint(0, 2), randint(0, 2)

    # Use strategy to choose a and b
    a, b = strategy(x, y)

    # Referee decides if Alice and Bob win or lose
    if (a != b) == (x & y):
        return 1  # Win
    return 0  # Lose
```

```{code-cell} ipython3
#create the chsh circuit based on the best strategy we presented

from qiskit import QuantumCircuit
from numpy import pi


def chsh_circuit(x, y):
    """Creates a `QuantumCircuit` that implements the best CHSH strategy.
    Args:
        x (int): Alice's bit (must be 0 or 1)
        y (int): Bob's bit (must be 0 or 1)
    Returns:
        QuantumCircuit: Circuit that, when run, returns Alice and Bob's
            answer bits.
    """

#we apply the general unitary gates qc.u  with parameters phi=0 and lambda=pi
    
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.barrier()

    # Alice
    if x == 0:
        qc.u(pi/4, 0,pi,0)
    else:
        qc.u(-pi / 4, 0, pi, 0)
    qc.measure(0, 0)

    # Bob
    if y == 0:
        qc.u(0,0,pi, 1)
    else:
        qc.u(pi / 2,0,pi, 1)
    qc.measure(1, 1)

    return qc
```

```{code-cell} ipython3
#run the circuit on simulator for fixed pair (x,y)

from qiskit_aer.primitives import Sampler

sampler = Sampler()

def quantum_strategy(x, y):
    """Carry out the best strategy for the CHSH game.
    Args:
        x (int): Alice's bit (must be 0 or 1)
        y (int): Bob's bit (must be 0 or 1)
    Returns:
        (int, int): Alice and Bob's answer bits (respectively)
    """
    # `shots=1` runs the circuit once
    result = sampler.run(chsh_circuit(x, y), shots=1).result()
    statistics = result.quasi_dists[0].binary_probabilities()
    bits = list(statistics.keys())[0]
    a, b = bits[0], bits[1]
    return a, b
```

```{code-cell} ipython3

NUM_GAMES = 1000
TOTAL_SCORE = 0

for _ in range(NUM_GAMES):
    TOTAL_SCORE += chsh_game(quantum_strategy)

print("Fraction of games won:", TOTAL_SCORE / NUM_GAMES)
```

```{code-cell} ipython3

```

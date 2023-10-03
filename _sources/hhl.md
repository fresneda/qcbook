---
title: HHL algorithm
---

# Overview

We will discuss a hybrid classical-quantum algorithm to solve the linear
system of equations $$Ax=b\,,\label{eq:linear-sytem}$$ where
$A\in\mathcal{M}_{N}\left(\mathbb{C}\right)$ is a complex matrix of
order $N$, following the HHL algorithm.

Without lost of generality, we will assume that
$\left\Vert b\right\Vert =1$ where
$\left\Vert b\right\Vert =\sqrt{b^{*}b}$ is the usual hermitian norm in
$\mathbb{C}^{N}$ ($b^{*}$ is the conjugate transpose of $b$), and that
$A$ is hermitian, $A=A^{*}$. The first assumption can be made by simply
rescaling of the solution by $\left\Vert b\right\Vert$, while the second
assumption be realized in terms of the system $$\left(\begin{array}{cc}
0 & A\\
A^{*} & 0
\end{array}\right)\left(\begin{array}{c}
0\\
x
\end{array}\right)=\left(\begin{array}{c}
b\\
0
\end{array}\right)$$ which is equivalent to $Ax=b$ and has an hermitian
coefficient matrix.

Let us further assume that $A$ is invertible, that is, the linear system
$Ax=b$ has the unique solution $x=A^{-1}b$. This is equivalent to the
assumption that all eigenvalues of $A$ are non-zero. By the spectral
theorem for hermitian matrices, the eigenvalues
$\lambda_{0}$,\...,$\lambda_{N-1}$ of $A$ lie in the compact set
$\left[\lambda_{min},\lambda_{max}\right]\subset\mathbb{R}$ with
$\lambda_{min}>0$. In addition, $A$ has an orthonormal basis of
eigenvectors
$\left\{ \left|v_{0}\right\rangle ,...,\left|v_{N-1}\right\rangle \right\}$,
which allows $A$ to be written in terms of the orthogonal projection
operators
$$A=\sum_{i=0}^{N-1}\lambda_{i}\left|v_{i}\right\rangle \left\langle v_{i}\right|\,.$$
Likewise, the inverse matrix $A^{-1}$ can be written as
$$A^{-1}=\sum_{i=0}^{N-1}\lambda_{i}^{-1}\left|v_{i}\right\rangle \left\langle v_{i}\right|\,,$$
so for $b$ given by the linear combination
$$b=\sum_{i=0}^{N-1}b_{i}\left|v_{i}\right\rangle$$ the solution to the
system $Ax=b$ is given by
$$x=A^{-1}b=\sum_{i=0}^{N-1}\frac{b_{i}}{\lambda_{i}}v_{i}\label{eq:solution}$$
Thus, the problem of finding the solution of $Ax=b$ reduces to the
problem of diagonalizing $A$, or rather, of finding its eigenvalues.

# Pre-processing

In order to express this problem in terms of an $n$-qubit quantum
system, we need to extend $A$ to a $2^{n}\times2^{n}$ matrix by
completing it with the identity, so the linear system becomes
$$\left(\begin{array}{cc}
A & 0\\
0 & I
\end{array}\right)\left(\begin{array}{c}
x\\
0
\end{array}\right)=\left(\begin{array}{c}
b\\
0
\end{array}\right)$$ For instance, if
$A\in M_{3}\left(\mathbb{C}\right)$, then we consider the alternative
system $$\left(\begin{array}{cc}
A & 0\\
0 & 1
\end{array}\right)\left(\begin{array}{c}
x\\
0
\end{array}\right)=\left(\begin{array}{c}
b\\
0
\end{array}\right)\,,$$ where the coefficient matrix is now of order
$2^{2}$.

Another key adaptation for a quantum system is the necessity of having
the eigenvalues in the interval $\left(0,1\right]$. This can be
accomplished by rescaling and shifting the spectrum of $A$. However, it
is important to accomplish this without prior knowledge of the point
spectrum of $A$. This can be achieved by considering the spectral radius
$\rho\left(A\right)$ of $A$, which is defined by
$$\rho\left(A\right)=\max\left\{ \left|\lambda_{1}\right|,...\text{\ensuremath{\left|\lambda_{n}\right|}}\right\} \,.$$
Because $\rho\left(A\right)\leq\left\Vert A\right\Vert$ for any operator
norm $\left\Vert \cdot\right\Vert$, we can choose a convenient matrix
norm, such as the maximum absolute row sum,
$$\rho\equiv\left\Vert A\right\Vert _{\infty}=\max_{1\leq i\leq N}\left\{ \sum_{j=1}^{N}\left|a_{ij}\right|\right\}$$
so that
$spec\left(A\right)=\left[\lambda_{min},\lambda_{max}\right]\subset\left[-\rho,\rho\right]$.
Consider the matrix $$\tilde{A}=\Omega\left(A+\omega I\right)$$ where
$\Omega>0$ and $\omega$ are real numbers. The eigenvalues of $\tilde{A}$
are the $\tilde{\lambda}$ roots of the characteristic polynomial
$$\det\left(\tilde{A}-\tilde{\lambda}I\right)=\det\left(\Omega\left(A+\omega I\right)-\tilde{\lambda}I\right)=\Omega^{N}\det\left(A+\omega-\frac{\tilde{\lambda}}{\Omega}\right)\,.$$
The last determinant tells us that
$$\lambda=\frac{\tilde{\lambda}}{\Omega}-\omega\Leftrightarrow\tilde{\lambda}=\Omega\left(\omega+\lambda\right)\label{eq:eigenvalue}$$
are the eigenvalues of $A$, so we find a relation between the spectrum
of $A$ and that of $\tilde{A}$. Let us now impose the condition that the
spectrum of $\tilde{A}$ is in $\left(0,1\right]$: $$\begin{aligned}
\tilde{\lambda}\left(\lambda=-\rho\right) & =\Omega\left(\omega-\rho\right)=\frac{1}{2^{k}}\\
\tilde{\lambda}\left(\lambda=\rho\right) & =\Omega\left(\omega+\rho\right)=1-\frac{1}{2^{k}}\end{aligned}$$
for integer $k>1$. The solution for any $k$ is
$$\Omega=\left(\frac{1}{2}-\frac{1}{2^{k}}\right)\rho^{-1}\,,\,\,\omega=\frac{1}{1-2^{1-k}}\rho$$

For instance, consider the Hadamard matrix
$$H=\frac{1}{\sqrt{2}}\left(\begin{array}{cc}
1 & 1\\
1 & -1
\end{array}\right)\,.$$ It is hermitian and
$\left\Vert H\right\Vert _{\infty}=\sqrt{2}$, so
$spec\left(H\right)\subset\left[-\sqrt{2},\sqrt{2}\right]$, which is
true, of course, since its eigenvalues are $\pm1$. For $k=2$, we have
$\Omega=\sqrt{2}/8$ and $\omega=2\sqrt{2}$, so
$$\tilde{H}=\Omega\left(H+\omega I\right)=\frac{1}{8}\left(\begin{array}{cc}
5 & 1\\
1 & 3
\end{array}\right)\,.$$ One can check that the eigenvalues of
$\tilde{H}$ are $\frac{1}{8}\left(4\pm\sqrt{2}\right)$, so
$spec\left(\tilde{H}\right)\subset\left[\frac{1}{4},\frac{3}{4}\right]\subset\left(0,1\right]$.

Finally, it is important to notice that the eigenvector of $A$ belonging
to the eigenvalue $\lambda$ is also an eigenvector of $\tilde{A}$
belonging to $\tilde{\lambda}$. In fact,
$$\tilde{A}\left|v\right\rangle =\text{\ensuremath{\tilde{\lambda}\left|v\right\rangle \Leftrightarrow A\left|v\right\rangle =\frac{1}{\Omega}\left(\tilde{\lambda}-\Omega\omega\right)\left|v\right\rangle =\lambda\left|v\right\rangle }}$$
Therefore, from the solution
$$\tilde{x}=\tilde{A}^{-1}b=\sum_{i}\frac{b_{i}}{\tilde{\lambda}_{i}}\left|v_{i}\right\rangle$$
we get the solution
$$x=A^{-1}b=\sum_{i}\frac{b_{i}}{\lambda_{i}}\left|v_{i}\right\rangle$$
by simple replacement of each $\tilde{\lambda}_{i}$ for $\lambda_{i}$,
according to ([\[eq:eigenvalue\]](#eq:eigenvalue){reference-type="ref"
reference="eq:eigenvalue"}).

# hybrid HHL

Suppose all the necessary pre-processing has been deatlt with, and we
are to apply HHL to the linear system
([\[eq:linear-sytem\]](#eq:linear-sytem){reference-type="ref"
reference="eq:linear-sytem"}) where $A$ is an hermitian matrix of order
$2^{n}$, whose eigenvalues lie in the interval $\left(0,1\right]$.

Since we can only implement unitary transformations, we consider instead
the unitary matrix $U=e^{2\pi iA}$, whose eigenvalues are
$e^{2\pi i\lambda_{i}}$, $i=0,...,N-1$. So let's say we have prepared a
$n$-qubit state $\left|v\right\rangle$ which is an eigenstate of $A$
with eigenvalue $\lambda$, so
$$U\left|v\right\rangle =e^{2\pi i\lambda}\left|v\right\rangle \,.$$ We
then apply the Quantum Phase Estimation (QPE) algorithm to the state
$\left|v\right\rangle \otimes\left|0\right\rangle ^{\otimes t}$, where
$\left|0\right\rangle ^{\otimes t}$ is a second quantum register with
$t$ qubits. The result is
$$\left|v\right\rangle \otimes\left|0\right\rangle ^{\otimes t}\overset{QPE}{\mapsto}\left|v\right\rangle \otimes\left|\tilde{\lambda}\right\rangle$$
where $\left|\tilde{\lambda}\right\rangle$ is the state which outputs
$\lambda$ when measured with the accuracy of $t$ qubits, that is,
$\tilde{\lambda}$ is an approximation of $\lambda$ with $t$ bits of
accuracy. So one would measure $\lambda$ exactly if $2^{t}\lambda$ is an
integer. If we apply QPE to the state
$\left|b\right\rangle \left|0\right\rangle$, then we should arrive at
$$\left|b\right\rangle \left|0\right\rangle ^{\otimes t}=\sum_{i=0}^{N-1}b_{i}\left|v_{i}\right\rangle \left|0\right\rangle ^{\otimes t}\overset{QPE}{\mapsto}\sum_{i=0}^{N-1}b_{i}\left|v_{i}\right\rangle \left|\tilde{\lambda}_{i}\right\rangle$$
where $\tilde{\lambda}_{i}$ is the $t$-bit approximation for the
eigenvalue $\lambda_{i}$.

The next step is to get the reciprocals of the eigenvalues, which appear
in the solution of the linear system
([\[eq:solution\]](#eq:solution){reference-type="ref"
reference="eq:solution"}). The idea is to introduce an auxiliar qubit
$\left|0\right\rangle$ and perform a $y$-rotation $R_{y}$ by an
appropriate angle $\theta$. Since
$$R_{y}\left(\theta\right)\left|0\right\rangle =\cos\frac{\theta}{2}\left|0\right\rangle +\sin\frac{\theta}{2}\left|1\right\rangle$$
we need $\theta$ to be $2\arcsin f\left(\lambda\right)$ (we could also
work with $\arccos$), where
$f\left(\lambda\right)=\frac{\Lambda}{\lambda}$ and $\Lambda$ is a
scaling parameter. In this way,
$$R_{y}\left(2\arcsin\frac{\Lambda}{\lambda}\right)=\sqrt{1-\frac{\Lambda^{2}}{\lambda^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda}\left|1\right\rangle$$
and by measuring the auxiliar qubit we can find the reciprocal of the
eigenvalue, up to a scaling factor. Since all eigenvalues belong to
$\left[\lambda_{min},\lambda_{max}\right]$,
$\lambda\in\left[\lambda_{min},\lambda_{max}\right]$, and so
$\frac{\Lambda}{\lambda}\in\left[\frac{\Lambda}{\lambda_{max}},\frac{\Lambda}{\lambda_{min}}\right]$.
In order to use the function $\arccos$ must choose the scale $\Lambda$
so that
$\left[\frac{\Lambda}{\lambda_{max}},\frac{\Lambda}{\lambda_{min}}\right]$
fits the domain of $\arcsin$. This can by accomplished by rescaling by
$\Lambda\leq\lambda_{min}$, so we get
$\left[\frac{\Lambda}{\lambda_{max}},\frac{\Lambda}{\lambda_{min}}\right]=\left[\frac{\lambda_{min}}{\lambda_{max}},1\right]\subset\left(0,1\right]$.
This means
$\arcsin\left(\frac{\Lambda}{\lambda}\right)\in\left(0,\frac{\pi}{2}\right]$,
or equivalently,
$\frac{2}{\pi}\arcsin\left(\frac{\Lambda}{\lambda}\right)\in\left(0,1\right]$.

Therefore, if we could apply rotations
$R_{y}\left(2\arcsin\frac{\Lambda}{\lambda_{i}}\right)$ to the auxiliar
qubit, conditioned on the register
$\left|\tilde{\lambda}_{i}\right\rangle$, we would get the desired
result of inverting the eigenvalues, and arrive at the state
$$\sum_{i=0}^{N-1}b_{i}\left(\sqrt{1-\frac{\Lambda^{2}}{\lambda_{i}^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \right)\left|v_{i}\right\rangle \left|\tilde{\lambda}_{i}\right\rangle$$
The final step is to apply the inverse QPE to bring back the register
$\left|\tilde{\lambda}_{i}\right\rangle$ to
$\left|0\right\rangle ^{\otimes t}$:
$$\sum_{i=0}^{N-1}b_{i}\left(\sqrt{1-\frac{\Lambda^{2}}{\lambda_{i}^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \right)\left|v_{i}\right\rangle \left|\tilde{\lambda}_{i}\right\rangle \overset{QPE^{\dagger}}{\mapsto}\sum_{i=0}^{N-1}b_{i}\left(\sqrt{1-\frac{\Lambda^{2}}{\lambda_{i}^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \right)\left|v_{i}\right\rangle \left|0\right\rangle ^{\otimes t}$$
Now, measurement values of $1$ of the auxiliar qubit will result in the
state
$$\sum_{i=0}^{N-1}b_{i}\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \left|v_{i}\right\rangle \left|0\right\rangle ^{\otimes t}$$

A solution to the problem of efficiently implementing a quantum circuit
to perform eigenvalue inversion does not seem to be a hand. In
([ref1](https://arxiv.org/abs/1207.2485)), the eigenvalue inversion is
done using a combination of Newton's method and bissection, while
([ref2](https://arxiv.org/pdf/2009.04484.pdf)) and
([ref3](https://arxiv.org/pdf/1805.12445.pdf)) use piecewise polynomial
approximation to approximate $\arcsin$. We will adopt a third approach,
which is a hybrid approach, and involves first measuring the QPE circuit
multiple times to get a probability distribution of the eigenvalues.
Having computed the eigenvalues, we can apply the rotations controlled
rotations having the eigenvalues as control. To see how this works, let
the first measured eigenvalue be $\lambda_{1}$, which is stored in
$\left|\lambda_{1}\right\rangle$. So we apply a controlled rotation
$Ry\left(2\arcsin\frac{\Lambda}{\lambda_{1}}\right)$ to the target
qubit, whose control bits are the bits of the measured eigenvalue
$\lambda_{1}$. Repeat for each measured eigenvalue, and we are done.

Putting everything together, we have $$\begin{aligned}
\left|0\right\rangle \left|b\right\rangle \left|0\right\rangle ^{\otimes t}\overset{QPE}{\mapsto}\sum_{i=0}^{N-1}b_{i}\left|0\right\rangle \left|v_{i}\right\rangle \left|\tilde{\lambda}_{i}\right\rangle  & \overset{R_{y}}{\mapsto}\sum_{i=0}^{N-1}b_{i}\left(\sqrt{1-\frac{\Lambda^{2}}{\lambda_{i}^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \right)\left|v_{i}\right\rangle \left|\tilde{\lambda}_{i}\right\rangle \\
 & \overset{QPE^{\dagger}}{\mapsto}\sum_{i=0}^{N-1}b_{i}\left(\sqrt{1-\frac{\Lambda^{2}}{\lambda_{i}^{2}}}\left|0\right\rangle +\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \right)\left|v_{i}\right\rangle \left|0\right\rangle ^{\otimes t}\\
 & \overset{measurement}{\mapsto}\sum_{i=0}^{N-1}b_{i}\frac{\Lambda}{\lambda_{i}}\left|1\right\rangle \left|v_{i}\right\rangle \left|0\right\rangle ^{\otimes t}\end{aligned}$$

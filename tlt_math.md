### Two-Level Targeter (TLT)

By dividing the corrections process into different parts, the size of the Jacobian and, therefore, the computational load in calculating its matrix inverse is decreased. Additionally, monotonic convergence behavior is more frequently observed in multi-stage corrections processes than in standard multiple shooting algorithms. In on-board applications where communication with ground stations may be limited or lost, the predictability in convergence gained from a staged targeting algorithm is valuable. One such corrections strategy with two levels is called a *Two-Level Targeter (TLT)*. As originally formulated, the TLT is decomposed into Level-I which corrects for position continuity and Level-II which targets velocity continuity. To model maneuvers at patch points by chemical propulsion systems that occur over time intervals significantly smaller than the arc propagation times, the thrusting acceleration over a small duration may be represented as an instantaneous velocity change, labelled a Delta-V or impulsive maneuver. A two-level corrections algorithm that uses this assumption is named an impulsive TLT.

#### 3.4.1 Impulsive Level-I

In Level-I of an impulsive TLT algorithm, position continuity within the convergence tolerance is achieved for each arc in series. An initially discontinuous set of arcs is propagated starting with the first patch point. The spacecraft trajectory is propagated from the initial patch point, $p_o$, at time $t_o$, to the next patch point, $p_p$, at time $t_p$. The design vector,

$$
X = v^+_o
$$

> (3.53)

consists of the outgoing velocity vector at $p_o$. This quantity is varied within a single shooting algorithm to achieve position continuity at $p_p$, mathematically described by the constraint vector,

$$
F = p^−_p - p^+_p = 0
$$

> (3.54)

<p align="center">
<img src="https://i.imgur.com/KL1Xp48.png">
</p>

In practical terms, Level-I represents an impulsive maneuver to target a new upstream position. The remaining component necessary to formulate an iterative Newton algorithm is the Jacobian,

$$
DF = \frac{∂p^+_p}{∂v^+_o}
$$

> (3.55)

that models the outgoing position at $p_p$ as a constant. The STM is used to solve an iterative single shooting algorithm. The initial position, $p^+_o$, and patch point times are constant, resulting in $\delta p^+_o = 0$ and $\delta t_o = dt_p = 0$. By employing the STM and simplifying, the variations at $p_o$ and $p_p$ are related by

$$
\begin{align}
\begin{bmatrix} d p^+_p \\ d v^+_o \end{bmatrix} = Φ(t_p, t_o) = \begin{bmatrix} A_op & B_op \\ C_op & D_op \end{bmatrix}
\end{align}
$$

> (3.56); the version of LaTeX equations that GitHub .md files use is a bit odd, apologies for the odd stacking

This relationship leads to determination of the Jacobian,

$$
\frac{∂p^+_p}{∂v^+_o} = B_op
$$

> > (3.57)

which is substituted into the previous equation to yield

$$
d v^+_o = −B^{−1}_op d p^+_p
$$

> (3.58)

Each iteration in the single shooting algorithm is denoted a *local iteration*. Once the $F$ constraint is within the specified convergence tolerance, the algorithm is repeated for the next patch point in the series. The previous $p_p$ becomes the new $p_o$, and the next patch point is renamed $p_p$. Once all the patch points have reached convergence through these sequential single shooting processes, Level-I is finished for a single *global iteration*.

#### 3.4.2 Impulsive Level-II

Following completion of Level-I for a global iteration, the algorithm executes one multiple shooting correction in Level-II. The goal of Level-II is to reduce velocity discontinuities by shifting patch point positions and times. For a trajectory comprised of $k$ patch points, velocity continuity at the interior patch points is enforced by creating a constraint vector, i.e.,

$$
F = \begin{bmatrix} v^2 − v^+_2 \\ \vdots \\ v^−_k−1 − v^+_k−1 \end{bmatrix}
$$

As apparent in the provided content, all patch points may be considered continuous in position, expressed as

$$
p_i = p^+_i = p^−_i
$$

and, therefore, no superscript is necessary for the position vectors. Thus, the design variable vector is formulated as

$$
X = \begin{bmatrix} p_1 \\ t_1 \\ \vdots \\ p_k \\ t_k \end{bmatrix}
$$

<p align="center">
<img src="https://imgur.com/NZFgrKC">
</p>

Yes! Let me explain duality in T-splines using mathematical notation:

Basic Definitions :

Let $S$ be a spline space

Let $S^*$ be its dual space (linear functionals on $S$)

For $f \in S$ and $\lambda \in S^*$, we denote the duality pairing as $\langle \lambda, f \rangle$

B-spline Basis and Dual Functionals :

B-spline basis: ${N_i}_{i \in I}$

Dual functionals: ${\lambda_i}_{i \in I}$

Biorthogonality: $\langle \lambda_i, N_j \rangle = \delta_{ij}$ where $\delta_{ij}$ is Kronecker delta

Blossoming Connection :

For a spline $f$ of degree $n$, its blossom $f^b$ is an n-affine symmetric function

The dual functional $\lambda_i$ can be expressed using the blossom: $\lambda_i(f) = f^b(u_{i,1}, ..., u_{i,n})$ where $u_{i,j}$ are knots associated with basis function $i$

T-spline Case :

Local knot vectors: $(u_{i,1}, ..., u_{i,p+1})$ and $(v_{i,1}, ..., v_{i,q+1})$

T-spline basis function: $T_i(u,v) = N_i^u(u)N_i^v(v)$ where $N_i^u$ and $N_i^v$ are univariate B-splines

Dual Functionals for T-splines : $\lambda_i(f) = f^b(u_{i,1}, ..., u_{i,p+1}; v_{i,1}, ..., v_{i,q+1})$

Linear Independence :

Basis functions ${T_i}$ are linearly independent if: $\sum_i c_i T_i = 0 \implies c_i = 0 ;\forall i$

Using duals: $\det(\langle \lambda_i, T_j \rangle)_{i,j} \neq 0$

Analysis-suitable Conditions : For T-splines to be analysis-suitable:

No intersecting T-junction extensions

Dual functionals form a dual basis

Partition of unity: $\sum_i \lambda_i(1) = 1$

Refinement Using Duals : For a refinement operation: $T_i = \sum_j \alpha_{ij} \tilde{T} j$ The dual coefficients satisfy: $\alpha {ij} = \langle \lambda_i, \tilde{T}_j \rangle$

Geometric Properties : For a T-spline surface: $S(u,v) = \sum_i P_i T_i(u,v)$ The control points can be recovered: $P_i = \lambda_i(S)$

Smoothness Analysis :

$C^k$ continuity can be characterized through dual functionals: $\langle \lambda_i, D^\alpha T_j \rangle = 0$ for $|\alpha| \leq k$ at appropriate points

Key Properties:

Duality Pairing : $\langle \lambda, \alpha f + \beta g \rangle = \alpha\langle \lambda, f \rangle + \beta\langle \lambda, g \rangle$

Biorthogonality : $\langle \lambda_i, T_j \rangle = \delta_{ij}$

Coefficient Extraction : For $f = \sum_i c_i T_i$: $c_i = \langle \lambda_i, f \rangle$

Partition of Unity : $\sum_i T_i = 1 \iff \sum_i \lambda_i = \delta$ (Dirac functional)

Support Properties : $supp(\lambda_i) \subseteq$ supp$(T_i)$

This mathematical framework provides:

Clear theoretical foundation

Connection between blossoming and duality

Tools for analysis and implementation

Basis for refinement algorithms

Framework for extending to other spline types

The dual perspective unifies many aspects of T-splines and provides powerful tools for both analysis and implementation

For a B-spline of degree

$N_i^n(t) = \frac{t - t_i}{t_{i+n} - t_i}N_i^{n-1}(t) + \frac{t_{i+n+1} - t}{t_{i+n+1} - t_{i+1}}N_{i+1}^{n-1}(t)$

Its blossom $ b_i^n(u_1,...,u_n)$ satisfies:

• Symmetry: $ b*i^n(u_1,...,u_n) = b_i^n(u*{\pi(1)},...,u\_{\pi(n)}) $ for any permutation $ \pi$
•Multi-affine property
• Diagonal property: $ b_i^n(t,...,t) = N_i^n(t)$

Dual Functionals and B-splines :

The dual functional $ \lambda_i$ for a B-spline basis function satisfies:

$\langle \lambda_i, N_j \rangle = \delta_{ij}$

Can be expressed using blossoming:

$\lambda_i(f) = f^b(t_{i+1},...,t_{i+n})$

where $ f^b $ is the blossom of $ f $ and $ t_i $ are knots

Connection Through de Boor Algorithm :

De Boor points $ d_i^r $ can be expressed using blossoming:

$d_i^r = b(t_{i+1},...,t_{i+r},\underbrace{t,...,t}_{n-r})$

The dual functional relates to these points:

$\lambda_i(f) = d_i^n$ (final de Boor point)

Recursive Relations :

For blossoming:

$ b(u_1,...,u_n) = \alpha b(u_2,...,u_n,t_L) + (1-\alpha)b(u_2,...,u_n,t_R)$

For dual functionals:

$\lambda_i = \alpha\lambda_{i-1} + (1-\alpha)\lambda_i$

where $ \alpha $ is determined by knot differences

Properties in Terms of Support :

Support of B-spline:
$\text{supp}(N_i^n) = [t_i, t_{i+n+1}]$

Support of dual functional:
$\text{supp}(\lambda_i) \subseteq [t_{i+1}, t_{i+n}]$

Linear Independence :

Through blossoming:
$\det(b_i(t_{j+1},...,t_{j+n}))_{i,j} \neq 0$

Through duals:
$\det(\langle \lambda_i, N_j \rangle)_{i,j} \neq 0$

Coefficient Extraction :

For a spline $ f = \sum_i c_i N_i:$

$c_i = \lambda_i(f) = f^b(t_{i+1},...,t_{i+n})$

Extension to T-splines :

For T-spline basis function:
$T_i(u,v) = N_i^u(u)N_i^v(v)$

Dual functional:
$\lambda_i(f) = f^b(u_{i,1},...,u_{i,p+1};v_{i,1},...,v_{i,q+1})$

Smoothness Conditions:

Through blossoming:
$C^k $ continuity $ \iff b(u_1,...,u_{n-k-1},\underbrace{t,...,t}_{k+1}) \text{ continuous in } t$

$\text{Through duals:}
\langle \lambda_i, D^\alpha f \rangle = 0 \text{ for } |\alpha| \leq k$

Geometric Interpretation :

Control points:
$P_i = \lambda_i(S) = S^b(t_{i+1},...,t_{i+n})$

De Boor points:
$d_i^r(t) = S^b(t_{i+1},...,t_{i+r},\underbrace{t,...,t}_{n-r})$

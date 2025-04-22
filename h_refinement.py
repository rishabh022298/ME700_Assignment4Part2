####################################################
# THERMAL EXPANSION - LINEAR ELASTICITY (P2: h-REFINEMENT)
####################################################

from dolfinx import log, default_scalar_type
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import numpy as np
import ufl
from ufl import sym, grad, Identity, tr, inner
from mpi4py import MPI
from dolfinx import fem, mesh
import matplotlib.pyplot as plt

# MATERIAL PARAMETERS FOR STEEL
E = default_scalar_type(210.0e9)                                    # Young's modulus
nu = default_scalar_type(0.3)                                       # Poisson's ratio
alpha = default_scalar_type(12e-6)                                  # Coefficient of thermal expansion

# Geometry
L = 20.0  # Length of beam
W = 1.0   # Width and height

# Refinement levels (can be adjusted)
h_refinement_levels = [(1, 5, 5), (2, 5, 5), (4, 5, 5), (8, 5, 5), (16, 5, 5), (20, 5, 5)]

max_disp_h = []

T_applied = 200.0
u_analytical_max = alpha * T_applied * L


def solve_problem(h_level, degree):
    # Mesh creation
    domain = mesh.create_box(MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, W, W]], h_level, mesh.CellType.hexahedron)

    # Function space
    V = fem.functionspace(domain, ("Lagrange", degree, (domain.geometry.dim, )))

    # Boundary markers
    def left(x): return np.isclose(x[0], 0)
    fdim = domain.topology.dim - 1
    left_facets = mesh.locate_entities_boundary(domain, fdim, left)
    facet_tag = mesh.meshtags(domain, fdim, left_facets, np.full_like(left_facets, 1))
    u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
    left_dofs = fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
    bcs = [fem.dirichletbc(u_bc, left_dofs, V)]

    # Constants
    B = fem.Constant(domain, default_scalar_type((0, 0, 0)))
    T = fem.Constant(domain, default_scalar_type(T_applied))

    # Material constants
    mu = fem.Constant(domain, E / (2 * (1 + nu)))
    lmbda = fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))

    # Strain and stress
    def epsilon(u):
        return sym(grad(u)) - alpha * Identity(len(u)) * T

    def sigma(u):
        return 2.0 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(len(u))

    # Weak form
    v = ufl.TestFunction(V)
    u = fem.Function(V)
    dx = ufl.Measure("dx", domain=domain, metadata={"quadrature_degree": 4})
    F_form = inner(sigma(u), grad(v)) * dx - inner(v, B) * dx

    # Solver
    problem = NonlinearProblem(F_form, u, bcs)
    solver = NewtonSolver(domain.comm, problem)
    solver.atol = 1e-8
    solver.rtol = 1e-8
    solver.convergence_criterion = "incremental"
    solver.max_it = 100

    log.set_log_level(log.LogLevel.OFF)
    num_its, converged = solver.solve(u)
    assert converged
    u.x.scatter_forward()

    # Compute maximum displacement magnitude
    Vs = fem.functionspace(domain, ("Lagrange", degree))
    mag = fem.Function(Vs)
    expr = fem.Expression(ufl.sqrt(sum(u[i]**2 for i in range(len(u)))), Vs.element.interpolation_points())
    mag.interpolate(expr)
    return np.max(mag.x.array)

# h-refinement plot
fixed_degree = 2
mesh_sizes = [h[0] for h in h_refinement_levels]
for h_level in h_refinement_levels:
    max_disp = solve_problem(h_level, fixed_degree)
    max_disp_h.append(max_disp)

plt.figure(figsize=(8, 6))
plt.plot(mesh_sizes, max_disp_h, marker='o', linewidth=2.5, label="Numerical")
plt.axhline(y=u_analytical_max, color='r', linestyle='--', linewidth=2, label="Analytical")
plt.xlabel("Number of elements in x-direction", fontsize=14)
plt.ylabel("Max displacement magnitude", fontsize=14)
plt.title("h-refinement: Max Displacement vs Mesh Size", fontsize=16)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig("P2_h_refinement.png", dpi=300)
plt.close()

F = GF(
    21888242871839275222246405745257275088696311157297823662689037894645226208583
)


def build_M_old(k, s, base_field=F):
    if k <= 0 or s <= 0:
        raise ValueError("k and s must be positive integers")

    base_vars = ["t", "a", "c", "x"]
    l_vars = [f"L_{i}_{j}" for i in range(k) for j in range(4)]
    r_vars = [f"R_{i}_{j}" for i in range(4) for j in range(k)]
    names = base_vars + l_vars + r_vars

    R0 = PolynomialRing(base_field, names=names)
    P = R0.fraction_field()
    gens = P.gens_dict()

    t = gens["t"]
    a = gens["a"]
    c = gens["c"]
    b = c^2 / a      #to enforce ab = c^2
    x = gens["x"]

    L = Matrix(P, k, 4, [gens[f"L_{i}_{j}"] for i in range(k) for j in range(4)])
    R = Matrix(P, 4, k, [gens[f"R_{i}_{j}"] for i in range(4) for j in range(k)])

    U = Matrix(P, 
        [[a, c, -x, P(0)],
         [c, b, P(0), x ],
         [P(0), P(0), b, c],
         [P(0), P(0), c, a]])

    # Coefficient arrays now have length 2*s - 1.
    gamma = [t**m for m in range(2 * s - 1)]
 
    Gamma = Matrix(P, s, s, lambda i, j: gamma[i + j])
    
    G = L * U * R

    M = G.tensor_product(Gamma)
    M_flat = Matrix(
        P,
        k * k,
        2 * s - 1,
        lambda row, ell: G[row // k, row % k] * gamma[ell],
    )
    M_flat_polys = [M_flat[i, j] for i in range(M_flat.nrows()) for j in range(M_flat.ncols())]
    vars = P.gens()
    jacobian_flat = Matrix(
        P,
        len(M_flat_polys),
        len(vars),
        lambda i, j: M_flat_polys[i].derivative(vars[j]),
    )

    return {
        "ring": P,
        "M": M,
        "M_flat": M_flat,
        "M_flat_polys": M_flat_polys,
        "jacobian_flat": jacobian_flat,
        "vars": vars,
        "G": G,
        "Gamma": Gamma,
        "gamma": gamma,
        "L": L,
        "R": R,
    }


def build_M_new(k, s, base_field=F, with_scale=False):
    if k <= 0 or s <= 0:
        raise ValueError("k and s must be positive integers")

    base_vars = ["t", "u", "a", "c", "x"]
    l_vars = [f"L_{i}_{j}" for i in range(k) for j in range(2)]
    r_vars = [f"R_{i}_{j}" for i in range(2) for j in range(k)]
    names = base_vars + l_vars + r_vars

    R0 = PolynomialRing(base_field, names=names)
    P = R0.fraction_field()
    gens = P.gens_dict()

    t = gens["t"]
    u = gens["u"]
    a = gens["a"]
    c = gens["c"]
    b = c^2 / a     # to enforce ab = c^2
    x = gens["x"]

    L = Matrix(P, k, 2, [gens[f"L_{i}_{j}"] for i in range(k) for j in range(2)])
    R = Matrix(P, 2, k, [gens[f"R_{i}_{j}"] for i in range(2) for j in range(k)])

    U0 = Matrix(P, [[a, c], [c, b]])
    U1 = Matrix(P, [[P(0), -x], [x, P(0)]])

    # Coefficient arrays now have length 2*s - 1.
    gamma0 = [t**m for m in range(2 * s - 1)]
    gamma1 = [P(0)] + [m * t**(m - 1) for m in range(1, 2 * s - 1)]
    if with_scale:
        gamma1 = [gamma1[i] + u * gamma0[i] for i in range(0, 2 * s - 1)]

    Gamma0 = Matrix(P, s, s, lambda i, j: gamma0[i + j])
    Gamma1 = Matrix(P, s, s, lambda i, j: gamma1[i + j])

    G0 = L * U0 * R
    G1 = L * U1 * R

    M = G0.tensor_product(Gamma1) + G1.tensor_product(Gamma0)
    M_flat = Matrix(
        P,
        k * k,
        2 * s - 1,
        lambda row, ell: G0[row // k, row % k] * gamma1[ell]
        + G1[row // k, row % k] * gamma0[ell],
    )
    M_flat_polys = [M_flat[i, j] for i in range(M_flat.nrows()) for j in range(M_flat.ncols())]
    vars = P.gens()
    jacobian_flat = Matrix(
        P,
        len(M_flat_polys),
        len(vars),
        lambda i, j: M_flat_polys[i].derivative(vars[j]),
    )

    return {
        "ring": P,
        "M": M,
        "M_flat": M_flat,
        "M_flat_polys": M_flat_polys,
        "jacobian_flat": jacobian_flat,
        "vars": vars,
        "G0": G0,
        "G1": G1,
        "Gamma0": Gamma0,
        "Gamma1": Gamma1,
        "gamma0": gamma0,
        "gamma1": gamma1,
        "L": L,
        "R": R,
    }


s = 3

for k in range(4, 7):
    d = build_M_old(k, s, base_field=F)
    P = d["ring"]
    g = P.gens_dict()
    M = d["M"]
    J = d["jacobian_flat"]
    sample = {var: F.random_element() for var in P.gens()}
    print("old embedding, k = {}, s = {}, dim_OT = {}".format(k, s, J.subs(sample).rank()))


for k in range(4, 7):
    d = build_M_new(k, s, base_field=F, with_scale=True)
    P = d["ring"]
    g = P.gens_dict()
    M = d["M"]
    J = d["jacobian_flat"]
    sample = {var: F.random_element() for var in P.gens()}
    print("new embedding, k = {}, s = {}, dim_OT = {}".format(k, s, J.subs(sample).rank()))

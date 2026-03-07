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

def build_M_antisym(k, s, base_field=F):
    if k <= 0 or s <= 0:
        raise ValueError("k and s must be positive integers")

    base_vars = ["t", "a", "b", "c", "x"]
    l_vars = [f"L_{i}_{j}" for i in range(k) for j in range(4)]
    names = base_vars + l_vars

    R0 = PolynomialRing(base_field, names=names)
    P = R0.fraction_field()
    gens = P.gens_dict()

    t = gens["t"]
    a = gens["a"]
    c = gens["c"]
    b = gens["b"]      
    d = a*b/c         #to enforce ab = cd
    x = gens["x"]

    L = Matrix(P, k, 4, [gens[f"L_{i}_{j}"] for i in range(k) for j in range(4)])
    R = Matrix(P, 4, k, [gens[f"L_{j}_{i}"] for i in range(4) for j in range(k)])
    # transposed matrix to preserve antisymmetry

    o = P(0) #zero
    U = Matrix(P, 
        [[o, -x, a, c],
         [x, o, d, b ],
         [-a, -d, o, o],
         [-c, -b, o, o]])

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


s = 3
if False:
    for k in range(4, 7):
        d = build_M_antisym(k, s, base_field=F)
        P = d["ring"]
        g = P.gens_dict()
        M = d["M"]
        J = d["jacobian_flat"]
        sample = {var: F.random_element() for var in P.gens()}
        dim_OT = J.subs(sample).rank()
        assert(dim_OT == 2*k - 2)
        print("antisym embedding, k = {}, s = {}, dim_OT = {}".format(k, s, dim_OT))

## code above confirms that what we obtain has generic form antisymmetric of rk 2 \otimes toeplitz of rank 1
## now let's avoid divisions and just add up m terms of this form to check when they generate everything

def build_sum(k, s, m, base_field=F):
    if k <= 0 or s <= 0 or m <= 0:
        raise ValueError("k, s and m must be positive integers")

    t_vars = [f"t_{r}" for r in range(m)]
    a_vars = [f"a_{r}_{i}" for r in range(m) for i in range(k)]
    b_vars = [f"b_{r}_{i}" for r in range(m) for i in range(k)]
    names = t_vars + a_vars + b_vars

    P = PolynomialRing(base_field, names=names)
    gens = P.gens_dict()

    G_terms = []
    gamma_terms = []
    Gamma_terms = []

    for r in range(m):
        a = vector(P, [gens[f"a_{r}_{i}"] for i in range(k)])
        b = vector(P, [gens[f"b_{r}_{i}"] for i in range(k)])

        # a \wedge b = a^t b - b^t a (antisymmetric, rank <= 2)
        G = a.column() * b.row() - b.column() * a.row()
        gamma = [gens[f"t_{r}"] ** ell for ell in range(2 * s - 1)]
        Gamma = Matrix(P, s, s, lambda i, j: gamma[i + j])

        G_terms.append(G)
        gamma_terms.append(gamma)
        Gamma_terms.append(Gamma)

    M = sum(G_terms[r].tensor_product(Gamma_terms[r]) for r in range(m))
    M_flat = Matrix(
        P,
        k * k,
        2 * s - 1,
        lambda row, ell: sum(
            G_terms[r][row // k, row % k] * gamma_terms[r][ell] for r in range(m)
        ),
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
        "G_terms": G_terms,
        "Gamma_terms": Gamma_terms,
        "gamma_terms": gamma_terms,
    }

s = 4
k = 5

for m in range(1, 12):
    d = build_sum(k, s, m, base_field=F)
    P = d["ring"]
    g = P.gens_dict()
    M = d["M"]
    J = d["jacobian_flat"]
    sample = {var: F.random_element() for var in P.gens()}
    dim_2m_rk = J.subs(sample).rank()
    print("k = {}, s = {}, m = {}, dim: {}".format(k, s, m, dim_2m_rk))

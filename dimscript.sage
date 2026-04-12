OK = "\033[92m{}\033[0m"
ERR = "\033[91m{}\033[0m"
EQ = OK.format("==")
NEQ = ERR.format("!=")

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

# s = 3
# if False:
#     for k in range(4, 7):
#         d = build_M_antisym(k, s, base_field=F)
#         P = d["ring"]
#         g = P.gens_dict()
#         M = d["M"]
#         J = d["jacobian_flat"]
#         sample = {var: F.random_element() for var in P.gens()}
#         dim_OT = J.subs(sample).rank()
#         assert(dim_OT == 2*k - 2)
#         print("antisym embedding, k = {}, s = {}, dim_OT = {}".format(k, s, dim_OT))

## code above confirms that what we obtain has generic form antisymmetric of rk 2 \otimes toeplitz of rank 1
## now let's avoid divisions and just add up m terms of this form to check when they generate everything

def build_sum(k, s, m, n, base_field=F):
    if k <= 0 or s <= 0 or m <= 0 or n <= 0:
        raise ValueError("k, s, m and n must be positive integers")

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

    M = sum(Gamma_terms[r].tensor_product(G_terms[r]) for r in range(m))
    if n > M.nrows() or n > M.ncols():
        raise ValueError("n must be at most k*s")
    M = M.matrix_from_rows_and_columns(range(n), range(n))
    n_rows = M.nrows()
    n_cols = M.ncols()
    seen = set()
    unique_nonzero = []
    for i in range(n_rows):
        for j in range(n_cols):
            val = M[i, j]
            if val != 0 and val not in seen and (-val) not in seen:
                seen.add(val)
                unique_nonzero.append(val)
    M_flat = Matrix(P, len(unique_nonzero), 1, unique_nonzero)

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

def build_sum_rank_2_antisymm(k, s, m, n, base_field=F):
    if k <= 0 or s <= 0 or m <= 0 or n <= 0:
        raise ValueError("k, s, m and n must be positive integers")

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

    M = sum(Gamma_terms[r].tensor_product(G_terms[r]) for r in range(m))
    if n > M.nrows() or n > M.ncols():
        raise ValueError("n must be at most k*s")

    s = set()
    for el in (M[i, j] for i in range(n) for j in range(n)):
        if not -el in s:
            s.add(el)
    s.discard(0)

    M_flat_polys = list(s)
    
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
        "M_flat_polys": M_flat_polys,
        "jacobian_flat": jacobian_flat,
        "vars": vars,
        "G_terms": G_terms,
        "Gamma_terms": Gamma_terms,
        "gamma_terms": gamma_terms,
    }

    


import itertools as iter

for (s, k) in sorted(iter.product(range(2, 6), range(2, 6)), key=lambda x: x[0] * x[1]):
    for n in range((k * s - k + 1) >> 1 << 1, k * s + 1, 2):
        if n <= k: 
            continue
        print("--")
        if n > 16:
            continue
        m_upper_bound = (n + 3) // 2
        m_lower_bound = 1 #  max(m_upper_bound - k - 2, 1)
        for m in range(m_lower_bound, m_upper_bound + 1):
            d = build_sum_rank_2_antisymm(k, s, m, n, base_field=F)
            P = d["ring"]
            g = P.gens_dict()
            M = d["M"]
            J = d["jacobian_flat"]
            sample = {var: F.random_element() for var in P.gens()}
            dim_2m_rk = J.subs(sample).rank()
            full = len(d["M_flat_polys"])
            
            expected_full = (k - 1) * n - (k - 1) * k / 2
           
            expected_dim = (2 * k - 3 + 1) * m - max(2 * m + k - n, 0) * (2 * m + k - n - 1) / 2 + max(0, 2 * m - n) * (2 * m - n + 1) / 2 
            
            efeq = EQ if expected_full == full else NEQ
            edeq = EQ if expected_dim == dim_2m_rk else NEQ                    

            print(f"at: k={k} s={s} n={n: ^2} m={m: ^2} : full={str(full).rjust(2, ' ')} {efeq} {str(expected_full).ljust(2, ' ')}=expected_full, dim={str(dim_2m_rk).rjust(2, ' ')} {edeq} {str(expected_dim).ljust(2, ' ')}=expected_dim")



# def build_sum_rank_1(k, s, m, n, base_field=F):
#     if k <= 0 or s <= 0 or m <= 0 or n <= 0:
#         raise ValueError("k, s, m and n must be positive integers")
# 
#     t_vars = [f"t_{r}" for r in range(m)]
#     a_vars = [f"a_{r}_{i}" for r in range(m) for i in range(k)]
#     b_vars = [f"b_{r}_{i}" for r in range(m) for i in range(k)]
#     names = t_vars + a_vars + b_vars
# 
#     P = PolynomialRing(base_field, names=names)
#     gens = P.gens_dict()
# 
#     G_terms = []
#     gamma_terms = []
#     Gamma_terms = []
# 
#     for r in range(m):
#         a = vector(P, [gens[f"a_{r}_{i}"] for i in range(k)])
#         b = vector(P, [gens[f"b_{r}_{i}"] for i in range(k)])
# 
#         # a \wedge b = a^t b - b^t a (antisymmetric, rank <= 2)
#         G = a.column() * b.row()
#         gamma = [gens[f"t_{r}"] ** ell for ell in range(2 * s - 1)]
#         Gamma = Matrix(P, s, s, lambda i, j: gamma[i + j])
# 
#         G_terms.append(G)
#         gamma_terms.append(gamma)
#         Gamma_terms.append(Gamma)
# 
#     M = sum(Gamma_terms[r].tensor_product(G_terms[r]) for r in range(m))
#     if n > M.nrows() or n > M.ncols():
#         raise ValueError("n must be at most k*s")
# 
#     M_flat_polys = list(set(M[i, j] for i in range(n) for j in range(n)))
#     
#     vars = P.gens()
#     jacobian_flat = Matrix(
#         P,
#         len(M_flat_polys),
#         len(vars),
#         lambda i, j: M_flat_polys[i].derivative(vars[j]),
#     )
# 
#     return {
#         "ring": P,
#         "M": M,
#         "M_flat_polys": M_flat_polys,
#         "jacobian_flat": jacobian_flat,
#         "vars": vars,
#         "G_terms": G_terms,
#         "Gamma_terms": Gamma_terms,
#         "gamma_terms": gamma_terms,
#     }
# 
#     
# 
# 
# for s in range(2, 5):
#     print("--")
#     for k in range(2, 5):
#         print("--")
#         for n in range(k * s - k + 1, k * s + 1):
#             if n > 20:
#                 continue
#             m_upper_bound = n
#             m_lower_bound = max(m_upper_bound - k - 2, 1)
#             for m in range(m_lower_bound, m_upper_bound + 1):
#                 d = build_sum_rank_1(k, s, m, n, base_field=F)
#                 P = d["ring"]
#                 g = P.gens_dict()
#                 M = d["M"]
#                 J = d["jacobian_flat"]
#                 sample = {var: F.random_element() for var in P.gens()}
#                 dim_2m_rk = J.subs(sample).rank()
#                 full = len(d["M_flat_polys"])
#                 
#                 expected_dim = 2 * m * k - (max(k + m - n, 0) ** 2)
#                 
#                 if expected_dim != dim_2m_rk: 
#                     print(f"Failed at: k={k}, s={s}, n={n}, m={m}, full={full}, dim={dim_2m_rk} != expected_dim={expected_dim}")
#                 else: 
#                     print(f"OK     at: k={k}, s={s}, n={n}, m={m}, full={full}, dim={dim_2m_rk},   expected_dim={expected_dim}")
# 


# for s in range(1, 6):
#     # print("Enter s = {}".format(s))
#     for k in range(2, 6):
#         m_upper_bound = k * s // 2
#         m_lower_bound = max(m_upper_bound - k, 1)
#         for m in range(1, m_upper_bound):
#             if m > 14:
#                 continue
#             n = 2 * m + 2
#             d = build_sum(k, s, m, n, base_field=F)
#             P = d["ring"]
#             g = P.gens_dict()
#             M = d["M"]
#             J = d["jacobian_flat"]
#             sample = {var: F.random_element() for var in P.gens()}
#             dim_2m_rk = J.subs(sample).rank()
#             full = len(d["M_flat_polys"])
#             if dim_2m_rk + 1 != full:
#                 print("Found error at: k = {}, s = {}, m = {}, dim: {}, full: {}".format(k, s, m, dim_2m_rk, full))



# def build_M_block_rank(P, t, block_rank, k):
#     gens = P.gens_dict()
#     L = Matrix(P, k, block_rank, [
#             gens[f"L_{t}_{i - block_rank}_{j}"] if i + 1 > block_rank else (F(1) if i == j else F(0)) 
#         for i in range(k) for j in range(block_rank)
#     ])
#     R = Matrix(P, block_rank, k, [gens[f"R_{t}_{i}_{j}"] for i in range(block_rank) for j in range(k)])
#     return L * R
# 
# def build_M_toeplitz_rank(P, t, toeplitz_rank, s):
#     gens = P.gens_dict()
#     total = Matrix(P, s, s, lambda i, j: 0)
#     for k in range(toeplitz_rank):
#         t = gens[f"t_{t}_{k}"]
#         gamma = [t**m for m in range(2 * s - 1)]
#         total += Matrix(P, s, s, lambda i, j: gamma[i + j])
# 
#     return total
# 
# def build_M_from_rank(P, t, block_rank, k, toeplitz_rank, s):
#     if k <= 0 or s <= 0:
#         raise ValueError("k and s must be positive integers")
#     return build_M_toeplitz_rank(P, t, toeplitz_rank, s).tensor_product(build_M_block_rank(P, t, block_rank, k))
# 
# def build_sum_by_ranks(block_rank, k, toeplitz_rank, s, m, base_field=F):
#     names = []
#     for t in range(m):
#         names += [f"L_{t}_{i - block_rank}_{j}" for i in range(block_rank, k) for j in range(block_rank)]
#         names += [f"R_{t}_{i}_{j}" for i in range(block_rank) for j in range(k)]
#         names += [f"t_{t}_{i}" for i in range(toeplitz_rank)]
#     
#     R0 = PolynomialRing(base_field, names=names)
#     P = R0
#     
#     return sum([
#         build_M_from_rank(P, t, block_rank, k, toeplitz_rank, s)
#         for t in range(m)
#     ]), P
# 
# def block_toeplitz_matrix_eqs(M, k, s, r, elements):
#     if elements.dimensions() != (k, 2 * r):
#         raise Exception(f"elements should have shape (k, 2r), current: {elements.shape}, expected: ({k}, {r})")
#     size = k * s
#     
#     eqs = []
#     for i in range(k):
#         for j in range(r + 1):
#             eqs.append(M[i, size - 1 - i - r + j] - elements[i, j])
#         for j in range(r + 1, 2 * r):
#             eqs.append(M[i + j - r, size - 1 - i] - elements[i, j])
#     return eqs
# 
# k = 3
# s = 3
# r = 2
# M, P = build_sum_by_ranks(1, k, 1, s, r)
# 
# EQS = block_toeplitz_matrix_eqs(M, k, s, r, Matrix(F, k, 2 * r, lambda i, j: F.random_element()))
# ideal = Ideal(*EQS)
# print(ideal.dimension())
# print(ideal.variety())

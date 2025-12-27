# This (very simple) script computes costs of Miller loop and final exponent under various assumptions.
# Assumptions / choices that we have:
# 1. Permutation program is safe (bool)
# 2. Allowed batching level
# Basically, multiplication gates can be batching using RLC.
# batching level 0:
# no batching allowed
# batching level 1:
# allowed to batch ``a * b = c``, ``a * d = e`` as ``a * (b + td) = c + te``
# this doesn't introduce new witness
# [UNSUPPORTED] batching level 2:
# allowed to batch arbitrary multiplications using cross-product trick
# this doesn't look very safe but definitely would be very profitable

pseudo_binary_encoding = [0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0,
0, 0, 1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, 1, 1,
]
Q = 21888242871839275222246405745257275088696311157297823662689037894645226208583
R = 21888242871839275222246405745257275088548364400416034343698204186575808495617
E = (Q ** 12 - 1) // R

class Cfg:
    ### lookup_size parameter is bitsize of the lookup
    ### 1 means no lookup allowed
    def __init__(self, lookup_size, batching_level):
        self.lookup_size = lookup_size
        self.batching_level = batching_level
        self.n_wtns = 0
        self.n_rk = 0

    ### cost of the rangecheck of an already allocated element
    def pay_range_check(self, n):
        # allocate new witness in a permutation
        if self.lookup_size > 1:
            self.n_rk += n
            self.n_wtns += n
        else:
            # it is a single bit, so we just do normal all-accept ADP that has rank 1
            self.n_rk += n

    ### non-optimized multiplication
    def pay_mul(self, n):
        self.n_rk += 2 * n
        self.n_wtns += n

    def pay_batch_mul(self, width, n):
        if self.batching_level == 0:
            self.pay_mul(width * n)
        elif self.batching_level == 1:
            self.n_rk += 2 * n
            self.n_wtns += width * n
        else:
            assert False

    ### validate that value is in the base field
    def pay_base_field_check(self, n):
        n_limbs = div_ceil(254, self.lookup_size)
        self.n_wtns += (n_limbs - 1) * n
        self.pay_range_check(n_limbs * n)

    def pay_is_on_curve(self, n):
        # x2 = x * x
        # x3 = x * x2
        # these two are batchable
        # single constraint to check (x3 - B) = y * y 
        self.pay_batch_mul(2, n)
        self.pay_mul(n)
        self.n_wtns -= n # as we do not actually pay for result, just constrain


    ### cost of the initialization
    def pay_prelude(self):
        # lookups cost (num_lookuped_element + table_size) rank
        # and (2 * num_lookuped_element + table_size) witness
        # so in the prelude we pay 1 table size in both witness and rank
        if self.lookup_size > 1:
            self.n_rk += 2 ** self.lookup_size
            self.n_wtns += 2 ** self.lookup_size
        # need to normalize inputs
        self.pay_base_field_check(8) # 2 G1 points and 1 G2 point
        # find element x satisfying defining equation of Fq^12
        self.n_wtns += 12
        self.n_rk += 24 # approx, actually it is maybe 28 when we use batching_level = 0?
        self.pay_is_on_curve(3)

    ### doubles the point and simultaneously computes linefunc in another point
    ### f <- f^2 * linefunc(T, T, P)
    ### T <- 2T
    def pay_miller_double(self):
        # slope is 3x^2 / 2y
        # x_2 = x * x
        # s * y = 3/2 x_2
        # y_inv * y = 1
        # last 2 are batched
        self.pay_mul(1)
        self.pay_batch_mul(2, 1) # 2 new witnesses introduced which is accidentally correct
        # introduce new point
        self.n_wtns += 2
        self.pay_is_on_curve(1)
        # evaluate (y_2T - y_T) - s * (x_2T - x_T)
        # simultaneously check that they do not coincide (x_2T - x_T) * x_inv = 1
        # once again 2 new witnesses introduced - eval and x_inv
        self.pay_batch_mul(2, 1)
        # f2 = f * linefunc
        # f_new = f * f2
        # can be batched, lol
        self.pay_batch_mul(2, 1)

    ### f <- f * linefunc(T, Q, P), T <- T + Q
    def pay_miller_add(self):
        # slope is (y2-y1)/(x2 - x1)
        # simultaneously check that x2-x1 != 0,
        # which also validates that points are different
        self.pay_batch_mul(2, 1)
        # introduce new point
        self.n_wtns += 2
        self.pay_is_on_curve(1)
        # evaluate slope (y3 - y1) - s * (x3 - x1)
        # simultaneously evaluate that x3 - x1 is nonzero
        self.pay_batch_mul(2, 1)
        # validate also that x2 - x1 is nonzero
        self.pay_mul(1)
        # multiply f by linefunc output
        self.pay_mul(1)

    def pay_miller_loop(self):
        for val in pseudo_binary_encoding[63::-1]:
            self.pay_miller_double()
            if val != 0:
                self.pay_miller_add()
        # frobenius twists of Q
        self.pay_base_field_check(4) # to parse two coordinates into G2 point
        # after that, Frobenius is linear
        # multiply 4 elements by 12 coeffs of 1st frobenius and 12 coeffs of 2nd one
        self.pay_batch_mul(24, 4)
        # actually this definitely can be double-batched using cross-products...
        self.pay_miller_add()
        self.pay_miller_add()

    ### pay exponentiation in Fq12
    ### method - compute all x^2^k, then multiply ones that
    ### we need by using running product 
    def pay_exp(self, e):
        n_ones = 0
        log_e = 0
        while e > 0:
            n_ones += e % 2 # count bits equal to 1
            log_e += 1
            e //= 2
        self.pay_mul(log_e - n_ones)
        self.pay_batch_mul(2, n_ones)

    def ciphertext_n_bytes(self):
        return ((self.n_rk + 1) ** 2 * self.n_wtns * 32)

def div_ceil(a, b):
    return (a + b - 1) // b

cfg = Cfg(lookup_size=1, batching_level=0)
cfg.pay_prelude()
cfg.pay_miller_loop()
cfg.pay_miller_loop()
cfg.pay_exp(E)

print("{x} TB".format(x = cfg.ciphertext_n_bytes() // 2**40))
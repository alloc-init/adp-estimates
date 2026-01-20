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
    def __init__(self):
        self.n_wtns = 0
        self.n_gates = 0

    ### c = sqrt(ab) gate
    def mul_sq(self):
        self.n_gates += 1
        self.n_wtns += 1

    # Extension finding gadget, outputs s_i basis of the extension
    def extension_basis_finder(self):
        self.n_gates += 22
        self.n_wtns += 11

    # Checks that value is in the base field. Both bitcheck and frobcheck cost the same
    def basefield_check(self):
        self.n_wtns += 253
        self.n_gates += 254

    # Gadget that takes G1 point (x, y) and outputs u_i = x s_i, v_i = y s_i
    # Performs Frobcheck of x, y and also checks that (x, y) is on curve
    def miller_bases(self):
        self.basefield_check() # x
        self.basefield_check() # y
        self.n_gates += 1 # y * y = (x_3 + b) * 1
        self.n_wtns += 2 * 22 # introduce all u_i, v_i for i > 1
        self.n_gates += 2 * 44 # constrain them + certify them using the squaring trick
    
    # Outputs a twist by element of mu_{32}, which is a cofactor of the odd-order subgroup of (Fq^{12})^*
    # ---
    tmp = Q ** 12 - 1
    assert(tmp % 32 == 0)
    assert(tmp % 64 != 0)
    # ---
    def twist(self):
        self.n_wtns += 5 # t, t2, t4, t8, t16 
        self.n_gates += 5 # tt = t2, t2t2 = t4, t4t4 = t8, t8t8 = t16, t16t16 = 1

    def miller_iter(self, pseudo_bit):
        self.twist()
        self.mul_sq()
        if pseudo_bit != 0:
            self.twist()
            self.mul_sq()
        self.mul_sq()
        self.mul_sq()
    
    def miller_loop(self):
        self.miller_iter(1)
        for v in pseudo_binary_encoding[63::-1]:
            self.miller_iter(v)

    def odd_exp_fixpow(self, N):
        n_bits = 0
        while N > 0:
            N >>= 1
            n_bits += 1
        for i in range(n_bits):
            self.mul_sq()

    def odd_exp_varpow(self, bitwidth):
        for i in range(bitwidth):
            self.mul_sq()
            self.mul_sq()

    # This skips finding extension basis as we will do it in the main phase
    def verify_pairing_eq(self):
        self.miller_bases()
        self.miller_bases()
        self.miller_loop()
        self.miller_loop()
        self.miller_loop()
        self.odd_exp_varpow(130*7)
        self.mul_sq()
        self.mul_sq()
        self.mul_sq()
        self.odd_exp_fixpow(E)        
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.mul_sq()
        self.mul_sq()
        self.mul_sq()
        self.mul_sq()
        self.n_gates += 1

    # compute 5-th power
    def sbox(self):
        self.n_wtns += 3
        self.n_gates += 5

    # width 4 permutation, certify *dropped* outputs
    def poseidon(self):
        width = 4
        n_full_rounds = 8
        n_part_rounds = 54 # approx, need to check the script for Fq specifically
        for i in range(n_full_rounds * width + n_part_rounds):
            self.sbox()
        for i in range(width - 1):
            self.basefield_check()
    
    def verify_arith(self):
        self.n_gates += 1
        self.n_wtns += 254
        self.n_gates += 254 # z^2 fully fits

        self.n_gates += 4
        self.n_wtns += 130 * 5 # computing z^3 in a non-reduced form \sum a_i 2^{64 i}
        self.n_gates += 130 * 5 # 5 limbs (because z is 2 limb and z^2 is 4 limb), each limb is 128 + < 10

        self.n_gates += 4
        self.n_wtns += 130 * 7
        self.n_gates += 130 * 7 # computing z^4

        # multiplying p_i' by phi^i
        # p_0'
        self.odd_exp_varpow(254)
        # p_1'
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(126)
        # p_2'
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        # p_3'
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(130 * 5)

        for i in range(4):
            self.odd_exp_varpow(254) # v_i psi^i
            self.basefield_check() # check v_i
        for i in range(6):
            self.mul_sq()

        # checking v1 v2 - v3 = z (v0 vZ + u + z)
        
        # v1v2 - v3
        self.n_gates += 4
        self.n_wtns += 130 * 7
        self.n_gates += 130 * 7

        # v0vZ + u + z
        self.n_gates += 4
        self.n_wtns += 130 * 7
        self.n_gates += 130 * 7

        # varpows
        self.odd_exp_varpow(130 * 7)
        self.odd_exp_varpow(130 * 7)
        self.odd_exp_varpow(128)

    def verify_full(self):
        (w, g) = (self.n_wtns, self.n_gates)
        self.extension_basis_finder()
        print("Extension:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
        self.basefield_check() #check u
        for i in range(4):
            self.basefield_check() # check p_i
        print("Glue:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
 
        self.poseidon() # z = Hash(pp, u, P.x, P.y)
        print("Hash:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
        
        self.basefield_check() # check z
        self.verify_arith()
        print("Arithmetic:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
 
        self.verify_pairing_eq() # P, Q checked by miller_bases routine
        print("Pairing equation:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
 

    def ciphertext_n_bytes(self):
        return ((2 * self.n_gates + 1) ** 2 * self.n_wtns * 32)

def div_ceil(a, b):
    return (a + b - 1) // b

cfg = Cfg()
cfg.verify_full()
print("{x} TB".format(x = cfg.ciphertext_n_bytes() // 2**40))
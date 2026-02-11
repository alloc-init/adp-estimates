import math

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
        self.odd_exp_varpow(128) # multiplying by z
        self.mul_sq()
        self.mul_sq()
        self.mul_sq()
        self.odd_exp_fixpow(E)        
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
        self.odd_exp_varpow(254)
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
        self.basefield_check() # obtain z via truncation of poseidon output
        # checking vA^2 - vC = vH vZ + u
        # in the exponent        
        self.odd_exp_varpow(256 * 2)
        self.odd_exp_varpow(256)
        self.odd_exp_varpow(256 * 2)
        self.odd_exp_varpow(256)


    def verify_full(self):
        (w, g) = (self.n_wtns, self.n_gates)
        self.extension_basis_finder()
        print("Extension:")
        print("wtns = {}, gates = {}".format(self.n_wtns - w, self.n_gates - g))
        (w, g) = (self.n_wtns, self.n_gates)
        self.basefield_check() #check u
        for i in range(4):
            self.basefield_check() # check v_A, v_C, v_H, v_Z
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

cfg = Cfg()
cfg.verify_full()
print("wtns: {}, n_gates: {}".format(cfg.n_wtns, cfg.n_gates))
print("Ciphertext size: {x} TB".format(x = cfg.ciphertext_n_bytes() // 2**40))
print()


print("* base costs")
# Benchmarked bn254 fp12 multiplication using the mcl library on M3 MacBook Air
ns_per_field_mul = 700 
print(f"- bn254 fp12 mul cost = {ns_per_field_mul}ns (mcl library on m3 macbook air)")

# Since witness is in fp12 and ciphertext is in fp, we don't have to compute full fp12 multiplication to compute E.
# It should be significantly cheaper to compute fp x fp12 than fp12 x fp12.
# - ChatGPT thinks it is 4.5x cheaper.
ns_per_scalar_mul = ns_per_field_mul / 4.5
print(f"- bn254 fp x fp12 mul cost ~= {ns_per_scalar_mul:.2f}ns (4.5x cheaper than fp12 mul)")

ns_per_inv = 2500
print(f"- bn254 fp12 inv cost = {ns_per_inv}ns (mcl library on m3 macbook air)")

print()

##########################################
# Decryption main costs                  #
##########################################
# 1. compute E <- Mhat(1|x)
# 2. solve for t: det(E - t * x0 ...) = 0
print("* serial costs")

# In order to compute E, you have to compute (2*n+1)^2 inner products with the witness, which is exactly the size of M.
M_size = (2 * cfg.n_gates + 1) ** 2 * cfg.n_wtns
print(f"- number of scalar muls to compute E = {M_size}")

mins_to_compute_e = int(M_size * ns_per_scalar_mul / 1_000_000_000 / 60)
print(f"- serial time to compute E = {mins_to_compute_e} min")

# Determinant computed with gaussian elimination uses n^3/3 multiplications.
det_num_field_mults = (2*cfg.n_gates+1) ** 3 / 3
print(f"- number of fp12 mults for computing determinant = n^3 / 3 = {det_num_field_mults}")
# Determinant computed with gaussian elimination uses n(n-1)/2 inversions.
det_num_invs = cfg.n_gates * (cfg.n_gates-1) / 2
print(f"- number of fp12 inversions for computing determinant = n(n-1)/2 = {det_num_invs}")
det_total_min = int((det_num_field_mults * ns_per_field_mul + det_num_invs * ns_per_inv) / 1_000_000_000 / 60)
print(f"- serial time to compute determinant = {det_total_min} min")

print()

print("* cost to compute parallel decryption in 10m")

# * A Parallel Algorithm for Calculation of Large Determinants with High Accuracy for GPUs and MPI clusters
# - reports no advantages for gpu
# - large cpu cluster speedup is linear with number of cores
# - https://arxiv.org/abs/1308.1536

print("- determinant compute speedup linear with number of cores")
print("- computing E is embarassingly parallel")
print(f"- how many cores does it take to compute both E and det in 10m?")
print(f"    - E_minutes / n + det_minutes / n = 10")
ncores = int((det_total_min + mins_to_compute_e)/10)
print(f"    - n = {ncores}")
print(f"- how much storage does each core require?")
print(f"    - ciphertext_size / n = {cfg.ciphertext_n_bytes() / ncores / 2**30:.2f} gb")
print(f"- how many 500gb/s ciphertext storage locations do you need?")
print(f"    - since E is so parallel, you can stream M to the shards, you have as much time to send the ciphertext as it takes to compute E")
parallel_mins_to_compute_e = mins_to_compute_e/ncores
total_shards = ncores / 256
print(f"    - assume each shard has 256 cpu cores, so you have {total_shards} total shards")
print(f"    - with n={ncores}, each shard computing E takes {parallel_mins_to_compute_e:.2f} min")
print(f"    - assume 500gb/s hpc-level network")
print(f"    - total ciphertext size is {cfg.ciphertext_n_bytes() // 2**30} GB")
print(f"    - so you need {(cfg.ciphertext_n_bytes() // 2**30)} / (n_data_locs * 500) = {parallel_mins_to_compute_e:.2f} * 60")
print(f"        - {(cfg.ciphertext_n_bytes() // 2**30)} / ({parallel_mins_to_compute_e:.2f} * 60 * 500) = n_data_locs")
n_data_locs = math.ceil((cfg.ciphertext_n_bytes() // 2**30) / (parallel_mins_to_compute_e * 60 * 500))
print(f"        -            = {n_data_locs} separate locations")
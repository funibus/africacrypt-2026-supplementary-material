from dilithium_py.ml_dsa.ml_dsa import ML_DSA
import random as rd


class MLDSARK(ML_DSA):

    def __init__(self, parameter_set):
        super().__init__(parameter_set)
        
        self.beta = 2 * self.tau * self.eta
        self.len_th = 256*3*self.k//8

    def _sk_size(self) -> int:
        if self.eta == 2:
            s_bytes = 128
        else:
            s_bytes = 256
        s1_len = s_bytes * self.l
        s2_len = s_bytes * self.k
        t0_len = 416 * self.k
        return 2 * 32 + 64 + s1_len + s2_len + t0_len + self.len_th
    
    def _sig_size(self) -> int:
        if self.gamma_1 == 131072:
            s_bytes = 18
        else:
            assert self.gamma_1 == 524288
            s_bytes = 20
        return self.c_tilde_bytes + self.l * 32 * s_bytes + self.omega + self.k + self.len_th

    def _expand_vector_from_seed(self, rho_prime):
        s1_elements = [
            self.R.rejection_bounded_poly(rho_prime, i, self.eta) for i in range(self.l)
        ]
        s2_elements = [
            self.R.rejection_bounded_poly(rho_prime, i, self.eta)
            for i in range(self.l, self.l + self.k)
        ]

        s1 = self.M.vector(s1_elements)
        s2 = self.M.vector(s2_elements)
        return s1, s2
  
    def _unpack_hint_t1(self, packed_bytes):
        matrix = self._hint_decode(packed_bytes, self.k, self.R.n)
        matrix = [self.R(r) for r in matrix]
        return self.M(matrix).transpose()

    def _encode_hint(self, th):
        packed_bytes = bytearray()

        # Function to encode each coefficient as 3 bits
        def encode_coeff(c):
            if (c % self.R.q) == 0:
                return 0b000  # 0 -> 000
            elif (c - 1) % self.R.q == 0:
                return 0b001  # 1 -> 001
            elif (c + 1) % self.R.q == 0:
                return 0b010  # -1 -> 010
            elif (c + 1022) % self.R.q == 0:
                return 0b011  # -1022 -> 011
            elif (c + 1023) % self.R.q == 0:
                return 0b100  # -1023 -> 100
            elif (c - 1022) % self.R.q == 0:
                return 0b101  # 1022 -> 101    
            else:
                raise ValueError(f"Coefficient must be in [0, 1, -1, -1022, -1023, 1022], got {c}")

        # Loop through rows 
        rows = [th[i, 0] for i in range(th.dim()[0])]
        for poly in rows:
            coeff_bits = []  # Collect encoded coefficients

            # Encode each coefficient as 3 bits
            for coeff in poly.coeffs:
                coeff_bits.append(encode_coeff(coeff))

            # Pack the coefficients into bytes
            # Each byte can hold up to 8 coefficients packed as 3 bits (24 bits per 3 bytes)
            bit_buffer = 0  # Temporary storage for bits
            bit_count = 0   # Count of bits in the buffer

            for coeff_bit in coeff_bits:
                bit_buffer = (bit_buffer << 3) | coeff_bit  # Add 3 bits to the buffer
                bit_count += 3

                # Once we have at least 8 bits, pack them into bytes
                while bit_count >= 8:
                    packed_bytes.append((bit_buffer >> (bit_count - 8)) & 0xFF)  # Extract top 8 bits
                    bit_count -= 8

            # If any bits remain in the buffer, pack them into the last byte (padded with zeros)
            if bit_count > 0:
                packed_bytes.append((bit_buffer << (8 - bit_count)) & 0xFF)

        return bytes(packed_bytes)
    
    def _hint_decode(self, packed_bytes, num_rows, num_coeffs_per_poly):
        decoded_polys = []

        # Function to decode each 3-bit value into a coefficient
        def decode_coeff(encoded):
            # q1 is (q-1)/alpha, which is 1023 in our case
            if encoded == 0b000:
                return 0  # 000 -> 0
            elif encoded == 0b001:
                return 1  # 001 -> 1
            elif encoded == 0b010:
                return -1  # 010 -> -1
            elif encoded == 0b011:
                return -1022  # 011 -> -1022
            elif encoded == 0b100:
                return -1023  # 100 -> -1023
            elif encoded == 0b101:
                return 1022  # 101 -> 1022
            else:
                raise ValueError(f"Invalid encoded value: {encoded}")

        # Iterate over each row in the packed bytes
        bit_buffer = 0    # Temporary storage for the bits
        bit_count = 0     # Number of bits in the buffer
        byte_index = 0    # Current byte being processed

        for _ in range(num_rows):
            poly_coeffs = []

            # We have num_coeffs_per_poly coefficients per polynomial
            for _ in range(num_coeffs_per_poly):
                # Refill the bit buffer if there are not enough bits
                while bit_count < 3:
                    bit_buffer = (bit_buffer << 8) | packed_bytes[byte_index]
                    bit_count += 8
                    byte_index += 1

                # Extract the top 3 bits (next coefficient) from the bit buffer
                encoded_coeff = (bit_buffer >> (bit_count - 3)) & 0b111
                poly_coeffs.append(decode_coeff(encoded_coeff))

                # Remove the 3 bits from the bit buffer
                bit_count -= 3

            # Append the decoded polynomial to the result
            decoded_polys.append(poly_coeffs[:num_coeffs_per_poly])  # Ensure correct slicing

        return decoded_polys
    
    def _pack_sk(self, rho, K, tr, s1, s2, t0, th):
        s1_bytes = s1.bit_pack_s(2*self.eta)
        s2_bytes = s2.bit_pack_s(2*self.eta)
        t0_bytes = t0.bit_pack_t0()
        th_bytes = self._encode_hint(th)
        return rho + K + tr + s1_bytes + s2_bytes + t0_bytes + th_bytes
    
    def _unpack_sk(self, sk):
        if self.eta == 2:
            s_bytes = 128
        else:
            s_bytes = 96
        s1_len = s_bytes * self.l
        s2_len = s_bytes * self.k
        t0_len = 416 * self.k
        if len(sk) != self._sk_size():
            raise ValueError("sk packed bytes is of the wrong length")

        # Split bytes between seeds and vectors
        sk_seed_bytes, sk_vec_bytes = sk[:128], sk[128:]

        # Unpack seed bytes
        rho, K, tr = (
            sk_seed_bytes[:32],
            sk_seed_bytes[32:64],
            sk_seed_bytes[64:128],
        )

        # Unpack vector bytes
        s1_bytes = sk_vec_bytes[:s1_len]
        s2_bytes = sk_vec_bytes[s1_len : s1_len + s2_len]
        t0_bytes = sk_vec_bytes[s1_len + s2_len:s1_len + s2_len + t0_len]

        # Unpack bytes to vectors 
        s1 = self.M.bit_unpack_s(s1_bytes, self.l, 2*self.eta)
        s2 = self.M.bit_unpack_s(s2_bytes, self.k, 2*self.eta)
        t0 = self.M.bit_unpack_t0(t0_bytes, self.k)

        # Add tH
        th = sk[-self.len_th:]
        return (rho, K, tr, s1, s2, t0, th)
    
    def _unpack_h(self, h_bytes):
        h_len = len(h_bytes) - self.len_th
        h = super()._unpack_h(h_bytes[:h_len])
        th = self._unpack_hint_t1(h_bytes[-self.len_th:])
        return h, th
    
    def _unpack_sig(self, sig_bytes):
        c_tilde = sig_bytes[:32]
        z_bytes = sig_bytes[32 : -(self.len_th +self.k + self.omega)]
        h_bytes = sig_bytes[-(self.len_th +self.k + self.omega) :]

        z = self.M.bit_unpack_z(z_bytes, self.l, self.gamma_1)
        h, th = self._unpack_h(h_bytes)
        return c_tilde, z, (h, th)
    
    def _pack_sig(self, c_tilde, z, h):
        h, th = h
        h_pack = super()._pack_h(h)
        return c_tilde + z.bit_pack_z(self.gamma_1) + h_pack + th

    def _get_null_hint_t1(self):
        matrix = [[self.R([0]*256) for _ in range(self.k)]]
        return self.M(matrix).T
   
    def _keygen_internal(self, zeta: bytes) -> tuple[bytes, bytes]:
        """
        Generates a public-private key pair from a seed following
        Algorithm 6 (FIPS 204)
        """
        # Expand with an XOF (SHAKE256)
        seed_domain_sep = zeta + bytes([self.k]) + bytes([self.l])
        seed_bytes = self._h(seed_domain_sep, 128)

        # Split bytes into suitable chunks
        rho, rho_prime, K = seed_bytes[:32], seed_bytes[32:96], seed_bytes[96:]

        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Generate the error vectors s1 ∈ R^l, s2 ∈ R^k
        s1, s2 = self._expand_vector_from_seed(rho_prime)
        s1_hat = s1.to_ntt()

        # Matrix multiplication
        t = (A_hat @ s1_hat).from_ntt() + s2

        t1, t0 = t.power_2_round(self.d)
        th = self._get_null_hint_t1()
        # Pack up the bytes
        pk = self._pack_pk(rho, t1)
        tr = self._h(pk, 64)

        sk = self._pack_sk(rho, K, tr, s1, s2, t0, th)
        return pk, sk

    def rand_sk(self, sk, mu):
        """
        Update a secrete key sk
        """
        rho, K,_, s1, s2, t0, _ = self._unpack_sk(sk)
        rho_prime = self._h(rho+mu, 128)
        
        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Generate the error vectors s1 ∈ R^l, s2 ∈ R^k
        s1_prime, s2_prime = self._expand_vector_from_seed(rho_prime)
        s1_sum = s1+s1_prime
        s2_sum = s2+s2_prime
        # Matrix multiplication
        t = (A_hat @ s1.to_ntt()).from_ntt() + s2
        r = (A_hat @ s1_prime.to_ntt()).from_ntt() + s2_prime
        t_sum = t + r
        
        t1_sum, t0_sum = t_sum.power_2_round(self.d)
        tr_prime = self._h(self._pack_pk(rho, t1_sum), 64)
        # compute th
        t1, _ = t.power_2_round(self.d)
        r1, _ = r.power_2_round(self.d)
        th = t1_sum - (t1.scale(1<<self.d) + r1.scale(1<<self.d)).power_2_round(self.d)[0]
        sk_prime = self._pack_sk(rho, K, tr_prime, s1_sum, s2_sum, t0_sum, th)
        return sk_prime
    
    def rand_pk(self, pk, mu):
        """
        Update a public master key pk
        """
        rho, t1 = self._unpack_pk(pk)
        rho_prime = self._h(rho+mu, 128)
        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Generate the error vectors s1 ∈ R^l, s2 ∈ R^k
        s1_prime, s2_prime = self._expand_vector_from_seed(rho_prime)
        s1_prime_hat = s1_prime.to_ntt()

        # Matrix multiplication
        r = (A_hat @ s1_prime_hat).from_ntt() + s2_prime
        r1, _ = r.power_2_round(self.d)
        t1_sum = (t1.scale(1<<self.d) + r1.scale(1<<self.d)).power_2_round(self.d)[0]
        pk_prime = self._pack_pk(rho, t1_sum)
        return pk_prime
    
    def _sign_internal(
        self, sk: bytes, m: bytes, rnd: bytes, external_mu: bool = False
    ) -> bytes:
        """
        Deterministic algorithm to generate a signature for a formatted message
        M' following Algorithm 7 (FIPS 204)

        When `external_mu` is `True`, the message `m` is interpreted instead as
        the pre-hashed message `mu = prehash_external_mu()`
        """
        # unpack the secret key
        rho, k, tr, s1, s2, t0, th = self._unpack_sk(sk)

        # Precompute NTT representation
        s1_hat = s1.to_ntt()
        s2_hat = s2.to_ntt()
        t0_hat = t0.to_ntt()

        # Generate matrix A ∈ R^(kxl) in the NTT domain
        A_hat = self._expand_matrix_from_seed(rho)

        # Set seeds and nonce (kappa)
        if external_mu:
            mu = m
        else:
            mu = self._h(tr + m, 64)
        rho_prime = self._h(k + rnd + mu, 64)

        kappa = 0
        alpha = self.gamma_2 << 1
        while True:
            y = self._expand_mask_vector(rho_prime, kappa)
            y_hat = y.to_ntt()
            w = (A_hat @ y_hat).from_ntt()

            # increment the nonce
            kappa += self.l

            # NOTE: there is an optimisation possible where both the high and
            # low bits of w are extracted here, which speeds up some checks
            # below and requires the use of make_hint_optimised() -- to see the
            # implementation of this, look at the signing algorithm for
            # dilithium. We include this slower version to mirror the FIPS 204
            # document precisely.
            # Extract out only the high bits
            w1 = w.high_bits(alpha)

            # Create challenge polynomial
            w1_bytes = w1.bit_pack_w(self.gamma_2)
            c_tilde = self._h(mu + w1_bytes, self.c_tilde_bytes)
            c = self.R.sample_in_ball(c_tilde, self.tau)
            c_hat = c.to_ntt()

            # NOTE: unlike FIPS 204 we start again as soon as a vector
            # fails the norm bound to reduce any unneeded computations.
            c_s1 = s1_hat.scale(c_hat).from_ntt()
            z = y + c_s1
            if z.check_norm_bound(self.gamma_1 - self.beta):
                continue

            c_s2 = s2_hat.scale(c_hat).from_ntt()
            r0 = (w - c_s2).low_bits(alpha)
            if r0.check_norm_bound(self.gamma_2 - self.beta):
                continue

            c_t0 = t0_hat.scale(c_hat).from_ntt()
            if c_t0.check_norm_bound(self.gamma_2):
                continue

            h = (-c_t0).make_hint(w - c_s2 + c_t0, alpha)
            if h.sum_hint() > self.omega:
                continue

            h = (h, th)
            return self._pack_sig(c_tilde, z, h)
        
    def _verify_internal(self, pk: bytes, m: bytes, sig: bytes) -> bool:
        """
        Internal function to verify a signature sigma for a formatted message M'
        following Algorithm 8 (FIPS 204)
        """
        rho, t1 = self._unpack_pk(pk)
        
        try:
            c_tilde, z, (h, th) = self._unpack_sig(sig)
            t1 = t1 + th
        except ValueError:
            print("Failed to unpack signature")
            return False

        if h.sum_hint() > self.omega:
            return False

        if z.check_norm_bound(self.gamma_1 - self.beta):
            return False

        A_hat = self._expand_matrix_from_seed(rho)

        tr = self._h(pk, 64)
        mu = self._h(tr + m, 64)
        c = self.R.sample_in_ball(c_tilde, self.tau)

        # Convert to NTT for computation
        c = c.to_ntt()
        z = z.to_ntt()

        t1 = t1.scale(1 << self.d)
        t1 = t1.to_ntt()

        Az_minus_ct1 = (A_hat @ z) - t1.scale(c)
        Az_minus_ct1 = Az_minus_ct1.from_ntt()

        w_prime = h.use_hint(Az_minus_ct1, 2 * self.gamma_2)
        w_prime_bytes = w_prime.bit_pack_w(self.gamma_2)

        return c_tilde == self._h(mu + w_prime_bytes, self.c_tilde_bytes)
    

MLDSARK_PARAMS44     =  {
         "d": 13,  # number of bits dropped from t
        "tau": 39,  # number of ±1 in c
        "gamma_1": 131072,  # coefficient range of y: 2^17
        "gamma_2": 95232,  # low order rounding range: (q-1)/88
        "k": 4,  # Dimensions of A = (k, l)
        "l": 4,  # Dimensions of A = (k, l)
        "eta": 2,  # Private key range
        "omega": 80,  # Max number of ones in hint
        "c_tilde_bytes": 32,
        "oid": (2, 16, 840, 1, 101, 3, 4, 3, 17),
    }

MLDSARK44 = MLDSARK(MLDSARK_PARAMS44)


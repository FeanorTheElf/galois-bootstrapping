#pragma once
#include <vector>
#include <assert.h>
#include <ostream>
#include "seal/util/uintarithsmallmod.h"
#include "seal/seal.h"
#include "polyarith.h"

/**
 * Contains operations to work with the slot structure of
 * the plaintext space.
*/

extern size_t kswitch_counter;

 void log_relin();
 void log_galois();

inline uint64_t log2_ceil(uint64_t x) {
	if (x == 0) {
		return 0;
	}
	unsigned long result = 0;
	while (x != 0) {
		x = x >> 1;
		result += 1;
	}
	return result - 1;
};

inline uint64_t log2_exact(uint64_t x) {
	unsigned long result = log2_ceil(x);
	if (((uint64_t)1 << result) == x) {
		return result;
	}
	else {
		throw std::invalid_argument("Not a power of 2");
	}
};

inline uint64_t highest_dividing_power2(uint64_t x) {
	if (x == 0) {
		return std::numeric_limits<uint64_t>::digits;
	}
	unsigned long result = 0;
	while ((x & 1) == 0) {
		result += 1;
		x = x >> 1;
	}
	return result;
};

void test_g_automorphisms();

//tex:
//Represents the ring decomposition
//$$(\mathbb{Z}/p^e\mathbb{Z})[X]/(X^N + 1) = (\mathbb{Z}/p^e\mathbb{Z})[X]/(F_1) \times ... \times (\mathbb{Z}/p^e\mathbb{Z})[X]/(F_n)$$
//where $n$ is the number of slots.
class SlotRing {

	// the number of slots
	uint64_t n;
	// the rank of a slot (as free Z/p^eZ-module)
	uint64_t d;
	// the prime base of the characteristic; scalars are in Z/p^eZ
	uint64_t p;
	// the exponent of the characteristic; scalars are in Z/p^eZ
	uint64_t e;

	//tex:
	//first generator of the automorphism group $(\mathbb{Z}/2N\mathbb{Z})^*$ of order $N/2$
	uint64_t g1;

	//tex:
	//second generator of the automorphism group $(\mathbb{Z}/2N\mathbb{Z})^*$ of order $2$
	uint64_t g2;

	//tex:
	//whether the slot group (i.e. $(\mathbb{Z}/2N\mathbb{Z})^*/\langle p \rangle)$ is cyclic
    //(isomorphic to $\mathbb{Z}/n\mathbb{Z}$) or isomorphic to $\mathbb{Z}/2\mathbb{Z} \times \mathbb{Z}/\frac n 2 \mathbb{Z}$
	bool slot_group_cyclic;

	//tex:
	//whether the slot moduli polynomials (or equivalently the unit vector polynomials) are $d$-sparse, i.e. polynomials in $X^d$;
	//note that we require them to be at least $d/2$-sparse, but in some occasions, $d$-sparse can yield a performance gain
	bool d_sparse;

	// base-2 logarithm of 2N
	uint64_t log2N;

	//tex:
    //base-$(g_1, g_2)$ logarithm of $p$ in $(\mathbb{Z}/2N\mathbb{Z})^*$
	std::tuple<size_t, size_t> p_log;

	// modulus of the scalar ring
	seal::Modulus mod_pe;

	//tex:
    //modulus of the "index ring", i.e. the ring $\mathbb{Z}/2N\mathbb{Z}$ which provides 
	//the structure of the multiplicative group $$\langle X \rangle \subseteq ((\mathbb{Z}/p^e\mathbb{Z})[X]/(X^N + 1))^*$$
	seal::Modulus mod_2N;

	// the polynomial modulus of the whole ring, i.e. X^N + 1
	PolyModulus mod_2N_cyclotomic;

	// the minimal polynomial of a generator of one slot; currently, this matches the first slot, i.e. f mod slot_modulus is the slot 0 part of f
	PolyModulus slot_modulus;

	//tex:
    //the "unit vectors" $e_i$, i.e. elements $e_i \in (\mathbb{Z}/p^e\mathbb{Z})[X]/(X^N + 1)$ that have one slot set to 1 and the others to 0, i.e.
    //$$e_i \equiv 1 \mod F_i \quad \mathrm{and} \quad e_i \equiv 0 \mod F_j$$
    //where $j \neq i$. The $e_i$ are stored in the same order as the slots;
	//Note that if the slot group is not cyclic, this contains first the unit vectors corresponding to $\mathbb{Z}/\frac n 2\mathbb{Z}$, and then their $g_2$-automorphism
	//conjugate
	std::vector<poly> unit_vectors;

	void init_slot_group();
	void init_p_log();
	void init_slot_modulus(poly base_unit_vector);

	void apply_frobenius(poly& x, size_t iters, const seal::Modulus* mod = nullptr) const;

	//tex:
	//Computes the automorphism $X -> X^{g_1^i}$ inplace.
	//If the given modulus is nullptr, we use $p^e$ as the modulus.
	void g1_automorphism(poly& x, size_t iters, const seal::Modulus* mod = nullptr) const;

	//tex:
	//Computes the automorphism $X -> X^{g_2}$ inplace.
	//If the given modulus is nullptr, we use $p^e$ as the modulus.
	void g2_automorphism(poly& x, const seal::Modulus* mod = nullptr) const;

public:

	/**
	 * Represents a subring of the main ring. Mainly used for a single slot.
	*/
	struct SubringView {
		const seal::Modulus& scalar_mod;
		const PolyModulus& poly_mod;

		poly generator();
	};

	//tex:
	//Encapsulates a "rotation" of the slots. This refers to the map
	//$$ S^n \to S^n, \quad (x_i)_i \mapsto (x_{\lfloor i / n' \rfloor n' + i - s \mod n'})_i$$
	//In the simplest case, we have $n' = n$ and so this is $x_i \mapsto x_{i - s \mod n}$.
	class Rotation {

		//tex:
       //elements at the positions where fmask is 1 will be moved "forward" by $s$, i.e. $\text{$g_1$-auto}(\cdot, s)$
		poly fmask;

		//tex:
        //elements at the positions where bmask = 1 will be moved "backwards" by $s$, i.e. via $\text{$g_1$-auto}(\cdot, \mathrm{ord}(g_1) - s)$;
		//Note that in the case of noncyclic $(\mathbb{Z}/2N\mathbb{Z})^*/\langle p\rangle$ and maximal block_size, they will also
		//change the orbit of the action of $g_2$, i.e. we apply $\text{$g_2$-auto}(\cdot)$.
		poly bmask;
		size_t s;
		const SlotRing& slot_ring;
		size_t block_size;
		mutable std::unique_ptr<seal::Plaintext> fmask_plain;
		mutable std::unique_ptr<seal::Plaintext> bmask_plain;

		inline Rotation(poly fmask, poly bmask, size_t s, const SlotRing& slot_ring, size_t block_size)
			: fmask(std::move(fmask)), bmask(std::move(bmask)), s(s), slot_ring(slot_ring), block_size(block_size), fmask_plain(nullptr), bmask_plain(nullptr) {}

		size_t effective_block_size() const;

	public:
		Rotation(const Rotation&) = default;
		Rotation(Rotation&&) = default;
		~Rotation() = default;

		poly operator()(const poly& x) const;
		void apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const;
		std::tuple<uint32_t, uint32_t> galois_elements() const;
		bool is_identity() const;

		inline const poly& get_forward_mask() const { return fmask; }
		inline const poly& get_backward_mask() const { return bmask; }
		const std::tuple<uint64_t, uint64_t> get_forward_g1_g2_decomp() const;
		const std::tuple<uint64_t, uint64_t> get_backward_g1_g2_decomp() const;

		friend SlotRing;
	};

	/**
	 * Encapsulates the slot-wise Frobenius map
	*/
	class Frobenius {
		size_t iters;
		const SlotRing& slot_ring;

		inline Frobenius(size_t iters, const SlotRing& slot_ring)
			: iters(iters), slot_ring(slot_ring) {}

	public:
		Frobenius(const Frobenius&) = default;
		Frobenius(Frobenius&&) = default;
		~Frobenius() = default;

		poly operator()(poly x) const;
		void apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const;
		uint32_t galois_element() const;
		bool is_identity() const;

		const std::tuple<uint64_t, uint64_t> get_g1_g2_decomp() const;

		friend SlotRing;
	};

	/**
	 * Encapsulates a raw automorphism of the base ring that is induced by a
	 * Galois automorphism of the corresponding cyclotomic field.
	*/
	class RawAuto {
		//tex:
		//the automorphism will be of the form $X \mapsto X^k$ with $k \perp N$.
		//This element stores $(i, j)$ such that $k = g_1^ig_2^j$.
		std::tuple<size_t, size_t> g1_g2_decomp;
		const SlotRing& slot_ring;

		inline RawAuto(const SlotRing& slot_ring, std::tuple<size_t, size_t> g1_g2_decomp)
			: slot_ring(slot_ring), g1_g2_decomp(g1_g2_decomp) {}

	public:
		RawAuto(const RawAuto&) = default;
		RawAuto(RawAuto&&) = default;
		~RawAuto() = default;

		poly operator()(poly x) const;
		void apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const;
		uint32_t galois_element() const;
		bool is_identity() const;

		friend SlotRing;
	};

	SlotRing(
		uint64_t log2N, 
		uint64_t p, 
		uint64_t e, 
		poly base_unit_vector
	);

	SlotRing(
		uint64_t log2N,
		uint64_t p,
		uint64_t e,
		uint64_t g2,
		poly base_unit_vector
	);

	SlotRing(const SlotRing&) = default;
	SlotRing(SlotRing&&) = default;
	~SlotRing() = default;

	size_t slot_group_dims() const noexcept;
	size_t slot_group_len() const noexcept;
	size_t slot_rank() const noexcept;

	Rotation block_rotate(size_t slots, size_t block_size) const;

	/**
	 * Computes the isomorphism that moves every slot by the given number of slots.
	 * Note that relatively heavy precomputations are required here, this is why the result
	 * is a function that can (should) be reused.
	*/
	inline Rotation rotate(size_t slots) const {
		return block_rotate(slots, slot_group_len());
	}

	inline Frobenius frobenius(size_t iters) const {
		return Frobenius(iters, *this);
	}

	inline RawAuto raw_auto(size_t g1_power, size_t g2_power) const {
		return RawAuto(*this, std::make_tuple(g1_power, g2_power));
	}
	uint64_t g1_ord() const noexcept;
	uint64_t g2_ord() const noexcept;
	bool g1_generated_slotgroup() const noexcept;

	SubringView R() const noexcept;
	SubringView slot() const noexcept;

	poly extract_slot_value(poly x, size_t slot) const;
	poly from_slot_value(const poly& x, size_t slot) const;
	const poly& slot_one(size_t slot) const noexcept;
	uint64_t prime() const noexcept;
	size_t exponent() const noexcept;
	uint64_t N() const noexcept;
	const seal::Modulus& index_mod() const noexcept;

	//tex:
	//Checks that x_power_n is a polynomial in $(\mathbb{Z}/p^e\mathbb{Z})[X^d]$, i.e. has only
	//monomials of the form $X^{kd}$.
	//If $d \ | \ n$, this means that in the ring $(\mathbb{Z}/p^e\mathbb{Z})[X]/(X^n - \text{x_power_n})$,
	//only the powers $X^{kd}$ of X have nonzero trace.
	static bool is_poly_sparse(const poly& x_power_n, size_t d);

	//tex:
    //Returns whether the factors of $X^N + 1$ are in $(\mathbb{Z}/p^e\mathbb{Z})[X^{d/2}]$, i.e.
	//are polynomials in $X^{d/2}$. In this case, also the "unit vectors" $e_i$
	//are polynomials in $X^{d/2}$, which simplifies many computations.
    //Currently, only the case that this is true is supported.
	bool is_d_half_sparse() const;
	bool is_d_sparse() const;

	//tex:
    //Returns the subring generated by $X^k$ where $k = N / 2^{\text{index_log2}}$
	SlotRing power_x_subring(size_t index_log2) const;
	void in_power_x_subring(const SlotRing& subring, const poly& el, poly& result) const;
	void from_power_x_subring(const SlotRing& subring, const poly& el, poly& result) const;
	SlotRing change_exponent(uint64_t new_exp) const;

	friend void test_g_automorphisms();
};

std::shared_ptr<SlotRing> small_test_parameters();
std::shared_ptr<SlotRing> small_p_257_test_parameters();
std::shared_ptr<SlotRing> p_127_test_parameters();
std::shared_ptr<SlotRing> p_257_test_parameters();
void test_g_automorphisms();
void test_block_rotate();
void test_rotate_noncyclic();

inline size_t SlotRing::slot_group_dims() const noexcept
{
	return 1;
}

inline size_t SlotRing::slot_group_len() const noexcept
{
	return n;
}

inline size_t SlotRing::slot_rank() const noexcept
{
	return d;
}

inline SlotRing::SubringView SlotRing::R() const noexcept
{
	return { mod_pe, mod_2N_cyclotomic };
}

inline SlotRing::SubringView SlotRing::slot() const noexcept
{
	return { mod_pe, slot_modulus };
}

inline uint64_t SlotRing::g1_ord() const noexcept
{
	return N() / 2;
}

inline uint64_t SlotRing::g2_ord() const noexcept
{
	return 2;
}

inline bool SlotRing::g1_generated_slotgroup() const noexcept
{
	return slot_group_cyclic;
}

inline uint64_t SlotRing::prime() const noexcept
{
	return p;
}

inline size_t SlotRing::exponent() const noexcept
{
	return e;
}

inline const seal::Modulus& SlotRing::index_mod() const noexcept
{
	return mod_2N;
}

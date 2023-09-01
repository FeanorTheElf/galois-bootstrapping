#pragma once
#include "seal/seal.h"
#include "transform.h"
#include "contextchain.h"

inline void debug_decrypt_and_print(const seal::SEALContext& context, const seal::Ciphertext& ct, const seal::SecretKey& sk, const SlotRing& slot_ring) {
	std::cout << std::endl << "=============== Out ==============" << std::endl << std::endl;
	seal::Decryptor dec(context, sk);
	seal::Plaintext plain;
	dec.decrypt(ct, plain);
	poly result_poly(plain.data(), plain.data() + plain.coeff_count());
	result_poly.resize(slot_ring.N());
	for (size_t i = 0; i < 1; ++i) {
		print_poly(std::cout, slot_ring.extract_slot_value(result_poly, i)) << std::endl << std::endl;
	}
}

inline void debug_decrypt_and_print(const ContextChain& context_chain, const seal::Ciphertext& ct, const seal::SecretKey& sk, const SlotRing& slot_ring) {
	std::cout << std::endl << "=============== Out ==============" << std::endl << std::endl;
	const size_t context_index = context_chain.get_context_index(ct.parms_id());
	seal::SecretKey used_sk;
	context_chain.convert_sk(sk, used_sk);
	seal::Decryptor dec(context_chain.target_context(), used_sk);
	seal::Ciphertext switched_ciphertext = ct;
	context_chain.multiply_switch_inplace(switched_ciphertext, context_chain.size() - context_index - 1);
	seal::Plaintext plain;
	std::cout << "Noise budget: " << dec.invariant_noise_budget(switched_ciphertext) << std::endl;
	dec.decrypt(switched_ciphertext, plain);
	poly result_poly(plain.data(), plain.data() + plain.coeff_count());
	result_poly.resize(slot_ring.N());
	for (size_t i = 0; i < 2; ++i) {
		print_poly(std::cout, slot_ring.extract_slot_value(result_poly, i)) << std::endl << std::endl;
	}
}

inline void debug_log_noise_budget(const ContextChain& context_chain, const seal::Ciphertext& ct, const seal::SecretKey& sk, const SlotRing& slot_ring) {
	const size_t context_index = context_chain.get_context_index(ct.parms_id());
	seal::SecretKey used_sk;
	context_chain.convert_sk(sk, used_sk);
	seal::Decryptor dec(context_chain.target_context(), used_sk);
	seal::Ciphertext switched_ciphertext = ct;
	context_chain.multiply_switch_inplace(switched_ciphertext, context_chain.size() - context_index - 1);
	seal::Plaintext plain;
	std::cout << "Current noise budget: " << dec.invariant_noise_budget(switched_ciphertext) << std::endl;
}

inline void debug_decrypt_and_print_raw(const ContextChain& context_chain, const seal::Ciphertext& ct, const seal::SecretKey& sk, const SlotRing& slot_ring) {
	std::cout << std::endl << "=============== Out ==============" << std::endl << std::endl;
	const size_t context_index = context_chain.get_context_index(ct.parms_id());
	seal::SecretKey used_sk;
	context_chain.convert_sk(sk, used_sk);
	seal::Decryptor dec(context_chain.target_context(), used_sk);
	seal::Ciphertext switched_ciphertext = ct;
	context_chain.multiply_switch_inplace(switched_ciphertext, context_chain.size() - context_index - 1);
	std::cout << "Noise budget: " << dec.invariant_noise_budget(switched_ciphertext) << std::endl;
	seal::Plaintext plain;
	dec.decrypt(switched_ciphertext, plain);
	poly result_poly(plain.data(), plain.data() + plain.coeff_count());
	result_poly.resize(slot_ring.N());
	print_poly(std::cout, result_poly) << std::endl << std::endl;
}

//#define DEBUG_PARAMS const seal::SecretKey& debug_sk
//#define DEBUG_PASS(sk) sk
//#define DEBUG_DEC_PRINT(context, el) debug_decrypt_and_print(context, el, debug_sk, *slot_ring)
//#define DEBUG_DEC_PRINT_RAW(context, el) debug_decrypt_and_print_raw(context, el, debug_sk, *slot_ring)
//#define DEBUG_LOG_NOISE_BUDGET(context, el) debug_log_noise_budget(context, el, debug_sk, *slot_ring)

#define DEBUG_PARAMS
#define DEBUG_PASS(x)
#define DEBUG_DEC_PRINT(x, y)
#define DEBUG_DEC_PRINT_RAW(x, y)
#define DEBUG_LOG_NOISE_BUDGET(x, y)

/**
 * Encapsulates the computation of the norm in each slot
*/
class SlotwiseNorm {

	const SlotRing& slot_ring;
	size_t subfield_index_log2;

public:
	SlotwiseNorm(const SlotRing& slot_ring, size_t subfield_index_log2);
	SlotwiseNorm(const SlotwiseNorm&) = default;
	SlotwiseNorm(SlotwiseNorm&&) = default;
	~SlotwiseNorm() = default;

	poly operator()(const poly& x) const;
	void apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, const seal::RelinKeys& rk, seal::Ciphertext& destination) const;
	std::vector<uint32_t> galois_elements() const;
};

/**
 * Encapsulates the computation of the trace in each slot
*/
class SlotwiseTrace {

	const SlotRing& slot_ring;
	size_t source_subfield_index_log2;
	size_t target_subfield_index_log2;

public:
	SlotwiseTrace(const SlotRing& slot_ring, size_t source_subfield_index_log2, size_t target_subfield_index_log2);
	SlotwiseTrace(const SlotwiseTrace&) = default;
	SlotwiseTrace(SlotwiseTrace&&) = default;
	~SlotwiseTrace() = default;

	poly operator()(const poly& x) const;
	void apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& destination) const;
	std::vector<uint32_t> galois_elements() const;
	size_t field_index() const;
};

/**
 * Abstract base class for all objects that can evaluate the digit extraction
 * polynomial. This includes mainly the "general" 3log(d) and the "special" 2log(d)
 * algorithms as described in the paper. 
*/
class PolyEvaluator {

public:
	PolyEvaluator() = default;
	PolyEvaluator(const PolyEvaluator&) = delete;
	PolyEvaluator(PolyEvaluator&&) = delete;
	virtual ~PolyEvaluator() = default;

	virtual poly operator()(const poly&) const = 0;
	virtual void apply_ciphertext(const seal::Ciphertext& in, const seal::SEALContext& context, const seal::Evaluator& eval, const seal::GaloisKeys& gk, const seal::RelinKeys& rk, seal::Ciphertext& dst) const = 0;
	virtual std::vector<uint32_t> galois_elements() const = 0;
};

//tex:
//Evaluates the digit extraction polynomial in the $p = 127$ case;
//The formula is $$N(\text{evaluation_element} - x) + \text{correction_poly}(x) + x^{2^\text{log2_exponent}}$$
class P127PolyEvaluator : public PolyEvaluator {

	const SlotRing& slot_ring;
	poly evaluation_element;
	size_t log2_exponent;
	SlotwiseNorm norm_op;
	// we do not include X^N in the correction poly to keep it short
	poly correction_poly;

public:
	P127PolyEvaluator(const SlotRing& slot_ring, poly evaluation_element, poly correction_poly, size_t log2_exponent);
	virtual ~P127PolyEvaluator() = default;

	virtual poly operator()(const poly& x) const;
	virtual void apply_ciphertext(const seal::Ciphertext& in, const seal::SEALContext& context, const seal::Evaluator& eval, const seal::GaloisKeys& gk, const seal::RelinKeys& rk, seal::Ciphertext& destination) const;
	virtual std::vector<uint32_t> galois_elements() const;
};

//tex:
//Evaluates the digit extraction polynomial in the case $p = 257$;
//The formula is $$x \cdot N(\text{evaluation_element} - x) + \text{correction_poly}(x)$$
//assuming that the poly to evaluate is monic and of degree $2^{\text{log2_exponent}} + 1$
class P257CorrectionPolyEvaluator : public PolyEvaluator {

	const SlotRing& slot_ring;
	poly evaluation_element;
	size_t log2_exponent;
	SlotwiseNorm norm_op;
	poly correction_poly;
	uint64_t constant_correction;

public:
	P257CorrectionPolyEvaluator(const SlotRing& slot_ring, poly evaluation_element, poly correction_poly, size_t log2_exponent);
	virtual ~P257CorrectionPolyEvaluator() = default;

	virtual poly operator()(const poly& x) const;
	virtual void apply_ciphertext(const seal::Ciphertext& in, const seal::SEALContext& context, const seal::Evaluator& eval, const seal::GaloisKeys& gk, const seal::RelinKeys& rk, seal::Ciphertext& destination) const;
	virtual std::vector<uint32_t> galois_elements() const;
};

class Bootstrapper;

class BootstrappingKey {

	seal::Ciphertext encrypted_sk;
	seal::GaloisKeys gk;
	seal::RelinKeys rk;

	bool is_valid_for(const seal::SEALContext& context) const;

public:
	BootstrappingKey() = default;
	BootstrappingKey(const BootstrappingKey&) = default;
	BootstrappingKey(BootstrappingKey&&) = default;
	~BootstrappingKey() = default;

	const seal::GaloisKeys& galois_keys() const;
	const seal::RelinKeys& relin_keys() const;

	friend Bootstrapper;
};

/**
 * Hardcoded digit extraction polynomial in the p = 127 case
*/
std::unique_ptr<PolyEvaluator> p_127_test_parameters_digit_extractor(const SlotRing& slot_ring);

/**
 * Hardcoded digit extraction polynomial in the p = 257 case
*/
std::unique_ptr<PolyEvaluator> p_257_test_parameters_digit_extractor(const SlotRing& slot_ring);

void test_galois_poly_evaluator();
void test_homomorphic_noisy_decrypt();
void test_slotwise_digit_extract();
void test_coeffs_to_slots();

class Bootstrapper {

public:

	ContextChain context_chain;
	seal::MemoryPoolHandle pool;
	std::shared_ptr<const SlotRing> slot_ring;

	std::unique_ptr<PolyEvaluator> digit_extract_poly;
	std::unique_ptr<SlotwiseTrace> trace_op;
	std::unique_ptr<CompiledSubringLinearTransform> slots_to_coefficients = nullptr;
	std::unique_ptr<CompiledSubringLinearTransform> coefficients_to_slots = nullptr;

	size_t poly_modulus_degree() const noexcept;
	const seal::Evaluator& bootstrapping_evaluator() const;
	const seal::SEALContext& bootstrapping_context() const;

public:

	/**
	 * Computes the "noisy decryption" `c0 + c1 * s`
	*/
	void homomorphic_noisy_decrypt(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const;
	
	/**
	 * Computes the digit extraction step in each slot.
	 * 
	 * Details: The input must be a ciphertext w.r.t. plaintext modulus `p^e` that contains
	 * a scalar value `x[i]` in each slot. The output is then a ciphertext w.r.t. plaintext modulus `p` that
	 * contains the scalar values `round(x[i] / p^(e - 1))` in each slot.
	*/
	void slotwise_digit_extract(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const;
	
	/**
	 * Computes a ciphertext encrypting a plaintext whose coefficients w.r.t. X^d are the values that were
	 * stored in the slots before. This requires that every slot only contains a scalar value.
	*/
	void slots_to_coeffs(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const;
	
	/**
	 * Computes a ciphertext with a scalar value in each slot.
	 * The value stored in the i-th slot is the coefficient of X^di in the encrypted input plaintext.
	 * The coefficients of powers of X not of the form di are ignored.
	*/
	void coeffs_to_slots(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const;

public:

	Bootstrapper(const seal::SEALContext& bootstrapped_context, std::shared_ptr<const SlotRing> slot_ring, std::unique_ptr<PolyEvaluator> digit_extract_poly, MemoryPoolHandle pool = MemoryManager::GetPool());
	Bootstrapper(const Bootstrapper&) = default;
	Bootstrapper(Bootstrapper&&) = default;
	~Bootstrapper() = default;

	void initialize();
	void create_bootstrapping_key(const seal::SecretKey& sk, BootstrappingKey& bk, bool allow_uninitialized_for_test = false) const;
	void create_secret_key(const seal::SecretKey& base_sk, seal::SecretKey& destination) const;

	friend void test_homomorphic_noisy_decrypt();
	friend void test_slotwise_digit_extract();
};

inline size_t Bootstrapper::poly_modulus_degree() const noexcept
{
	return context_chain.poly_modulus_degree();
}

inline const seal::Evaluator& Bootstrapper::bootstrapping_evaluator() const
{
	return context_chain.target_evaluator();
}

inline const seal::GaloisKeys& BootstrappingKey::galois_keys() const
{
	return gk;
}

inline const seal::RelinKeys& BootstrappingKey::relin_keys() const
{
	return rk;
}
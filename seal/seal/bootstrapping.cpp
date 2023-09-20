#include "bootstrapping.h"
#include "seal/util/scalingvariant.h"
#include "seal/util/polyarithsmallmod.h"

using namespace seal;

// See the issue https://github.com/microsoft/SEAL/issues/592
void fast_exponentiate(const Ciphertext& ciphertext, size_t exponent, const Evaluator& evaluator, const RelinKeys& rk, Ciphertext& result) {
	// we let result uninitialized until we have to, this is represented by result_is_one == true
	bool result_is_one = true;
	
	// a classical square-and-multiply algorithm
	for (int64_t i = std::numeric_limits<size_t>::digits - 1; i >= 0; --i) {
		if (!result_is_one) {
			evaluator.square_inplace(result);
			evaluator.relinearize_inplace(result, rk); log_relin();
		}
		if (((exponent >> i) & 1) == 1) {
			// multiplying ciphertext to the result here will cause it to be
			// squared exactly i times
			if (!result_is_one) {
				evaluator.multiply_inplace(result, ciphertext);
				evaluator.relinearize_inplace(result, rk); log_relin();
			}
			else {
				result = ciphertext;
				result_is_one = false;
			}
		}
	}
}

void Bootstrapper::homomorphic_noisy_decrypt(const Ciphertext& ciphertext, const BootstrappingKey& bk, Ciphertext& destination, MemoryPoolHandle pool DEBUG_PARAMS) const
{
	// ciphertext is validated in convert_ciphertext_plain()
	if (!bk.is_valid_for(context_chain.target_context())) {
		throw std::invalid_argument("Invalid bootstrapping key");
	}
	seal::Plaintext ciphertext_as_plaintext[2];
	context_chain.convert_ciphertext_plain(ciphertext, ciphertext_as_plaintext, pool);
	bootstrapping_evaluator().multiply_plain(bk.encrypted_sk, ciphertext_as_plaintext[1], destination, pool);
	bootstrapping_evaluator().add_plain_inplace(destination, ciphertext_as_plaintext[0]);
}

Bootstrapper::Bootstrapper(const seal::SEALContext& bootstrapped_context, std::shared_ptr<const SlotRing> slot_ring, std::unique_ptr<PolyEvaluator> digit_extract_poly, MemoryPoolHandle pool)
	:	context_chain(bootstrapped_context, slot_ring->exponent(), pool),
		slot_ring(slot_ring),
		digit_extract_poly(std::move(digit_extract_poly)),
		pool(pool)
{
	// context compatibility is checked by ContextChain

	if (slot_ring->is_d_sparse()) {
		trace_op = std::make_unique<SlotwiseTrace>(*slot_ring, 0, log2_exact(slot_ring->slot_rank()));
	}
	else {
		trace_op = std::make_unique<SlotwiseTrace>(*slot_ring, 0, log2_exact(slot_ring->slot_rank() / 2));
	}
}

void Bootstrapper::create_bootstrapping_key(const seal::SecretKey& sk, BootstrappingKey& bk, bool allow_uninitialized_for_test) const
{
	assert(coefficients_to_slots != nullptr || allow_uninitialized_for_test);
	// sk validation is done by convert_secret_key_plain()

	Plaintext sk_plain(MemoryManager::GetPool(mm_prof_opt::mm_force_new, true));
	context_chain.convert_sk_plain(sk, sk_plain);

	SecretKey bootstrapping_sk;
	context_chain.convert_sk(sk, bootstrapping_sk);

	KeyGenerator keygen(context_chain.target_context(), bootstrapping_sk);
	PublicKey pk;
	keygen.create_public_key(pk);

	Encryptor enc(context_chain.target_context(), pk);
	seal::Ciphertext enc_sk;
	enc.encrypt(sk_plain, enc_sk, pool);

	std::vector<uint32_t> galois_elements;
	auto digit_extract_galois_elements = digit_extract_poly->galois_elements();
	galois_elements.insert(galois_elements.end(), digit_extract_galois_elements.cbegin(), digit_extract_galois_elements.cend());

	if (coefficients_to_slots != nullptr) {
		auto lin_transform_galois_elements = coefficients_to_slots->galois_elements();
		galois_elements.insert(galois_elements.end(), lin_transform_galois_elements.cbegin(), lin_transform_galois_elements.cend());
		lin_transform_galois_elements = slots_to_coefficients->galois_elements();
		galois_elements.insert(galois_elements.end(), lin_transform_galois_elements.cbegin(), lin_transform_galois_elements.cend());
	}

	auto trace_op_galois_elements = trace_op->galois_elements();
	galois_elements.insert(galois_elements.end(), trace_op_galois_elements.cbegin(), trace_op_galois_elements.cend());
	std::sort(galois_elements.begin(), galois_elements.end());
	const auto end = std::unique(galois_elements.begin(), galois_elements.end());
	galois_elements.resize(end - galois_elements.begin());
	std::cout << "creating " << galois_elements.size() << " galois keys" << std::endl;

	GaloisKeys gk;
	keygen.create_galois_keys(galois_elements, gk);

	RelinKeys rk;
	keygen.create_relin_keys(rk);

	bk.encrypted_sk = std::move(enc_sk);
	bk.gk = std::move(gk);
	bk.rk = std::move(rk);
}

void Bootstrapper::create_secret_key(const seal::SecretKey& base_sk, seal::SecretKey& destination) const
{
	context_chain.convert_sk(base_sk, destination);
}

std::unique_ptr<PolyEvaluator> p_127_test_parameters_digit_extractor(const SlotRing& slot_ring)
{
	assert(slot_ring.prime() == 127);
	assert(slot_ring.exponent() == 3);
	assert(slot_ring.N() == (size_t)1 << 15);
	// these values are computed with SAGE and hardcoded;
	// Computing them requires factoring a polynomial over a finite field and
	// then Hensel-lifting the result up to the slot ring
	auto coefficients = { 24040, 1645466, 638722, 310629, 373759, 2040343, 1702357, 1012959, 1330654, 972879, 85062, 95640, 1896466, 262437, 471266, 1157648, 1336311, 1250040, 622841, 1418600, 854122, 448430, 1113661, 1076329, 67242, 1020027, 1315226, 1004280, 1890984, 510873, 1490450, 1664995, 1524537, 719791, 24100, 1731453, 1055624, 163912, 297023, 738333, 877689, 1121258, 744242, 323617, 1961469, 26253, 1441237, 414438, 1012890, 1345724, 1094738, 1430813, 1926035, 913882, 1431639, 65384, 702073, 1479718, 930179, 827350, 213339, 895486, 162364, 1534952, 1261135, 216259, 908225, 1390823, 811185, 773571, 858483, 116902, 1485134, 95760, 741014, 1719246, 1988449, 843628, 572664, 1474961, 1689397, 140558, 1288343, 1763928, 90415, 1540960, 694622, 501217, 1840446, 1280873, 564935, 833880, 1975345, 1308774, 460575, 1378421, 307978, 1778634, 791153, 860979, 1929116, 1816565, 387144, 780466, 575449, 1541672, 1408027, 1471407, 483539, 1082463, 1588493, 1881429, 154889, 1767856, 1831001, 555167, 109362, 1552978, 1099674, 1224673, 1596695, 1518968, 1194351, 155081, 1016728, 96717, 1447799, 1026354 };
	auto indices = { 508, 504, 500, 496, 492, 488, 484, 480, 476, 472, 468, 464, 460, 456, 452, 448, 444, 440, 436, 432, 428, 424, 420, 416, 412, 408, 404, 400, 396, 392, 388, 384, 380, 376, 372, 368, 364, 360, 356, 352, 348, 344, 340, 336, 332, 328, 324, 320, 316, 312, 308, 304, 300, 296, 292, 288, 284, 280, 276, 272, 268, 264, 260, 256, 252, 248, 244, 240, 236, 232, 228, 224, 220, 216, 212, 208, 204, 200, 196, 192, 188, 184, 180, 176, 172, 168, 164, 160, 156, 152, 148, 144, 140, 136, 132, 128, 124, 120, 116, 112, 108, 104, 100, 96, 92, 88, 84, 80, 76, 72, 68, 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0 };

	poly element;
	element.resize(512);
	auto index_it = indices.begin();
	auto coeff_it = coefficients.begin();
	for (; index_it != indices.end(); ++index_it, ++coeff_it) {
		element[*index_it] = *coeff_it;
	}
	
	poly correction_poly = { util::negate_uint_mod(63, slot_ring.R().scalar_mod) };
	size_t log2_exponent = 7;
	return std::make_unique<P127PolyEvaluator>(slot_ring, std::move(element), std::move(correction_poly), log2_exponent);
}

std::unique_ptr<PolyEvaluator> p_257_test_parameters_digit_extractor(const SlotRing& slot_ring)
{
	assert(slot_ring.prime() == 257);
	assert(slot_ring.exponent() <= 3);
	assert(slot_ring.N() == (size_t)1 << 15);
	// these values are computed with SAGE and hardcoded;
	// Computing them requires factoring a polynomial over a finite field and
	// then Hensel-lifting the result up to the slot ring
	auto coefficients = { 7115559, 16940412, 10241707, 9687872, 7201911, 16887470, 9372790, 1221264, 6275169, 14505080, 8819983, 14996464, 2406034, 5962143, 7813571, 4005088, 4998650, 3602369, 11674739, 4680484, 11873400, 6752675, 1240539, 3703884, 1647627, 7878849, 9589184, 539957, 15866923, 1039565, 12618957, 5853175, 864291, 14372211, 2076817, 5417046, 10501534, 2150062, 16916511, 378561, 4452782, 14364758, 9375360, 16639208, 8940516, 14678041, 4137957, 12331374, 8516209, 8083421, 9679391, 15762067, 16385806, 1336400, 12016292, 15577541, 15859213, 16426926, 5464333, 12789605, 8899653, 5446858, 10508216, 16479097, 5746263, 7952608, 4592847, 8471748, 6331452, 7836187, 6641651, 5342002, 3759910, 15072279, 9800438, 9391808, 12205187, 14590147, 15296897, 3101990, 6596933, 16169412, 8059777, 16644091, 10922757, 9584558, 6265917, 9342207, 13596842, 13144779, 14964596, 1520155, 11831766, 5055961, 15309233, 9321390, 169877, 15592961, 6079335, 11224218, 3642718, 3033371, 4996080, 1847059, 13589646, 2120250, 14092595, 1328433, 16499914, 2382647, 3506251, 15659781, 8704590, 7576874, 16481924, 13506121, 13873888, 2528623, 2780483, 15841223, 4194240, 985338, 10900398, 2961925, 1105100, 6917155, 4652471, 3964482, 8487425 };
	auto indices = { 255, 253, 251, 249, 247, 245, 243, 241, 239, 237, 235, 233, 231, 229, 227, 225, 223, 221, 219, 217, 215, 213, 211, 209, 207, 205, 203, 201, 199, 197, 195, 193, 191, 189, 187, 185, 183, 181, 179, 177, 175, 173, 171, 169, 167, 165, 163, 161, 159, 157, 155, 153, 151, 149, 147, 145, 143, 141, 139, 137, 135, 133, 131, 129, 127, 125, 123, 121, 119, 117, 115, 113, 111, 109, 107, 105, 103, 101, 99, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1, 0 };

	poly element;
	element.resize(256);
	auto index_it = indices.begin();
	auto coeff_it = coefficients.begin();
	for (; index_it != indices.end(); ++index_it, ++coeff_it) {
		element[*index_it] = slot_ring.R().scalar_mod.reduce(*coeff_it);
	}

	poly correction_poly = { 0, util::negate_uint_mod(3, slot_ring.R().scalar_mod) };
	size_t log2_exponent = 8;
	return std::make_unique<P257CorrectionPolyEvaluator>(slot_ring, std::move(element), std::move(correction_poly), log2_exponent);
}

P127PolyEvaluator::P127PolyEvaluator(const SlotRing& slot_ring, poly evaluation_element, poly correction_poly, size_t log2_exponent)
	: slot_ring(slot_ring), correction_poly(std::move(correction_poly)), log2_exponent(log2_exponent), norm_op(slot_ring, log2_exact(slot_ring.slot_rank()) - log2_exponent)
{
	if (evaluation_element.size() != slot_ring.slot_rank()) {
		throw std::invalid_argument("Elements must have slot_rank() coefficients.");
	}
	poly slotted_evaluation_element;
	SEAL_ITERATE(util::iter(size_t(0)), slot_ring.slot_group_len(), [&slotted_evaluation_element, &evaluation_element, &slot_ring](auto I) {
		poly_add(slotted_evaluation_element, slot_ring.from_slot_value(evaluation_element, I), slot_ring.R().scalar_mod);
	});
	this->evaluation_element = std::move(slotted_evaluation_element);
}

poly P127PolyEvaluator::operator()(const poly& x) const
{
	throw std::invalid_argument("Unimplemented");
}

/**
 * Homomorphically evaluates the polynomial at x. This uses Horner's schema
 * with very bad multiplicative depth, so use only for very smal polynomials!
 * (we only use it for degree-1-polynomial I think)
*/
void naive_add_poly_eval_inplace(Ciphertext& destination, const poly& poly, const Ciphertext& x, const Evaluator& eval, const RelinKeys& rk, const Modulus& plain_modulus) {
	if (poly.size() == 0) {
		return;
	}
	assert(poly[poly.size() - 1] != 0);

	Ciphertext correction;
	bool is_correction_set = false;
	Plaintext current;
	current.resize(1);
	// this loop will not go over the constant coefficient - we include that later
	SEAL_ITERATE(poly.rbegin(), poly.size() - 1, [&](auto& I) {
		current[0] = plain_modulus.reduce(I);
		if (is_correction_set) {
			eval.add_plain_inplace(correction, current);
			eval.multiply_inplace(correction, x);
			eval.relinearize_inplace(correction, rk); log_relin();
		}
		else if (current[0] != 0) {
			eval.multiply_plain(x, current, correction);
			is_correction_set = true;
		}
	});
	if (is_correction_set) {
		eval.add_inplace(destination, correction);
	}
	current[0] = poly[0];
	eval.add_plain_inplace(destination, current);
}

void P127PolyEvaluator::apply_ciphertext(const Ciphertext& in, const SEALContext& context, const Evaluator& eval, const GaloisKeys& gk, const RelinKeys& rk, Ciphertext& destination) const
{
	if (!is_metadata_valid_for(in, context)) {
		throw std::invalid_argument("Invalid ciphertext");
	}

	const Modulus& plain_modulus = context.first_context_data()->parms().plain_modulus();
	Plaintext evaluation_element(gsl::span(this->evaluation_element));
	util::modulo_poly_coeffs(
		util::ConstCoeffIter(this->evaluation_element.data()), 
		this->evaluation_element.size(), 
		plain_modulus,
		util::CoeffIter(evaluation_element.data())
	);
	
	// first, compute the norm of (evaluation_element - in); this is the result up to a correction delta
	Ciphertext tmp;
	eval.sub_plain(in, evaluation_element, tmp);
	eval.negate_inplace(tmp);
	norm_op.apply_ciphertext(tmp, eval, gk, rk, destination);

	// then compute the correction
	naive_add_poly_eval_inplace(destination, correction_poly, in, eval, rk, plain_modulus);

	fast_exponentiate(in, (size_t)1 << log2_exponent, eval, rk, tmp);
	eval.sub_inplace(destination, tmp);
}

std::vector<uint32_t> P127PolyEvaluator::galois_elements() const
{
	return norm_op.galois_elements();
}

SlotwiseNorm::SlotwiseNorm(const SlotRing& slot_ring, size_t subfield_index_log2) : slot_ring(slot_ring), subfield_index_log2(subfield_index_log2)
{
}

poly SlotwiseNorm::operator()(const poly& x) const
{
	poly current = x;
	poly copy;
	for (size_t i = 0; i < log2_exact(slot_ring.slot_rank()) - subfield_index_log2; ++i) {
		copy = slot_ring.frobenius((size_t)1 << i)(current);
		current = poly_mul_mod(
			current,
			copy,
			slot_ring.R().scalar_mod,
			slot_ring.R().poly_mod
		);
	}
	return current;
}

void SlotwiseNorm::apply_ciphertext(const Ciphertext& in, const Evaluator& eval, const GaloisKeys& gk, const RelinKeys& rk, Ciphertext& destination) const
{
	destination = in;
	seal::Ciphertext copy;
	for (size_t i = 0; i < log2_exact(slot_ring.slot_rank()) - subfield_index_log2; ++i) {
		eval.apply_galois(destination, static_cast<uint32_t>(seal::util::exponentiate_uint_mod(slot_ring.prime(), (size_t)1 << i, slot_ring.index_mod())), gk, copy); log_galois();
		eval.multiply_inplace(destination, copy);
		eval.relinearize_inplace(destination, rk); log_relin();
	}
}

std::vector<uint32_t> SlotwiseNorm::galois_elements() const
{
	std::vector<uint32_t> result;
	for (size_t i = 0; i < log2_exact(slot_ring.slot_rank()) - subfield_index_log2; ++i) {
		result.push_back(static_cast<uint32_t>(seal::util::exponentiate_uint_mod(slot_ring.prime(), (size_t)1 << i, slot_ring.index_mod())));
	}
	return result;
}

bool BootstrappingKey::is_valid_for(const seal::SEALContext& context) const
{
	return is_metadata_valid_for(gk, context) && is_metadata_valid_for(rk, context) && is_metadata_valid_for(encrypted_sk, context) && is_buffer_valid(gk) && is_buffer_valid(rk) && is_buffer_valid(encrypted_sk);
}

void Bootstrapper::slotwise_digit_extract(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const
{
	const size_t digits_to_remove = slot_ring->exponent() - 1;
	const size_t highest_digit_index = slot_ring->exponent() - 1;
	util::Pointer<Ciphertext> worktable = util::allocate<Ciphertext>(digits_to_remove, pool);

	// add modulus / 2 to perform rounding instead of flooring
	Plaintext shift;
	shift.resize(1);
	for (size_t i = 0; i < digits_to_remove; ++i) {
		shift.data()[0] *= slot_ring->prime();
		shift.data()[0] += slot_ring->prime() / 2;
	}
	destination = ciphertext;
	bootstrapping_evaluator().add_plain_inplace(destination, shift);

	SEAL_ITERATE(iter(worktable), digits_to_remove, [&destination](auto& I) {
		I = destination;
	});

	SEAL_ITERATE(iter(worktable, size_t(0)), digits_to_remove, [&](auto I) {
		const size_t current_index = std::get<1>(I);
		Ciphertext& current = std::get<0>(I);
		const size_t context_index = context_chain.size() - current_index - 1;
		const Evaluator& evaluator = context_chain.get_evaluator(context_index);
		const SEALContext& context = context_chain.get_context(context_index);
		RelinKeys rk;
		GaloisKeys gk;
		context_chain.convert_kswitchkey(bk.galois_keys(), gk, context_index);
		context_chain.convert_kswitchkey(bk.relin_keys(), rk, context_index);
		Ciphertext tmp;
		context_chain.divide_exact_switch_inplace(current, current_index);

		SEAL_ITERATE(util::iter(size_t(1)), highest_digit_index - current_index, [&](auto J) {
			digit_extract_poly->apply_ciphertext(current, context, evaluator, gk, rk, tmp);
			current = tmp;
			if (current_index + J < digits_to_remove) {
				context_chain.multiply_switch_inplace(tmp, current_index);
				bootstrapping_evaluator().sub_inplace(worktable[current_index + J], tmp);
			}
		});
		DEBUG_LOG_NOISE_BUDGET(context_chain, current);
		context_chain.multiply_switch_inplace(current, current_index);
		bootstrapping_evaluator().sub_inplace(destination, current);
	});

	context_chain.divide_exact_switch_inplace(destination, digits_to_remove);
}

void Bootstrapper::slots_to_coeffs(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const
{
	GaloisKeys base_gk;
	context_chain.convert_kswitchkey(bk.galois_keys(), base_gk, 0);
	slots_to_coefficients->apply_ciphertext(ciphertext, context_chain.get_evaluator(0), base_gk, destination);
}

void Bootstrapper::coeffs_to_slots(const seal::Ciphertext& ciphertext, const BootstrappingKey& bk, seal::Ciphertext& destination, seal::MemoryPoolHandle pool DEBUG_PARAMS) const
{
	seal::Ciphertext tmp;
	// first we remove all the coefficients "to discard" by applying the trace
	trace_op->apply_ciphertext(ciphertext, bootstrapping_evaluator(), bk.galois_keys(), tmp);
	// fix the constant factor introduced by trace
	const uint64_t error_factor = slot_ring->R().scalar_mod.reduce(trace_op->field_index());
	const uint64_t correction_factor = inv_mod(error_factor, slot_ring->R().scalar_mod);
	Plaintext factor;
	factor.resize(1);
	factor[0] = correction_factor;
	bootstrapping_evaluator().multiply_plain_inplace(tmp, factor, pool);
	coefficients_to_slots->apply_ciphertext(tmp, bootstrapping_evaluator(), bk.galois_keys(), destination);
}

void Bootstrapper::initialize()
{
	if (coefficients_to_slots == nullptr) {
		coefficients_to_slots = std::make_unique<CompiledSubringLinearTransform>(CompiledSubringLinearTransform::coeffs_to_slots(slot_ring));
		slots_to_coefficients = std::make_unique<CompiledSubringLinearTransform>(CompiledSubringLinearTransform::slots_to_coeffs(slot_ring));
	}
}

const seal::SEALContext& Bootstrapper::bootstrapping_context() const
{
	return context_chain.target_context();
}

void test_galois_poly_evaluator()
{
	{
		EncryptionParameters parms(scheme_type::bfv);
		// we must use p_127_test_parameters() here as the slot rank for small_test_parameters() is too small for fast poly evaluation
		SlotRing slot_ring = *p_127_test_parameters();
		parms.set_poly_modulus_degree(slot_ring.N());
		parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring.N()));
		parms.set_plain_modulus(127 * 127 * 127);
		std::unique_ptr<PolyEvaluator> digit_extractor = p_127_test_parameters_digit_extractor(slot_ring);
		SEALContext context(parms);

		KeyGenerator keygen(context);
		SecretKey sk = keygen.secret_key();
		PublicKey pk;
		keygen.create_public_key(pk);
		RelinKeys rk;
		keygen.create_relin_keys(rk);
		GaloisKeys gk;
		keygen.create_galois_keys(digit_extractor->galois_elements(), gk);

		poly data;
		poly_add(data, slot_ring.from_slot_value({ 5 * 127 * 127 + 2 * 127 - 7 }, 0), slot_ring.R().scalar_mod);
		poly_add(data, slot_ring.from_slot_value({ 127 * 127 + 66 }, 1), slot_ring.R().scalar_mod);
		Plaintext x_plain{ gsl::span<const uint64_t>(data) };

		Encryptor encryptor(context, pk);
		Ciphertext x_enc;
		encryptor.encrypt(x_plain, x_enc);

		Ciphertext result_enc;
		Evaluator eval(context);
		digit_extractor->apply_ciphertext(x_enc, context, eval, gk, rk, result_enc);

		Decryptor decryptor(context, sk);
		Plaintext result;
		decryptor.decrypt(result_enc, result);

		poly result_poly(result.data(), result.data() + result.coeff_count());
		result_poly.resize(slot_ring.N());
		poly expected;
		poly_add(expected, slot_ring.from_slot_value({ 1822697 /* computed with sage */ }, 0), slot_ring.R().scalar_mod);
		poly_add(expected, slot_ring.from_slot_value({ 66 /* computed with sage */ }, 1), slot_ring.R().scalar_mod);

		assert(result_poly == expected);
	}
	{
		EncryptionParameters parms(scheme_type::bfv);
		// we must use p_127_test_parameters() here as the slot rank for small_test_parameters() is too small for fast poly evaluation
		SlotRing slot_ring = *p_257_test_parameters();
		parms.set_poly_modulus_degree(slot_ring.N());
		parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring.N()));
		parms.set_plain_modulus(slot_ring.R().scalar_mod.value());
		std::unique_ptr<PolyEvaluator> digit_extractor = p_257_test_parameters_digit_extractor(slot_ring);
		SEALContext context(parms);

		KeyGenerator keygen(context);
		SecretKey sk = keygen.secret_key();
		PublicKey pk;
		keygen.create_public_key(pk);
		RelinKeys rk;
		keygen.create_relin_keys(rk);
		GaloisKeys gk;
		keygen.create_galois_keys(digit_extractor->galois_elements(), gk);

		poly data;
		poly_add(data, slot_ring.from_slot_value({ 5 * 257 * 257 + 6 * 257 - 3 }, 0), slot_ring.R().scalar_mod);
		Plaintext x_plain{ gsl::span<const uint64_t>(data) };

		Encryptor encryptor(context, pk);
		Ciphertext x_enc;
		encryptor.encrypt(x_plain, x_enc);

		Ciphertext result_enc;
		Evaluator eval(context);
		digit_extractor->apply_ciphertext(x_enc, context, eval, gk, rk, result_enc);

		Decryptor decryptor(context, sk);
		Plaintext result;
		decryptor.decrypt(result_enc, result);

		poly result_poly(result.data(), result.data() + result.coeff_count());
		result_poly.resize(slot_ring.N());
		poly expected;
		poly_add(expected, slot_ring.from_slot_value({ 15521769 /* computed with sage */ }, 0), slot_ring.R().scalar_mod);

		assert(result_poly == expected);
	}
	std::cout << "test_galois_poly_evaluator(): success" << std::endl;
}

void test_homomorphic_noisy_decrypt() {
	{
		EncryptionParameters parms(scheme_type::bfv);
		std::shared_ptr<const SlotRing> slot_ring = p_127_test_parameters();
		parms.set_poly_modulus_degree(slot_ring->N());
		parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring->N()));
		parms.set_plain_modulus(slot_ring->prime());
		SEALContext context(parms, false, sec_level_type::none);

		KeyGenerator keygen(context);
		SecretKey sk = keygen.secret_key();
		PublicKey pk;
		keygen.create_public_key(pk);

		Encryptor encryptor(context, pk);
		Evaluator evaluator(context);

		std::vector<uint64_t> data = { 4 };
		Plaintext x_plain{ gsl::span<const uint64_t>(data) };

		Ciphertext x_enc;
		encryptor.encrypt(x_plain, x_enc);

		poly evaluation_element_fake;
		evaluation_element_fake.resize(slot_ring->slot_rank());
		Bootstrapper bootstrapper(context, slot_ring, std::make_unique<P127PolyEvaluator, const SlotRing&, poly, poly, size_t>(*slot_ring, std::move(evaluation_element_fake), { 1 }, 7));
		// dont initialize, as we do not require the linear transforms -> save time

		BootstrappingKey bk;
		bootstrapper.create_bootstrapping_key(sk, bk, true);

		Ciphertext noisy_enc_x;
		bootstrapper.homomorphic_noisy_decrypt(x_enc, bk, noisy_enc_x, MemoryManager::GetPool() DEBUG_PASS(sk));

		SecretKey bootstrapping_sk;
		bootstrapper.create_secret_key(sk, bootstrapping_sk);
		Decryptor decryptor(bootstrapper.bootstrapping_context(), bootstrapping_sk);

		Plaintext noisy_x;
		decryptor.decrypt(noisy_enc_x, noisy_x);

		auto div_rounded = [&slot_ring](uint64_t x) {
			return slot_ring->R().scalar_mod.reduce(x + (slot_ring->R().scalar_mod.value() / slot_ring->prime() / 2)) / (slot_ring->R().scalar_mod.value() / slot_ring->prime());
		};
		assert(div_rounded(noisy_x[0]) == 4);
		for (size_t i = 1; i < noisy_x.coeff_count(); ++i) {
			assert(div_rounded(noisy_x[i]) == 0);
		}
	} 
	{
		EncryptionParameters parms(scheme_type::bfv);
		std::shared_ptr<const SlotRing> slot_ring = p_257_test_parameters();
		parms.set_poly_modulus_degree(slot_ring->N());
		parms.set_coeff_modulus(CoeffModulus::Create(slot_ring->N(), { 40, 40, 40 }));
		parms.set_plain_modulus(slot_ring->prime());
		SEALContext context(parms, false, sec_level_type::none);

		KeyGenerator keygen(context);
		SecretKey sk = keygen.secret_key();
		PublicKey pk;
		keygen.create_public_key(pk);

		Encryptor encryptor(context, pk);
		Evaluator evaluator(context);

		std::vector<uint64_t> data = { 12 };
		Plaintext x_plain{ gsl::span<const uint64_t>(data) };

		Ciphertext x_enc;
		encryptor.encrypt(x_plain, x_enc);

		poly evaluation_element_fake;
		evaluation_element_fake.resize(slot_ring->slot_rank());
		Bootstrapper bootstrapper(context, slot_ring, std::make_unique<P127PolyEvaluator, const SlotRing&, poly, poly, size_t>(*slot_ring, std::move(evaluation_element_fake), { 1 }, 8));
		// dont initialize, as we do not require the linear transforms -> save time

		BootstrappingKey bk;
		bootstrapper.create_bootstrapping_key(sk, bk, true);

		Ciphertext noisy_enc_x;
		bootstrapper.homomorphic_noisy_decrypt(x_enc, bk, noisy_enc_x, MemoryManager::GetPool() DEBUG_PASS(sk));

		SecretKey bootstrapping_sk;
		bootstrapper.create_secret_key(sk, bootstrapping_sk);
		Decryptor decryptor(bootstrapper.bootstrapping_context(), bootstrapping_sk);

		Plaintext noisy_x;
		decryptor.decrypt(noisy_enc_x, noisy_x);

		auto div_rounded = [&slot_ring](uint64_t x) {
			return slot_ring->R().scalar_mod.reduce(x + (slot_ring->R().scalar_mod.value() / slot_ring->prime() / 2)) / (slot_ring->R().scalar_mod.value() / slot_ring->prime());
		};
		assert(div_rounded(noisy_x[0]) == 12);
		for (size_t i = 1; i < noisy_x.coeff_count(); ++i) {
			assert(div_rounded(noisy_x[i]) == 0);
		}
	}
	std::cout << "test_homomorphic_noisy_decrypt(): success" << std::endl;
}

void test_slotwise_digit_extract() {

	EncryptionParameters parms(scheme_type::bfv);
	// we must use p_127_test_parameters() here as the slot rank for small_test_parameters() is too small for fast poly evaluation
	std::shared_ptr<const SlotRing> slot_ring = p_127_test_parameters();
	parms.set_poly_modulus_degree(slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring->N()));
	parms.set_plain_modulus(slot_ring->prime());
	std::unique_ptr<PolyEvaluator> digit_extractor = p_127_test_parameters_digit_extractor(*slot_ring);
	SEALContext context(parms, false);

	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();

	poly data;
	poly_add(data, slot_ring->from_slot_value({ 5 * 127 * 127 + 2 * 127 - 7 * 127 }, 0), slot_ring->R().scalar_mod);
	poly_add(data, slot_ring->from_slot_value({ 111 * 127 * 127 - 8 * 127 + 5 * 127 }, 63), slot_ring->R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(data) };

	Bootstrapper bootstrapper(context, slot_ring, std::move(digit_extractor));
	// dont initialize, as we do not require the linear transforms -> save time

	SEALContext bootstrapping_context = bootstrapper.bootstrapping_context();
	SecretKey bootstrapping_sk;
	bootstrapper.create_secret_key(sk, bootstrapping_sk);
	KeyGenerator bootstrapping_keygen(bootstrapping_context, bootstrapping_sk);
	PublicKey bootstrapping_pk;
	bootstrapping_keygen.create_public_key(bootstrapping_pk);
	Encryptor encryptor(bootstrapping_context, bootstrapping_pk);

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);

	BootstrappingKey bk;
	bootstrapper.create_bootstrapping_key(sk, bk, true);

	Ciphertext result_enc;
	bootstrapper.slotwise_digit_extract(x_enc, bk, result_enc, MemoryManager::GetPool() DEBUG_PASS(sk));

	Decryptor decryptor(context, sk);
	Plaintext result;
	decryptor.decrypt(result_enc, result);

	SlotRing result_slot_ring = slot_ring->change_exponent(1);

	poly result_poly(result.data(), result.data() + result.coeff_count());
	result_poly.resize(result_slot_ring.N());
	poly expected;
	poly_add(expected, result_slot_ring.from_slot_value({ 5 }, 0), result_slot_ring.R().scalar_mod);
	poly_add(expected, result_slot_ring.from_slot_value({ 111 }, 63), result_slot_ring.R().scalar_mod);

	assert(result_poly == expected);
	std::cout << "test_slotwise_digit_extract(): success" << std::endl;
}

void test_coeffs_to_slots()
{
	EncryptionParameters parms(scheme_type::bfv);
	std::shared_ptr<const SlotRing> slot_ring = p_257_test_parameters();
	parms.set_poly_modulus_degree(slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring->N()));
	parms.set_plain_modulus(slot_ring->prime());
	SEALContext context(parms, false, sec_level_type::none);

	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();

	std::vector<uint64_t> data;
	data.resize(512);
	data[0] = 4;
	data[256] = 12;
	data[128] = 3;
	data[257] = 7;
	Plaintext x_plain{ gsl::span<const uint64_t>(data) };


	poly evaluation_element_fake;
	evaluation_element_fake.resize(slot_ring->slot_rank());
	Bootstrapper bootstrapper(context, slot_ring, std::make_unique<P127PolyEvaluator, const SlotRing&, poly, poly, size_t>(*slot_ring, std::move(evaluation_element_fake), { 1 }, 8));
	bootstrapper.initialize();

	SecretKey bootstrapping_sk;
	bootstrapper.create_secret_key(sk, bootstrapping_sk);
	KeyGenerator bootstrapping_keygen(bootstrapper.bootstrapping_context(), bootstrapping_sk);
	PublicKey pk;
	bootstrapping_keygen.create_public_key(pk);

	Encryptor encryptor(bootstrapper.bootstrapping_context(), pk);
	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);

	BootstrappingKey bk;
	bootstrapper.create_bootstrapping_key(sk, bk, true);

	Ciphertext noisy_enc_x;
	bootstrapper.coeffs_to_slots(x_enc, bk, noisy_enc_x, MemoryManager::GetPool() DEBUG_PASS(sk));

	Decryptor decryptor(bootstrapper.bootstrapping_context(), bootstrapping_sk);

	Plaintext noisy_x;
	decryptor.decrypt(noisy_enc_x, noisy_x);
	
	poly result_poly(noisy_x.data(), noisy_x.data() + noisy_x.coeff_count());
	result_poly.resize(slot_ring->N());
	assert(slot_ring->extract_slot_value(result_poly, 0)[0] == 4);
	assert(slot_ring->extract_slot_value(result_poly, 1)[0] == 12);
	for (size_t i = 2; i < slot_ring->slot_group_len(); ++i) {
		assert(is_zero(slot_ring->extract_slot_value(result_poly, i)));
	}

	std::cout << "test_coeffs_to_slots(): success" << std::endl;
}


SlotwiseTrace::SlotwiseTrace(const SlotRing& slot_ring, size_t source_subfield_index_log2, size_t target_subfield_index_log2) : slot_ring(slot_ring), source_subfield_index_log2(source_subfield_index_log2), target_subfield_index_log2(target_subfield_index_log2)
{
	assert(target_subfield_index_log2 > source_subfield_index_log2);
}

poly SlotwiseTrace::operator()(const poly& x) const
{
	poly current = x;
	poly copy;
	for (size_t i = log2_exact(slot_ring.slot_rank()) - target_subfield_index_log2; i < log2_exact(slot_ring.slot_rank()) - source_subfield_index_log2; ++i) {
		copy = slot_ring.frobenius((size_t)1 << i)(current);
		poly_add(
			current,
			copy,
			slot_ring.R().scalar_mod
		);
	}
	return current;
}

void SlotwiseTrace::apply_ciphertext(const Ciphertext& in, const Evaluator& eval, const GaloisKeys& gk, Ciphertext& destination) const
{
	destination = in;
	seal::Ciphertext copy;
	for (size_t i = log2_exact(slot_ring.slot_rank()) - target_subfield_index_log2; i < log2_exact(slot_ring.slot_rank()) - source_subfield_index_log2; ++i) {
		eval.apply_galois(destination, static_cast<uint32_t>(seal::util::exponentiate_uint_mod(slot_ring.prime(), (size_t)1 << i, slot_ring.index_mod())), gk, copy); log_galois();
		eval.add_inplace(destination, copy);
	}
}

std::vector<uint32_t> SlotwiseTrace::galois_elements() const
{
	std::vector<uint32_t> result;
	for (size_t i = log2_exact(slot_ring.slot_rank()) - target_subfield_index_log2; i < log2_exact(slot_ring.slot_rank()) - source_subfield_index_log2; ++i) {
		result.push_back(static_cast<uint32_t>(seal::util::exponentiate_uint_mod(slot_ring.prime(), (size_t)1 << i, slot_ring.index_mod())));
	}
	return result;
}

size_t SlotwiseTrace::field_index() const
{
	return (size_t)1 << (target_subfield_index_log2 - source_subfield_index_log2);
}

P257CorrectionPolyEvaluator::P257CorrectionPolyEvaluator(const SlotRing& slot_ring, poly evaluation_element, poly correction_poly, size_t log2_exponent)
	: slot_ring(slot_ring), constant_correction(0), log2_exponent(log2_exponent), norm_op(slot_ring, log2_exact(slot_ring.slot_rank()) - log2_exponent)
{
	if (evaluation_element.size() != slot_ring.slot_rank()) {
		throw std::invalid_argument("Elements must have slot_rank() coefficients.");
	}
	poly slotted_evaluation_element;
	SEAL_ITERATE(util::iter(size_t(0)), slot_ring.slot_group_len(), [&slotted_evaluation_element, &evaluation_element, &slot_ring](auto I) {
		poly_add(slotted_evaluation_element, slot_ring.from_slot_value(evaluation_element, I), slot_ring.R().scalar_mod);
	});
	this->evaluation_element = std::move(slotted_evaluation_element);

	if (correction_poly.size() != 0) {
		constant_correction = correction_poly[0];
		for (size_t i = 1; i < correction_poly.size(); ++i) {
			this->correction_poly.push_back(correction_poly[i]);
		}
	}
}

poly P257CorrectionPolyEvaluator::operator()(const poly& x) const
{
	throw std::invalid_argument("Unimplemented");
}

void P257CorrectionPolyEvaluator::apply_ciphertext(const seal::Ciphertext& in, const seal::SEALContext& context, const seal::Evaluator& eval, const seal::GaloisKeys& gk, const seal::RelinKeys& rk, seal::Ciphertext& destination) const
{
	if (!is_metadata_valid_for(in, context)) {
		throw std::invalid_argument("Invalid ciphertext");
	}
	const Modulus& plain_modulus = context.first_context_data()->parms().plain_modulus();
	Plaintext evaluation_element(gsl::span(this->evaluation_element));
	util::modulo_poly_coeffs(
		util::ConstCoeffIter(this->evaluation_element.data()),
		this->evaluation_element.size(),
		plain_modulus,
		util::CoeffIter(evaluation_element.data())
	);

	// first, compute the norm of (evaluation_element - in); this is the result up to a correction delta
	Ciphertext tmp;
	eval.sub_plain(in, evaluation_element, tmp);
	eval.negate_inplace(tmp);
	norm_op.apply_ciphertext(tmp, eval, gk, rk, destination);

	// then compute non-constant part of the correction
	naive_add_poly_eval_inplace(destination, correction_poly, in, eval, rk, plain_modulus);

	// multiply by x
	eval.multiply_inplace(destination, in);
	eval.relinearize_inplace(destination, rk); log_relin();

	// add constant part of correction
	Plaintext correction_constant;
	correction_constant.resize(1);
	correction_constant[0] = plain_modulus.reduce(constant_correction);
	eval.add_plain_inplace(destination, correction_constant);
}

std::vector<uint32_t> P257CorrectionPolyEvaluator::galois_elements() const
{
	return norm_op.galois_elements();
}

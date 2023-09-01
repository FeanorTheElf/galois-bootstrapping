#include "transform.h"
#include <fstream>
#include <string>
#include <numeric>


CompiledLinearTransform::CompiledLinearTransform(std::shared_ptr<const SlotRing> slot_ring, std::vector<poly> coefficients)
	: slot_ring(slot_ring), coefficients(std::move(coefficients))
{
	assert(slot_ring->g2_ord() == 2);
	assert(this->coefficients.size() == g1_subgroup_order() || this->coefficients.size() == 2 * g1_subgroup_order());
}

SlotRing::RawAuto CompiledLinearTransform::automorphism(size_t index) const
{
	const size_t g1_index = index % g1_subgroup_order();
	const size_t g2_index = (index % coefficients.size()) / g1_subgroup_order();
	return slot_ring->raw_auto(g1_index * (slot_ring->g1_ord() / g1_subgroup_order()), g2_index);
}

SlotRing::RawAuto CompiledLinearTransform::difference_automorphism(size_t from, size_t to) const
{
	const size_t from_g1_index = (from % g1_subgroup_order());
	const size_t from_g2_index = (from / g1_subgroup_order()) % g2_subgroup_order();
	const size_t to_g1_index = (to % g1_subgroup_order());
	const size_t to_g2_index = (to / g1_subgroup_order()) % g2_subgroup_order();
	size_t diff_g1_index = to_g1_index + g1_subgroup_order() - from_g1_index;
	if (diff_g1_index >= g1_subgroup_order()) {
		diff_g1_index -= g1_subgroup_order();
	}
	size_t diff_g2_index = to_g2_index + g2_subgroup_order() - from_g2_index;
	if (diff_g2_index >= g2_subgroup_order()) {
		diff_g2_index -= g2_subgroup_order();
	}
	return automorphism(diff_g1_index + diff_g2_index * g1_subgroup_order());
}

SlotRing::RawAuto CompiledLinearTransform::reverse_automorphism(size_t index) const
{
	return difference_automorphism(index, 0);
}

uint64_t CompiledLinearTransform::g1_subgroup_order() const
{
	return slot_ring->g1_ord();
}

uint64_t CompiledLinearTransform::g2_subgroup_order() const
{
	if (coefficients.size() == g1_subgroup_order()) {
		return 1;
	}
	else {
		assert(coefficients.size() == 2 * g1_subgroup_order());
		return 2;
	}
}

void CompiledLinearTransform::add_scaled_transform(const poly& scaling, const SlotRing::Rotation& rotation, const SlotRing::Frobenius& frobenius)
{
	{
		const uint64_t g1_power = (std::get<0>(rotation.get_forward_g1_g2_decomp()) + std::get<0>(frobenius.get_g1_g2_decomp())) % slot_ring->g1_ord();
		const uint64_t g2_power = (std::get<1>(rotation.get_forward_g1_g2_decomp()) + std::get<1>(frobenius.get_g1_g2_decomp())) % slot_ring->g2_ord();
		assert((g1_power * g1_subgroup_order()) % slot_ring->g1_ord() == 0);
		const size_t index = g1_power * g1_subgroup_order() / slot_ring->g1_ord() + g2_power * g1_subgroup_order();
		poly mask = slot_ring->raw_auto(g1_power, g2_power)(rotation.get_forward_mask());
		poly add = poly_mul_mod(scaling, mask, slot_ring->R().scalar_mod, slot_ring->R().poly_mod);

		poly_add(
			coefficients[index],
			add,
			slot_ring->R().scalar_mod
		);
	}
	{
		const uint64_t g1_power = (std::get<0>(rotation.get_backward_g1_g2_decomp()) + std::get<0>(frobenius.get_g1_g2_decomp())) % slot_ring->g1_ord();
		const uint64_t g2_power = (std::get<1>(rotation.get_backward_g1_g2_decomp()) + std::get<1>(frobenius.get_g1_g2_decomp())) % slot_ring->g2_ord();
		assert((g1_power * g1_subgroup_order()) % slot_ring->g1_ord() == 0);
		const size_t index = g1_power * g1_subgroup_order() / slot_ring->g1_ord() + g2_power * g1_subgroup_order();
		poly mask = slot_ring->raw_auto(g1_power, g2_power)(rotation.get_backward_mask());
		poly add = poly_mul_mod(scaling, mask, slot_ring->R().scalar_mod, slot_ring->R().poly_mod);
		
		poly_add(
			coefficients[index],
			add,
			slot_ring->R().scalar_mod
		);
	}
}

size_t CompiledLinearTransform::babystep_automorphism_count() const
{
	const size_t automorphism_count_log2 = log2_exact(coefficients.size());
	return (size_t)1 << (automorphism_count_log2 - automorphism_count_log2 / 2);
}

size_t CompiledLinearTransform::giantstep_automorphism_count() const
{
	const size_t automorphism_count_log2 = log2_exact(coefficients.size());
	return (size_t)1 << (automorphism_count_log2 / 2);
}

void CompiledLinearTransform::fix_coefficient_shift()
{
	for (size_t i = 0; i < giantstep_automorphism_count(); ++i) {
		for (size_t j = 0; j < babystep_automorphism_count(); ++j) {
			const auto automorphism = reverse_automorphism(i * babystep_automorphism_count());
			this->coefficients[i * babystep_automorphism_count() + j] = automorphism(this->coefficients[i * babystep_automorphism_count() + j]);
		}
	}
}

std::vector<poly> CompiledLinearTransform::compile_frobenius(
	const std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& sparse_transform_matrix,
	const SlotRing::SubringView& ring,
	uint64_t p,
	size_t N,
	const NegacyclicPowerTable& powertable
) {
	size_t d = ring.poly_mod.n;
	poly z = { 1 };

	const poly& x_power_n = ring.poly_mod.x_power_n;
	const seal::Modulus& mod_t = ring.scalar_mod;
	assert(SlotRing::is_poly_sparse(x_power_n, d / 2));

	uint64_t a = 0;
	if (d / 2 < x_power_n.size()) {
		a = x_power_n[d / 2];
	}
	const uint64_t b = x_power_n[0];

	// trace of x^0 = 1
	const uint64_t u0 = seal::util::barrett_reduce_64(d, mod_t);
	// trace of x^(d/2)
	const uint64_t u1 = seal::util::multiply_uint_mod(a, d / 2, mod_t);
	// trace of x^d
	const uint64_t u2 = seal::util::add_uint_mod(
		seal::util::multiply_uint_mod(
			seal::util::multiply_uint_mod(a, a, mod_t),
			d / 2,
			mod_t
		),
		seal::util::multiply_uint_mod(d, b, mod_t),
		mod_t
	);
	// trace of x^(3d/2)
	const uint64_t u3 = seal::util::multiply_uint_mod(
		seal::util::add_uint_mod(
			seal::util::exponentiate_uint_mod(a, 3, mod_t),
			seal::util::multiply_uint_mod(3, seal::util::multiply_uint_mod(a, b, mod_t), mod_t),
			mod_t
		),
		d / 2,
		mod_t
	);
	const std::array<uint64_t, 4> us = { u0, u1, u2, u3 };
	// factor to scale with at the end
	uint64_t global_factor = seal::util::negate_uint_mod(inv_mod(seal::util::sub_uint_mod(
		seal::util::multiply_uint_mod(u1, u3, mod_t),
		seal::util::multiply_uint_mod(u2, u2, mod_t),
		mod_t
	), mod_t), mod_t);

	seal::Modulus mod_2N(2 * N);
	poly tmp;
	std::vector<poly> results;
	for (size_t m = 0; m < d; ++m) {
		poly result;
		for (const auto& entry : sparse_transform_matrix) {
			size_t k = std::get<0>(entry.first);
			size_t l = std::get<1>(entry.first);
			assert(k <= d);
			assert(l <= d);
			int64_t delta = (2 * l) / d + 1;
			int64_t pm = seal::util::exponentiate_uint_mod(mod_2N.reduce(p), m, mod_2N);
			tmp.clear();
			poly_add(
				tmp,
				powertable[static_cast<int64_t>(k) + delta * pm * d / 2 - l * pm],
				mod_t,
				us[4 - delta]
			);
			poly_sub(
				tmp,
				powertable[static_cast<int64_t>(k) + (delta + 1) * pm * d / 2 - l * pm],
				mod_t,
				us[3 - delta]
			);
			uint64_t factor = seal::util::barrett_reduce_64(entry.second, mod_t);
			if (delta % 2 != 0) {
				factor = seal::util::negate_uint_mod(factor, mod_t);
			}
			poly_add(
				result,
				tmp,
				mod_t,
				factor
			);
		}
		poly_scale(result, global_factor, mod_t);
		results.push_back(std::move(result));
	}

	return results;
}

CompiledLinearTransform CompiledLinearTransform::scalar_slots_to_first_coefficients(std::shared_ptr<const SlotRing> slot_ring)
{
	const NegacyclicPowerTable slot_powertable(slot_ring->slot(), slot_ring->slot().generator(), slot_ring->N());
	auto basis_transform_matrix = [&slot_ring, &slot_powertable](std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& block, size_t block_row, size_t block_col) {
		const size_t power_of_x = block_col;
		const size_t out_slot = block_row;
		const size_t power_of_zeta = seal::util::multiply_uint_mod(
			power_of_x,
			inv_mod(std::get<0>(slot_ring->rotate(out_slot).galois_elements()), slot_ring->index_mod()),
			slot_ring->index_mod()
		);
		// this entry in the slot out_slot corresponds to the coset of x^power_of_x
		const poly entry = slot_powertable[power_of_zeta];
		for (size_t i = 0; i < slot_ring->slot_rank(); ++i) {
			if (entry[i] != 0) {
				block[std::make_tuple(i, 0)] = entry[i];
			}
		}
	};
	return CompiledLinearTransform::compile_slot_basis(slot_ring, basis_transform_matrix);
}

CompiledLinearTransform CompiledLinearTransform::first_coefficients_to_scalar_slots(std::shared_ptr<const SlotRing> slot_ring)
{
	auto basis_transform_matrix = [&slot_ring](std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& block, size_t block_row, size_t block_col) {
		for (size_t k = 0; k < slot_ring->slot_rank(); ++k) {
			const size_t power = seal::util::sub_uint_mod(
				block_row,
				seal::util::multiply_uint_mod(
					std::get<0>(slot_ring->rotate(block_col).galois_elements()),
					k,
					slot_ring->index_mod()
				),
				slot_ring->index_mod()
			);
			uint64_t value;
			if (power >= slot_ring->N()) {
				value = seal::util::negate_uint_mod(slot_ring->slot_one(block_col)[power - slot_ring->N()], slot_ring->R().scalar_mod);
			}
			else {
				value = slot_ring->slot_one(block_col)[power];
			}
			if (value != 0) {
				block[std::make_tuple(0, k)] = value;
			}
		}
	};
	return CompiledLinearTransform::compile_slot_basis(slot_ring, basis_transform_matrix);
}

CompiledLinearTransform CompiledLinearTransform::load_binary(std::shared_ptr<const SlotRing> slot_ring, std::istream& in)
{
	std::ios_base::sync_with_stdio(false);
	uint64_t version;
	uint64_t g1_subgroup_order;
	uint64_t g2_subgroup_order;
	in.read(reinterpret_cast<char*>(&version), sizeof(uint64_t));
	if (version != 3) {
		throw std::invalid_argument("wrong version");
	}
	in.read(reinterpret_cast<char*>(&g1_subgroup_order), sizeof(uint64_t));
	in.read(reinterpret_cast<char*>(&g2_subgroup_order), sizeof(uint64_t));
	assert(g1_subgroup_order == slot_ring->g1_ord());
	assert(g2_subgroup_order == 1 || g2_subgroup_order == 2);

	const size_t poly_len = slot_ring->N();
	const size_t coeff_count = g1_subgroup_order * g2_subgroup_order;
	std::vector<poly> coefficients;
	for (size_t i = 0; i < coeff_count; ++i) {
		coefficients.emplace_back();
		coefficients.back().resize(poly_len);
		in.read(reinterpret_cast<char*>(&coefficients.back()[0]), poly_len * sizeof(uint64_t));
	}
	return CompiledLinearTransform(slot_ring, std::move(coefficients));
}

void CompiledLinearTransform::save_binary(std::ostream& stream) const
{
	std::ios_base::sync_with_stdio(false);
	const uint64_t version = 3;
	const uint64_t g1_subgroup_order = this->g1_subgroup_order();
	stream.write(reinterpret_cast<const char*>(&version), sizeof(uint64_t));
	stream.write(reinterpret_cast<const char*>(&g1_subgroup_order), sizeof(uint64_t));
	size_t g2_subgroup_order = this->g2_subgroup_order();

	stream.write(reinterpret_cast<const char*>(&g2_subgroup_order), sizeof(uint64_t));
	assert(g1_subgroup_order * g2_subgroup_order == coefficients.size());

	size_t poly_len = slot_ring->N();
	for (size_t i = 0; i < coefficients.size(); ++i) {
		assert(coefficients[i].size() == poly_len);
		stream.write(reinterpret_cast<const char*>(&coefficients[i][0]), poly_len * sizeof(uint64_t));
	}
}

CompiledSubringLinearTransform CompiledLinearTransform::in_ring()&&
{
	return CompiledSubringLinearTransform(std::move(*this), slot_ring);
}

NegacyclicPowerTable::NegacyclicPowerTable(const SlotRing::SubringView& ring, poly generator, uint64_t half_order_generator)
	: ring(ring), generator(), N(half_order_generator)
{
	// the "identity" does not have to be 1 - it can be 1 in the current slot, and 0 in the others
	poly current = poly_exponentiate_mod(generator, 2 * half_order_generator, ring.scalar_mod, ring.poly_mod);
	content.reserve(N);
	for (size_t i = 0; i < N; ++i) {
		poly next = poly_mul_mod(current, generator, ring.scalar_mod, ring.poly_mod);
		content.push_back(std::move(current));
		current = std::move(next);
	}
	this->generator = std::move(generator);
}

poly NegacyclicPowerTable::operator[](int64_t i) const
{
	int64_t index = i % static_cast<int64_t>(2 * N);
	if (index < 0) {
		index += 2 * N;
	}
	if (index < static_cast<int64_t>(N)) {
		return content[index];
	}
	else {
		poly result = content[index - N];
		poly_scale(result, seal::util::negate_uint_mod(1, ring.scalar_mod), ring.scalar_mod);
		return result;
	}
}

void test_first_coeffs_to_scalar_slots()
{
	std::shared_ptr<SlotRing> slot_ring = std::make_shared<SlotRing>(p_257_test_parameters()->power_x_subring(7));
	std::shared_ptr<SlotRing> samller_slot_ring = std::make_shared<SlotRing>(slot_ring->power_x_subring(1));
	CompiledLinearTransform t = CompiledLinearTransform::first_coefficients_to_scalar_slots(samller_slot_ring);
	CompiledSubringLinearTransform transform(std::move(t), slot_ring);

	poly a = { 1, 0, 2, 0, 4, 0 };
	a.resize(slot_ring->N());
	poly b = transform(a);
	assert(slot_ring->extract_slot_value(b, 0)[0] == 1);
	assert(slot_ring->extract_slot_value(b, 1)[0] == 2);
	assert(slot_ring->extract_slot_value(b, 2)[0] == 4);
	for (size_t i = 3; i < slot_ring->slot_group_len(); ++i) {
		assert(slot_ring->extract_slot_value(b, i)[0] == 0);
	}
}

void test_compile_slot_basis()
{
	std::shared_ptr<SlotRing> slot_ring = small_test_parameters();
	auto matrix = [&slot_ring](std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& slotwise_matrix, size_t i, size_t j) {

		if (i == j) {
			for (size_t l = 0; l < slot_ring->slot_rank(); ++l) {
				slotwise_matrix[std::make_tuple(l, l)] = 1;
			}
		}
		else if (i == j + 16 && i < 32) {
				slotwise_matrix[std::make_tuple(0, 0)] = 1;
		}
		else if (j == i + 16 && j < 32) {
			slotwise_matrix[std::make_tuple(1, 1)] = 1;
		}
	};
	CompiledLinearTransform t = CompiledLinearTransform::compile_slot_basis(slot_ring, matrix);
	CompiledSubringLinearTransform transform = std::move(t).in_ring();

	poly a = slot_ring->from_slot_value({ 1, 1 }, 0);
	poly_add(a, slot_ring->from_slot_value({ 2, 0, 3 }, 1), slot_ring->R().scalar_mod);
	poly_add(a, slot_ring->from_slot_value({ 4, 0, 1 }, 16), slot_ring->R().scalar_mod);
	poly_add(a, slot_ring->from_slot_value({ 0, 1 }, 17), slot_ring->R().scalar_mod);

	poly result_poly = transform(a);

	poly expected;
	poly_add(expected, slot_ring->from_slot_value({ 1, 1 }, 0), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 2, 1, 3 }, 1), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 5, 0, 1 }, 16), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 2, 1 }, 17), slot_ring->R().scalar_mod);

	assert(result_poly == expected);
	std::cout << "test_compile_slot_basis(): success" << std::endl;
}

void test_apply_ciphertext()
{
	using namespace seal;

	std::shared_ptr<SlotRing> slot_ring = small_test_parameters();
	auto matrix = [&slot_ring](std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& slotwise_matrix, size_t i, size_t j) {
		
		if (i == j) {
			for (size_t l = 0; l < slot_ring->slot_rank(); ++l) {
				slotwise_matrix[std::make_tuple(l, l)] = 1;
			}
		}
		else if (i == j + 16 && i < 32) {
			slotwise_matrix[std::make_tuple(0, 0)] = 1;
		}
		else if (j == i + 16 && j < 32) {
			slotwise_matrix[std::make_tuple(1, 1)] = 1;
		}
	};
	CompiledLinearTransform t = CompiledLinearTransform::compile_slot_basis(slot_ring, matrix);
	CompiledSubringLinearTransform transform = std::move(t).in_ring();

	EncryptionParameters parms(scheme_type::bfv);
	parms.set_poly_modulus_degree(slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::Create(slot_ring->N(), { 40, 40 }));
	parms.set_plain_modulus(slot_ring->prime());
	SEALContext context(parms, false, sec_level_type::none);

	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	GaloisKeys gk;
	keygen.create_galois_keys(transform.galois_elements(), gk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);

	poly a = slot_ring->from_slot_value({ 1, 1 }, 0);
	poly_add(a, slot_ring->from_slot_value({ 2, 0, 3 }, 1), slot_ring->R().scalar_mod);
	poly_add(a, slot_ring->from_slot_value({ 4, 0, 1 }, 16), slot_ring->R().scalar_mod);
	poly_add(a, slot_ring->from_slot_value({ 0, 1 }, 17), slot_ring->R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(a) };

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);

	Ciphertext result_enc;
	transform.apply_ciphertext(x_enc, evaluator, gk, result_enc);

	Plaintext result;
	decryptor.decrypt(result_enc, result);

	poly result_poly(result.data(), result.data() + result.coeff_count());
	result_poly.resize(slot_ring->N());
	poly expected;
	poly_add(expected, slot_ring->from_slot_value({ 1, 1 }, 0), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 2, 1, 3 }, 1), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 5, 0, 1 }, 16), slot_ring->R().scalar_mod);
	poly_add(expected, slot_ring->from_slot_value({ 2, 1 }, 17), slot_ring->R().scalar_mod);

	assert(result_poly == expected);
	std::cout << "test_apply_ciphertext(): success" << std::endl;
}

void test_apply_ciphertext_subring()
{
	using namespace seal;

	std::shared_ptr<SlotRing> large_slot_ring = small_test_parameters();
	std::shared_ptr<SlotRing> slot_ring = std::make_shared<SlotRing>(large_slot_ring->power_x_subring(2));
	auto matrix = [&slot_ring](std::unordered_map<std::tuple<size_t, size_t>, uint64_t>& slotwise_matrix, size_t i, size_t j) {
		if (i == j) {
			for (size_t l = 0; l < 2; ++l) {
				slotwise_matrix[std::make_tuple(l, l)] = 1;
			}
		}
		else if (i == j + 16) {
			slotwise_matrix[std::make_tuple(0, 0)] = 1;
		}
		else if (j == i + 16) {
			slotwise_matrix[std::make_tuple(1, 1)] = 1;
		}
	};
	CompiledLinearTransform t = CompiledLinearTransform::compile_slot_basis(slot_ring, matrix);
	CompiledSubringLinearTransform transform(std::move(t), large_slot_ring);

	EncryptionParameters parms(scheme_type::bfv);
	parms.set_poly_modulus_degree(large_slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::Create(large_slot_ring->N(), { 40, 40 }));
	parms.set_plain_modulus(large_slot_ring->prime());
	SEALContext context(parms, false, sec_level_type::none);

	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	GaloisKeys gk;
	keygen.create_galois_keys(transform.galois_elements(), gk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);

	poly a = large_slot_ring->from_slot_value({ 1, 0, 0, 0, 1 }, 0);
	poly_add(a, large_slot_ring->from_slot_value({ 2 }, 1), large_slot_ring->R().scalar_mod);
	poly_add(a, large_slot_ring->from_slot_value({ 4 }, 16), large_slot_ring->R().scalar_mod);
	poly_add(a, large_slot_ring->from_slot_value({ 0, 0, 0, 0, 1 }, 17), large_slot_ring->R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(a) };

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);

	Ciphertext result_enc;
	transform.apply_ciphertext(x_enc, evaluator, gk, result_enc);

	Plaintext result;
	decryptor.decrypt(result_enc, result);

	poly result_poly(result.data(), result.data() + result.coeff_count());
	result_poly.resize(large_slot_ring->N());
	poly expected;
	poly_add(expected, large_slot_ring->from_slot_value({ 1, 0, 0, 0, 1 }, 0), large_slot_ring->R().scalar_mod);
	poly_add(expected, large_slot_ring->from_slot_value({ 2, 0, 0, 0, 1 }, 1), large_slot_ring->R().scalar_mod);
	poly_add(expected, large_slot_ring->from_slot_value({ 5 }, 16), large_slot_ring->R().scalar_mod);
	poly_add(expected, large_slot_ring->from_slot_value({ 2, 0, 0, 0, 1 }, 17), large_slot_ring->R().scalar_mod);

	assert(result_poly == expected);
	std::cout << "test_apply_ciphertext_subring(): success" << std::endl;
}

SlotRing::RawAuto CompiledSubringLinearTransform::automorphism(size_t index) const
{
	assert(slot_ring->slot_group_len() == subring_transform.slot_ring->slot_group_len());
	const size_t g1_index = index % subring_transform.g1_subgroup_order();
	const size_t g2_index = (index % subring_transform.coefficients.size()) / subring_transform.g1_subgroup_order();
	return slot_ring->raw_auto(g1_index * (subring_transform.slot_ring->g1_ord() / subring_transform.g1_subgroup_order()), g2_index);
}

SlotRing::RawAuto CompiledSubringLinearTransform::difference_automorphism(size_t from, size_t to) const
{
	const size_t from_g1_index = (from % subring_transform.g1_subgroup_order());
	const size_t from_g2_index = (from % subring_transform.coefficients.size()) / subring_transform.g1_subgroup_order();
	const size_t to_g1_index = (to % subring_transform.g1_subgroup_order());
	const size_t to_g2_index = (to % subring_transform.coefficients.size()) / subring_transform.g1_subgroup_order();
	size_t diff_g1_index = to_g1_index + subring_transform.g1_subgroup_order() - from_g1_index;
	if (diff_g1_index >= subring_transform.g1_subgroup_order()) {
		diff_g1_index -= subring_transform.g1_subgroup_order();
	}
	size_t diff_g2_index = to_g2_index + subring_transform.g2_subgroup_order() - from_g2_index;
	if (diff_g2_index >= subring_transform.g2_subgroup_order()) {
		diff_g2_index -= subring_transform.g2_subgroup_order();
	}
	return automorphism(diff_g1_index + diff_g2_index * subring_transform.g1_subgroup_order());
}

SlotRing::RawAuto CompiledSubringLinearTransform::reverse_automorphism(size_t index) const
{
	return difference_automorphism(index, 0);
}

size_t CompiledSubringLinearTransform::babystep_automorphism_count() const
{
	return subring_transform.babystep_automorphism_count();
}

size_t CompiledSubringLinearTransform::giantstep_automorphism_count() const
{
	return subring_transform.giantstep_automorphism_count();
}

CompiledSubringLinearTransform::CompiledSubringLinearTransform(CompiledLinearTransform&& transform, std::shared_ptr<const SlotRing> new_ring) : slot_ring(new_ring), subring_transform(std::move(transform))
{
}

CompiledSubringLinearTransform CompiledSubringLinearTransform::slots_to_coeffs(std::shared_ptr<const SlotRing> slot_ring)
{
	assert(slot_ring->is_d_half_sparse());
	if (slot_ring->is_d_sparse()) {
		std::shared_ptr<SlotRing> reduced_slot_ring = std::make_shared<SlotRing>(slot_ring->power_x_subring(log2_exact(slot_ring->slot_rank())));
		CompiledLinearTransform t = CompiledLinearTransform::scalar_slots_to_first_coefficients(reduced_slot_ring);
		return CompiledSubringLinearTransform(std::move(t), slot_ring);
	}
	else {
		std::shared_ptr<SlotRing> reduced_slot_ring = std::make_shared<SlotRing>(slot_ring->power_x_subring(log2_exact(slot_ring->slot_rank() / 2)));
		CompiledLinearTransform t = CompiledLinearTransform::scalar_slots_to_first_coefficients(reduced_slot_ring);
		return CompiledSubringLinearTransform(std::move(t), slot_ring);
	}
}

CompiledSubringLinearTransform CompiledSubringLinearTransform::coeffs_to_slots(std::shared_ptr<const SlotRing> slot_ring)
{
	assert(slot_ring->is_d_half_sparse());
	if (slot_ring->is_d_sparse()) {
		std::shared_ptr<SlotRing> reduced_slot_ring = std::make_shared<SlotRing>(slot_ring->power_x_subring(log2_exact(slot_ring->slot_rank())));
		CompiledLinearTransform t = CompiledLinearTransform::first_coefficients_to_scalar_slots(reduced_slot_ring);
		return CompiledSubringLinearTransform(std::move(t), slot_ring);
	}
	else {
		std::shared_ptr<SlotRing> reduced_slot_ring = std::make_shared<SlotRing>(slot_ring->power_x_subring(log2_exact(slot_ring->slot_rank() / 2)));
		CompiledLinearTransform t = CompiledLinearTransform::first_coefficients_to_scalar_slots(reduced_slot_ring);
		return CompiledSubringLinearTransform(std::move(t), slot_ring);
	}
}

poly CompiledSubringLinearTransform::operator()(const poly& x) const
{
	std::vector<poly> precomputed_values;
	precomputed_values.emplace_back(x);
	for (size_t i = 1; i < babystep_automorphism_count(); ++i) {
		// starting to remove lower digits has the effect of reusing rotations and recomputing frobenius,
		// which is faster since rotations take 2 key-switches
		size_t base_element_index = i - ((size_t)1 << highest_dividing_power2(i));
		assert(base_element_index < i);
		SlotRing::RawAuto automorphism_to_apply = difference_automorphism(base_element_index, i);
		precomputed_values.push_back(automorphism_to_apply(precomputed_values[base_element_index]));
	}

	poly result;
	for (size_t i = 0; i < giantstep_automorphism_count(); ++i) {
		poly current;
		for (size_t j = 0; j < babystep_automorphism_count(); ++j) {
			if (!is_zero(subring_transform.coefficients[i * babystep_automorphism_count() + j])) {
				poly coeff;
				slot_ring->from_power_x_subring(*subring_transform.slot_ring, subring_transform.coefficients[i * babystep_automorphism_count() + j], coeff);
				poly addition = poly_mul_mod(
					coeff,
					precomputed_values[j],
					slot_ring->R().scalar_mod,
					slot_ring->R().poly_mod
				);
				poly_add(
					current,
					addition,
					slot_ring->R().scalar_mod
				);
			}
		}
		SlotRing::RawAuto automorphism_to_apply = automorphism(i * babystep_automorphism_count());
		poly_add(
			result,
			automorphism_to_apply(current),
			slot_ring->R().scalar_mod
		);
	}
	return result;
}

void CompiledSubringLinearTransform::apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const
{
	if (coefficients_plain.size() == 0) {
		poly coeff;
		for (const auto& c : subring_transform.coefficients) {
			slot_ring->from_power_x_subring(*subring_transform.slot_ring, c, coeff);
			coefficients_plain.push_back(seal::Plaintext(coeff));
		}
	}
	else {
		assert(coefficients_plain.size() == subring_transform.coefficients.size());
	}

	std::vector<seal::Ciphertext> precomputed_values;
	precomputed_values.emplace_back(in);
	seal::Ciphertext tmp;
	for (size_t i = 1; i < babystep_automorphism_count(); ++i) {
		// starting to remove lower digits has the effect of reusing rotations and recomputing frobenius,
		// which is faster since rotations take 2 key-switches
		size_t base_element_index = i - ((size_t)1 << highest_dividing_power2(i));
		assert(base_element_index < i);
		SlotRing::RawAuto automorphism_to_apply = difference_automorphism(base_element_index, i);
		automorphism_to_apply.apply_ciphertext(precomputed_values[base_element_index], eval, gk, tmp);
		precomputed_values.push_back(tmp);
	}
	bool is_result_set = false;
	for (size_t i = 0; i < giantstep_automorphism_count(); ++i) {
		seal::Ciphertext current;
		bool is_current_set = false;

		for (size_t j = 0; j < babystep_automorphism_count(); ++j) {
			if (!is_zero(subring_transform.coefficients[i * babystep_automorphism_count() + j])) {
				if (is_current_set) {
					eval.multiply_plain(
						precomputed_values[j],
						coefficients_plain[i * babystep_automorphism_count() + j],
						tmp
					);
					eval.add_inplace(
						current,
						tmp
					);
				}
				else {
					eval.multiply_plain(
						precomputed_values[j],
						coefficients_plain[i * babystep_automorphism_count() + j],
						current
					);
					is_current_set = true;
				}
			}
		}
		if (is_current_set) {
			SlotRing::RawAuto automorphism_to_apply = automorphism(i * babystep_automorphism_count());
			automorphism_to_apply.apply_ciphertext(current, eval, gk, tmp);
			if (!is_result_set) {
				result = tmp;
				is_result_set = true;
			}
			else {
				eval.add_inplace(result, tmp);
			}
		}
	}
}

std::vector<uint32_t> CompiledSubringLinearTransform::galois_elements() const
{
	std::vector<uint32_t> result; 
	for (size_t i = 1; i < babystep_automorphism_count(); ++i) {
		// starting to remove lower digits has the effect of reusing rotations and recomputing frobenius,
		// which is faster since rotations take 2 key-switches
		size_t base_element_index = i - ((size_t)1 << highest_dividing_power2(i));
		assert(base_element_index < i);
		SlotRing::RawAuto automorphism_to_apply = difference_automorphism(base_element_index, i);
		result.push_back(automorphism_to_apply.galois_element());
	}
	for (size_t i = 1; i < giantstep_automorphism_count(); ++i) {
		const auto automorphism = this->automorphism(i * babystep_automorphism_count());
		result.push_back(automorphism.galois_element());
	}
	return result;
}

CompiledLinearTransform&& CompiledSubringLinearTransform::transform() &&
{
	return std::move(this->subring_transform);
}

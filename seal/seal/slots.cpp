#include "slots.h"

size_t kswitch_counter = 0;

void log_relin()
{
	++kswitch_counter;
}

void log_galois()
{
	++kswitch_counter;
}

poly SlotRing::SubringView::generator()
{
	poly result = { 0, 1 };
	poly_reduce_mod(result, scalar_mod, poly_mod);
	return result;
}

size_t order_mod_2N(uint64_t x, const seal::Modulus& mod_2N)
{
	uint64_t current = seal::util::barrett_reduce_64(x, mod_2N);
	for (size_t i = 0; i < log2_exact(mod_2N.value()); ++i) {
		if (current == 1) {
			return (uint64_t)1 << i;
		}
		else {
			current = seal::util::multiply_uint_mod(current, current, mod_2N);
		}
	}
	assert(false);
	return 0;
}

bool is_in_subgroup(uint64_t x, uint64_t gen, size_t order, const seal::Modulus& mod_2N)
{
	uint64_t current = 1;
	seal::util::MultiplyUIntModOperand mul_gen;
	mul_gen.set(mod_2N.reduce(gen), mod_2N);
	for (size_t i = 0; i < order; ++i) {
		if (current == x) {
			return true;
		}
		current = seal::util::multiply_uint_mod(current, mul_gen, mod_2N);
	}
	return false;
}

bool is_group_cyclic(uint64_t gen, size_t order, size_t ord_p, uint64_t p, const seal::Modulus& mod_2N)
{
	uint64_t gen_pow_n_half = seal::util::exponentiate_uint_mod(gen, order / 2, mod_2N);
	return !is_in_subgroup(gen_pow_n_half, p, ord_p, mod_2N);
}

void SlotRing::init_slot_group()
{
	d = order_mod_2N(p, mod_2N);
	n = N() / d;
	slot_group_cyclic = is_group_cyclic(g1, n, d, p, mod_2N);
}

void SlotRing::init_p_log()
{
	seal::util::MultiplyUIntModOperand mul_g1;
	mul_g1.set(g1, mod_2N);
	uint64_t current = 1;
	for (size_t i = 0; i < 2; ++i) {
		for (size_t j = 0; j < N() / 2; ++j) {
			if (current == mod_2N.reduce(p)) {
				p_log = std::make_tuple(j, i);
				return;
			}
			current = seal::util::multiply_uint_mod(current, mul_g1, mod_2N);
		}
		current = seal::util::multiply_uint_mod(current, g2, mod_2N);
	}
	assert(false);
}

void SlotRing::init_slot_modulus(poly base_unit_vector)
{
	assert(base_unit_vector.size() == N());
	if (slot_group_cyclic) {
		poly current = base_unit_vector;
		for (size_t i = 0; i < n; ++i) {
			unit_vectors.push_back(current);
			g1_automorphism(current, 1);
		}
	}
	else {
		poly current = base_unit_vector;
		for (size_t i = 0; i < n / 2; ++i) {
			unit_vectors.push_back(current);
			g1_automorphism(current, 1);
		}
		current = base_unit_vector;
		g2_automorphism(current);
		for (size_t i = 0; i < n / 2; ++i) {
			unit_vectors.push_back(current);
			g1_automorphism(current, 1);
		}
	}
	poly f;
	f.resize(N() + 1);
	f[0] = 1;
	f[N()] = 1;

	poly a = std::move(base_unit_vector);
	a[0] = seal::util::sub_uint_mod(a[0], 1, mod_pe);
	poly factor = std::get<2>(poly_eea(a, f, mod_pe));
	assert(factor[d] != 0);
	uint64_t lc_inv = inv_mod(factor[d], mod_pe);
	poly_scale(factor, seal::util::negate_uint_mod(lc_inv, mod_pe), mod_pe);
	assert(factor[d] == seal::util::negate_uint_mod(1, mod_pe));
	factor[d] = 0;
	poly_normalize(factor);
	// this should already be guaranteed, but check again
	assert(is_poly_sparse(factor, d / 2));
	slot_modulus.n = d;
	slot_modulus.x_power_n = std::move(factor);
}

void SlotRing::g1_automorphism(poly& x, size_t iters, const seal::Modulus* mod) const
{
	if (iters == 0) {
		return;
	}
	if (mod == nullptr) {
		mod = &mod_pe;
	}
	iters = iters % g1_ord();
	assert(x.size() == N());

#ifdef CONTRACT_TEST
	poly initial = x;
#endif

	unsigned long log_power_2_iters = highest_dividing_power2(iters);
	uint64_t power_2_iters = (uint64_t)1 << log_power_2_iters;
	uint64_t odd_iters = iters / power_2_iters;

	seal::util::MultiplyUIntModOperand g_mul;
	g_mul.set(seal::util::exponentiate_uint_mod(g1, odd_iters, mod_2N), mod_2N);
	seal::util::MultiplyUIntModOperand step_mul;
	step_mul.set(seal::util::exponentiate_uint_mod(g1, iters, mod_2N), mod_2N);

	const auto swap_entries = [&](uint64_t fst, uint64_t snd) {
		bool fst_red = fst >= N();
		bool snd_red = snd >= N();
		if (fst_red) {
			fst -= N();
		}
		if (snd_red) {
			snd -= N();
		}
		std::swap(x[fst], x[snd]);
		if (fst_red != snd_red) {
			if (fst != snd) {
				x[fst] = seal::util::negate_uint_mod(x[fst], *mod);
				x[snd] = seal::util::negate_uint_mod(x[snd], *mod);
			}
			else {
				x[fst] = seal::util::negate_uint_mod(x[fst], *mod);
			}
		}
	};

	for (size_t i = 0; i <= log2N - 2; ++i) {
		uint64_t base = (size_t)1 << i;
		// the log of number of entries until the sequence (base * (g1^iters)^i) % N starts repeating
		int64_t log_step_order;
		bool has_second_generator;
		// this strange logic here is because (Z/2Z)* and (Z/4Z)* are cyclic, but afterwards not anymore
		if (log2N <= i + 3) {
			log_step_order = log2N - 2 - i - log_power_2_iters;
			has_second_generator = false;
		}
		else {
			log_step_order = log2N - 3 - i - log_power_2_iters;
			has_second_generator = true;
		}

		if (log_step_order <= 0) {
			// in this case, have that X -> +/- X
			for (size_t k = 0; k < power_2_iters; ++k) {
				if ((seal::util::multiply_uint_mod(base, step_mul, mod_2N) >= N()) != (base >= N())) {
					x[base % N()] = seal::util::negate_uint_mod(x[base % N()], *mod);
				}
				if (has_second_generator) {
					base = seal::util::multiply_uint_mod(base, g2, mod_2N);
					if ((seal::util::multiply_uint_mod(base, step_mul, mod_2N) >= N()) != (base >= N())) {
						x[base % N()] = seal::util::negate_uint_mod(x[base % N()], *mod);
					}
				}
				base = seal::util::multiply_uint_mod(base, g_mul, mod_2N);
			}
		} else {
			for (size_t k = 0; k < power_2_iters; ++k) {
				// otherwise, we cycle through the sequence (base * (g1^iters)^i) % N
				// and swap/negate coefficients
				uint64_t current = base;
				for (size_t j = 0; j < ((size_t)1 << log_step_order); ++j) {
					current = seal::util::multiply_uint_mod(current, step_mul, mod_2N);
					swap_entries(base, current);
				}
				if (has_second_generator) {
					base = seal::util::multiply_uint_mod(base, g2, mod_2N);
					current = base;
					for (size_t j = 0; j < ((size_t)1 << log_step_order); ++j) {
						current = seal::util::multiply_uint_mod(current, step_mul, mod_2N);
						swap_entries(base, current);
					}
				}
				base = seal::util::multiply_uint_mod(base, g_mul, mod_2N);
			}
		}
	}

#ifdef CONTRACT_TEST
	for (uint64_t i = 0; i < N(); ++i) {
		uint64_t index = seal::util::multiply_uint_mod(i, step_mul, mod_2N);
		if (index < N()) {
			assert(initial[i] == x[index]);
		}
		else {
			assert(initial[i] == seal::util::negate_uint_mod(x[index - N()], *mod));
		}
	}
#endif
}

void SlotRing::g2_automorphism(poly& x, const seal::Modulus* mod) const
{
	if (mod == nullptr) {
		mod = &mod_pe;
	}
#ifdef CONTRACT_TEST
	poly initial = x;
#endif

	assert(x.size() == N());
	// we heavily use that g2 has order 2 in (Z/2NZ)*
	uint64_t partner = 0;
	for (uint64_t base = 0; base < N(); ++base) {
		assert(base == 0 || base != partner);
		assert(partner == seal::util::multiply_uint_mod(base, g2, mod_2N));
		assert(base == seal::util::multiply_uint_mod(partner, g2, mod_2N));
		if (partner >= N() && partner - N() == base) {
			x[base] = seal::util::negate_uint_mod(x[base], *mod);
		} else if (partner >= N() && base < (partner - N())) {
			std::swap(x[base], x[partner - N()]);
			x[base] = seal::util::negate_uint_mod(x[base], *mod);
			x[partner - N()] = seal::util::negate_uint_mod(x[partner - N()], *mod);
		} else if (partner < N() && base < partner) {
			std::swap(x[base], x[partner]);
		}
		partner = seal::util::add_uint_mod(partner, g2, mod_2N);
	}

#ifdef CONTRACT_TEST
	for (uint64_t i = 0; i < N(); ++i) {
		uint64_t index = seal::util::multiply_uint_mod(i, g2, mod_2N);
		if (index < N()) {
			assert(initial[i] == x[index]);
		}
		else {
			assert(initial[i] == seal::util::negate_uint_mod(x[index - N()], *mod));
		}
	}
#endif
}

void SlotRing::apply_frobenius(poly& x, size_t iters, const seal::Modulus* mod) const
{
	g1_automorphism(x, (iters * std::get<0>(p_log)) % g1_ord(), mod);
	if ((iters * std::get<1>(p_log)) % 2 != 0) {
		g2_automorphism(x);
	}
}

SlotRing::SlotRing(
	uint64_t log2N,
	uint64_t p,
	uint64_t e,
	uint64_t g2,
	poly base_unit_vector
) :
	log2N(log2N), p(p), e(e), mod_pe(seal::util::exponentiate_uint(p, e)), mod_2N((uint64_t)1 << log2N), g1(3), g2(g2), mod_2N_cyclotomic()
{
	mod_2N_cyclotomic.n = N();
	mod_2N_cyclotomic.x_power_n = { seal::util::negate_uint_mod(1, mod_pe) };
	init_slot_group();
	assert(order_mod_2N(g1, mod_2N) == N() / 2);
	assert(order_mod_2N(g2, mod_2N) == 2);
	d_sparse = is_poly_sparse(base_unit_vector, d);

	// assert that <g1, g2> is the whole index group
	assert(!is_in_subgroup(g2, g1, N() / 2, mod_2N));
	init_p_log();

	// currently, we assume that base_unit_vector is a polynomial in X^(d/2)
	assert(is_poly_sparse(base_unit_vector, d / 2));
	init_slot_modulus(std::move(base_unit_vector));
}

SlotRing::SlotRing(
	uint64_t log2N,
	uint64_t p,
	uint64_t e,
	poly base_unit_vector
) : SlotRing(log2N, p, e, ((uint64_t)1 << (log2N - 1)) - 1, base_unit_vector)
{}

poly SlotRing::extract_slot_value(poly x, size_t slot) const
{
	if (slot_group_cyclic) {
		const size_t shift = slot % n;
		g1_automorphism(x, g1_ord() - shift);
	}
	else {
		const size_t g1_shift = slot % (n / 2);
		const size_t g2_shift = (slot / (n / 2)) % 2;
		g1_automorphism(x, g1_ord() - g1_shift);
		if (g2_shift != 0) {
			g2_automorphism(x);
		}
	}
	poly_reduce_mod(x, mod_pe, slot_modulus);
	return x;
}

poly SlotRing::from_slot_value(const poly& x, size_t slot) const
{
	poly result = poly_mul_mod(x, unit_vectors[0], mod_pe, mod_2N_cyclotomic);
	if (slot_group_cyclic) {
		const size_t shift = slot % n;
		g1_automorphism(result, shift);
	}
	else {
		const size_t g1_shift = slot % (n / 2);
		const size_t g2_shift = (slot / (n / 2)) % 2;
		g1_automorphism(result, g1_shift);
		if (g2_shift != 0) {
			g2_automorphism(result);
		}
	}
#ifdef CONTRACT_TEST
	poly check = extract_slot_value(result, slot);
	poly_normalize(check);
	poly in = x;
	poly_reduce_mod(in, this->slot().scalar_mod, this->slot().poly_mod);
	poly_normalize(in);
	assert(check == in);
#endif
	return result;
}

const poly& SlotRing::slot_one(size_t slot) const noexcept
{
	return unit_vectors[slot % n];
}

uint64_t SlotRing::N() const noexcept
{
	return (uint64_t)1 << (log2N - 1);
}

bool SlotRing::is_poly_sparse(const poly& x_power_n, size_t d)
{
	if (d <= 1) {
		return true;
	}
	for (size_t i = 0; i < x_power_n.size(); ++i) {
		if ((x_power_n[i] != 0) && (i % d != 0)) {
			return false;
		}
	}
	return true;
}

bool SlotRing::is_d_half_sparse() const
{
	// we assume that this is always the case, and check it in the constructor via check_poly_sparse()
	return true;
}

bool SlotRing::is_d_sparse() const
{
	return d_sparse;
}

SlotRing SlotRing::power_x_subring(size_t index_log2) const
{
	const size_t index = (size_t)1 << index_log2;
	assert(is_d_half_sparse());
	assert(((size_t)1 << index_log2) <= d);
	//if (((size_t)1 << index_log2) > d / 2) {
	//	throw std::invalid_argument("Power-of-X subring is not trivially computable");
	//}
	poly new_base_unit_vector;
	new_base_unit_vector.resize(unit_vectors[0].size() / index);
	for (size_t i = 0; i < unit_vectors[0].size(); ++i) {
		if (i % index == 0) {
			new_base_unit_vector[i / index] = unit_vectors[0][i];
		}
		else {
			// this should already be guaranteed, but check anyway
			assert(unit_vectors[0][i] == 0);
		}
	}
	uint64_t new_index_modulus = (uint64_t)1 << (log2N - index_log2);
	SlotRing result(log2N - index_log2, p, e, g2 % new_index_modulus, std::move(new_base_unit_vector));

	poly check1;
	poly tmp = unit_vectors[0];
	g2_automorphism(tmp);
	in_power_x_subring(result, tmp, check1);
	poly check2 = result.unit_vectors[0];
	result.g2_automorphism(check2);
	assert(check1 == check2);

	return result;
}

void SlotRing::in_power_x_subring(const SlotRing& subring, const poly& el, poly& result) const
{
	assert(N() % subring.N() == 0);
	const size_t index = N() / subring.N();
	result.clear();
	result.resize(el.size() / index);
	for (size_t i = 0; i < el.size(); ++i) {
		if (i % index == 0) {
			result[i / index] = el[i];
		}
		else if (el[i] != 0) {
			throw std::invalid_argument("Given polynomial is not contained in the subring");
		}
	}
}

void SlotRing::from_power_x_subring(const SlotRing& subring, const poly& el, poly& result) const
{
	assert(N() % subring.N() == 0);
	const size_t index = N() / subring.N();
	result.clear();
	result.resize(el.size() * index);
	for (size_t i = 0; i < el.size(); ++i) {
		result[i * index] = el[i];
	}
}

SlotRing SlotRing::change_exponent(uint64_t new_exp) const
{
	assert(new_exp > 0);
	assert(new_exp <= e);
	seal::Modulus mod_new(seal::util::exponentiate_uint(p, new_exp));
	poly new_base_unit_vector;
	for (uint64_t x : unit_vectors[0]) {
		new_base_unit_vector.push_back(mod_new.reduce(x));
	}
	SlotRing result(log2N, p, new_exp, std::move(new_base_unit_vector));
	assert(result.d == d);
	assert(result.n == n);
	assert(result.p_log == p_log);
	return result;
}

SlotRing::Rotation SlotRing::block_rotate(size_t slot, size_t block_size) const
{
	assert(slot_group_len() % block_size == 0);
	poly forward_mask;

	// due to the order in which we store the unit vectors, this is also fine if !slot_group_cyclic
	size_t effective_block_size = block_size;
	if (!slot_group_cyclic && block_size == n) {
		effective_block_size = n / 2;
	}
	const size_t effective_block_count = slot_group_len() / effective_block_size;
	const size_t s = slot % effective_block_size;

	for (size_t i = 0; i < effective_block_count; ++i) {
		for (size_t j = 0; j < effective_block_size - s; ++j) {
			poly_add(forward_mask, unit_vectors[i * effective_block_size + j], mod_pe);
		}
	}
	poly backward_mask;
	for (size_t i = 0; i < effective_block_count; ++i) {
		for (size_t j = effective_block_size - s; j < effective_block_size; ++j) {
			poly_add(backward_mask, unit_vectors[i * effective_block_size + j], mod_pe);
		}
	}
	Rotation result(std::move(forward_mask), std::move(backward_mask), slot % block_size, *this, block_size);
	assert(result.effective_block_size() == effective_block_size);
	return result;
}

std::shared_ptr<SlotRing> small_test_parameters() {
	poly unit_vector;
	unit_vector.resize(512);
	auto values = { 1, 124, 4, 120, 11, 109, 29, 80, 76, 4, 72, 59, 13, 46, 94, 79, 15, 64, 78, 113, 92, 21, 71, 77, 121, 83, 38, 45, 120, 52, 68, 111, 84, 27, 57, 97, 87, 10, 77, 60, 17, 43, 101, 69, 32, 37, 122, 42, 80, 89, 118, 98, 20, 78, 69, 9, 60, 76, 111, 92, 19, 73, 73, 0, 73, 54, 19, 35, 111, 51, 60, 118, 69, 49, 20, 29, 118, 38, 80, 85, 122, 90, 32, 58, 101, 84, 17, 67, 77, 117, 87, 30, 57, 100, 84, 16, 68, 75, 120, 82, 38, 44, 121, 50, 71, 106, 92, 14, 78, 63, 15, 48, 94, 81, 13, 68, 72, 123, 76, 47, 29, 18, 11, 7, 4, 3, 1, 2 };
	auto it = values.end();
	for (size_t i = 0; i < 128; ++i) {
		--it;
		unit_vector[4 * i] = *it;
	}
	return std::make_shared<SlotRing>(10, 127, 1, std::move(unit_vector));
}

std::shared_ptr<SlotRing> small_p_257_test_parameters() {
	poly unit_vector;
	unit_vector.resize(512);
	auto values = { 142, 59, 77, 70, 87, 9, 125, 137, 218, 58, 6, 169, 177, 231, 210, 4, 27, 118, 154, 140, 174, 18, 250, 17, 179, 116, 12, 81, 97, 205, 163, 8, 54, 236, 51, 23, 91, 36, 243, 34, 101, 232, 24, 162, 194, 153, 69, 16, 108, 215, 102, 46, 182, 72, 229, 68, 202, 207, 48, 67, 131, 49, 138, 32, 216, 173, 204, 92, 107, 144, 201, 136, 147, 157, 96, 134, 5, 98, 19, 64, 175, 89, 151, 184, 214, 31, 145, 15, 37, 57, 192, 11, 10, 196, 38, 128, 93, 178, 45, 111, 171, 62, 33, 30, 74, 114, 127, 22, 20, 135, 76, 256, 186, 99, 90, 222, 85, 124, 66, 60, 148, 228, 254, 44, 40, 13, 152, 255 };
	auto it = values.end();
	for (size_t i = 0; i < 128; ++i) {
		--it;
		unit_vector[4 * i] = *it;
	}
	return std::make_shared<SlotRing>(10, 257, 1, std::move(unit_vector));
}

std::shared_ptr<SlotRing> p_127_test_parameters() {
	auto coefficients = { 1993024, 240875, 1207593, 404466, 20647, 348018, 811998, 401243, 1072206, 620779, 953549, 755478, 795511, 145718, 598024, 1753912, 1101795, 1868445, 1456233, 2048122, 1712610, 1283514, 1465554, 1560331, 67929, 1866101, 1288626, 1400697, 1138076, 534561, 784782, 480838, 1522444, 840918, 927040, 291915, 532443, 394526, 680455, 1554433, 1544666, 593498, 1618356, 1931381, 1791404, 1865968, 788282, 822881, 1782322, 190115, 1837873, 1709323, 1712134, 1718364, 451719, 239328, 644360, 1208324, 278256, 1940492, 1672965, 2046996, 1599175, 1599175, 1387, 1672965, 107891, 278256, 840059, 644360, 1809055, 451719, 330019, 1712134, 339060, 1837873, 1858268, 1782322, 1225502, 788282, 182415, 1791404, 117002, 1618356, 1454885, 1544666, 493950, 680455, 1653857, 532443, 1756468, 927040, 1207465, 1522444, 1567545, 784782, 1513822, 1138076, 647686, 1288626, 182282, 67929, 488052, 1465554, 764869, 1712610, 261, 1456233, 179938, 1101795, 294471, 598024, 1902665, 795511, 1292905, 953549, 1427604, 1072206, 1647140, 811998, 1700365, 20647, 1643917, 1207593, 1807508, 1993024, 32006 };
	auto indices = { 32512, 32256, 32000, 31744, 31488, 31232, 30976, 30720, 30464, 30208, 29952, 29696, 29440, 29184, 28928, 28672, 28416, 28160, 27904, 27648, 27392, 27136, 26880, 26624, 26368, 26112, 25856, 25600, 25344, 25088, 24832, 24576, 24320, 24064, 23808, 23552, 23296, 23040, 22784, 22528, 22272, 22016, 21760, 21504, 21248, 20992, 20736, 20480, 20224, 19968, 19712, 19456, 19200, 18944, 18688, 18432, 18176, 17920, 17664, 17408, 17152, 16896, 16640, 16128, 15872, 15616, 15360, 15104, 14848, 14592, 14336, 14080, 13824, 13568, 13312, 13056, 12800, 12544, 12288, 12032, 11776, 11520, 11264, 11008, 10752, 10496, 10240, 9984, 9728, 9472, 9216, 8960, 8704, 8448, 8192, 7936, 7680, 7424, 7168, 6912, 6656, 6400, 6144, 5888, 5632, 5376, 5120, 4864, 4608, 4352, 4096, 3840, 3584, 3328, 3072, 2816, 2560, 2304, 2048, 1792, 1536, 1280, 1024, 768, 512, 256, 0 };
	poly unit_vector;
	unit_vector.resize(1 << 15);
	auto index_it = indices.begin();
	auto coeff_it = coefficients.begin();
	for (; index_it != indices.end(); ++index_it, ++coeff_it) {
		unit_vector[*index_it] = *coeff_it;
	}
	return std::make_shared<SlotRing>(16, 127, 3, std::move(unit_vector));
}

std::shared_ptr<SlotRing> p_257_test_parameters()
{
	std::initializer_list<uint64_t> coefficients = { 7964572, 6261607, 16240935, 102099, 8906165, 9188787, 6652056, 6514830, 13560566, 4058859, 6803824, 2637246, 14964259, 8700452, 13885406, 12517446, 3289884, 7193548, 14438928, 14310157, 5483269, 14464749, 7681466, 12493558, 15081196, 16410594, 12121931, 2032437, 14600010, 10166868, 4700179, 13476574, 983850, 9275880, 13267162, 16806538, 13190359, 11592792, 7816384, 14932248, 13347910, 5088061, 14931724, 5704791, 15187352, 16367455, 12927940, 16583969, 13166475, 2450967, 2314130, 14807101, 16693874, 3439246, 7937160, 8764025, 7597636, 1868854, 15532357, 13962877, 1688878, 11856744, 13426589, 3975308, 8638757, 11867405, 4139446, 1645406, 2719681, 14214557, 15278080, 7294824, 6060464, 13542772, 10988902, 14604673, 12583496, 6989727, 13949465, 555184, 3553200, 14152308, 1356854, 588714, 6387435, 961982, 8029853, 7134849, 1050653, 14661907, 15973513, 16635621, 6987326, 13622481, 9026649, 3426195, 6011323, 13036503, 10929998, 396662, 3442686, 15482513, 13212403, 8408813, 11334031, 8307382, 3360145, 2311480, 6280586, 7839149, 11153619, 11748240, 11615815, 3330562, 14514936, 8478652, 16688123, 8347227, 3293264, 6218689, 2032504, 16192770, 1262895, 15951777, 3579793, 7381310, 2299274, 16841979 };
	auto indices = { 32512, 32256, 32000, 31744, 31488, 31232, 30976, 30720, 30464, 30208, 29952, 29696, 29440, 29184, 28928, 28672, 28416, 28160, 27904, 27648, 27392, 27136, 26880, 26624, 26368, 26112, 25856, 25600, 25344, 25088, 24832, 24576, 24320, 24064, 23808, 23552, 23296, 23040, 22784, 22528, 22272, 22016, 21760, 21504, 21248, 20992, 20736, 20480, 20224, 19968, 19712, 19456, 19200, 18944, 18688, 18432, 18176, 17920, 17664, 17408, 17152, 16896, 16640, 16384, 16128, 15872, 15616, 15360, 15104, 14848, 14592, 14336, 14080, 13824, 13568, 13312, 13056, 12800, 12544, 12288, 12032, 11776, 11520, 11264, 11008, 10752, 10496, 10240, 9984, 9728, 9472, 9216, 8960, 8704, 8448, 8192, 7936, 7680, 7424, 7168, 6912, 6656, 6400, 6144, 5888, 5632, 5376, 5120, 4864, 4608, 4352, 4096, 3840, 3584, 3328, 3072, 2816, 2560, 2304, 2048, 1792, 1536, 1280, 1024, 768, 512, 256, 0 };
	poly unit_vector;
	unit_vector.resize(1 << 15);
	auto index_it = indices.begin();
	auto coeff_it = coefficients.begin();
	for (; index_it != indices.end(); ++index_it, ++coeff_it) {
		unit_vector[*index_it] = *coeff_it;
	}
	return std::make_shared<SlotRing>(SlotRing(16, 257, 3, std::move(unit_vector)).change_exponent(2));
}

void test_g_automorphisms() {
	SlotRing r = *small_test_parameters();
	assert(r.slot_group_dims() == 1);
	assert(r.slot_group_len() == 64);
	assert(r.slot_rank() == 8);

	poly p;
	p.resize(512);
	p[1] = 1;
	r.g1_automorphism(p, 1);
	poly expected;
	expected.resize(512);
	expected[3] = 1;
	assert(expected == p);

	p.clear();
	p.resize(512);
	expected.clear();
	expected.resize(512);
	p[0] = 2;
	p[3] = 1;
	p[256] = 1;
	p[257] = 4;
	expected[0] = 2;
	expected[9] = 1;
	expected[256] = 126;
	expected[259] = 123;
	r.g1_automorphism(p, 1);
	assert(expected == p);

	p.clear();
	p.resize(512);
	expected.clear();
	expected.resize(512);
	p[1] = 1;
	p[2] = 2;
	p[3] = 4;
	p[5] = 8; // this is not in the group <g1>
	p[11] = 16; // this is in the group <g1>, but not in <g1^2>
	p[25] = 32; // this is in the group <g1^2>, but not in <g1^4>
	expected[497] = 126;
	expected[482] = 125;
	expected[467] = 123;
	expected[437] = 119;
	expected[347] = 111;
	expected[137] = 95;
	r.g1_automorphism(p, 3 * 4);
	assert(expected == p);

	std::cout << "test_g_automorphisms(): success" << std::endl;
}

void test_block_rotate()
{
	using namespace seal;

	SlotRing slot_ring = *small_test_parameters();
	EncryptionParameters parms(scheme_type::bfv);
	parms.set_poly_modulus_degree(slot_ring.N());
	parms.set_coeff_modulus(CoeffModulus::Create(slot_ring.N(), { 40, 40 }));
	parms.set_plain_modulus(slot_ring.R().scalar_mod.value());
	SEALContext context(parms, false, sec_level_type::none);

	SlotRing::Rotation rot = slot_ring.block_rotate(3, 32);

	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	GaloisKeys gk;
	keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);

	poly a = slot_ring.from_slot_value({ 1, 1 }, 0);
	poly_add(a, slot_ring.from_slot_value({ 2, 0, 3 }, 1), slot_ring.R().scalar_mod);
	poly_add(a, slot_ring.from_slot_value({ 4, 0, 1 }, 16), slot_ring.R().scalar_mod);
	poly_add(a, slot_ring.from_slot_value({ 0, 1 }, 17), slot_ring.R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(a) };

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);

	Ciphertext result_enc;
	rot.apply_ciphertext(x_enc, evaluator, gk, result_enc);

	Plaintext result;
	decryptor.decrypt(result_enc, result);

	poly result_poly(result.data(), result.data() + slot_ring.N());
	poly expected;
	poly_add(expected, slot_ring.from_slot_value({ 1, 1 }, 3), slot_ring.R().scalar_mod);
	poly_add(expected, slot_ring.from_slot_value({ 2, 0, 3 }, 4), slot_ring.R().scalar_mod);
	poly_add(expected, slot_ring.from_slot_value({ 4, 0, 1 }, 19), slot_ring.R().scalar_mod);
	poly_add(expected, slot_ring.from_slot_value({ 0, 1 }, 20), slot_ring.R().scalar_mod);

	assert(rot(a) == expected);
	assert(result_poly == expected);

	std::cout << "test_block_rotate(): success" << std::endl;
}

void test_rotate_noncyclic()
{
	using namespace seal;

	SlotRing slot_ring = *p_257_test_parameters();
	EncryptionParameters parms(scheme_type::bfv);
	parms.set_poly_modulus_degree(slot_ring.N());
	parms.set_coeff_modulus(CoeffModulus::Create(slot_ring.N(), { 40, 40, 40 }));
	parms.set_plain_modulus(slot_ring.R().scalar_mod.value());
	SEALContext context(parms, false, sec_level_type::none);


	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);

	poly a = slot_ring.from_slot_value({ 1, 1 }, 0);
	poly_add(a, slot_ring.from_slot_value({ 2, 0, 3 }, 63), slot_ring.R().scalar_mod);
	poly_add(a, slot_ring.from_slot_value({ 0, 1 }, 65), slot_ring.R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(a) };

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);
	{
		SlotRing::Rotation rot = slot_ring.rotate(3);

		GaloisKeys gk;
		keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk);

		Ciphertext result_enc;
		rot.apply_ciphertext(x_enc, evaluator, gk, result_enc);

		Plaintext result;
		decryptor.decrypt(result_enc, result);

		poly result_poly(result.data(), result.data() + slot_ring.N());
		poly expected;
		poly_add(expected, slot_ring.from_slot_value({ 1, 1 }, 3), slot_ring.R().scalar_mod);
		poly_add(expected, slot_ring.from_slot_value({ 2, 0, 3 }, 66), slot_ring.R().scalar_mod);
		poly_add(expected, slot_ring.from_slot_value({ 0, 1 }, 68), slot_ring.R().scalar_mod);

		assert(rot(a) == expected);
		assert(result_poly == expected);
	}
	{
		SlotRing::Rotation rot = slot_ring.rotate(67);

		GaloisKeys gk;
		keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk);

		Ciphertext result_enc;
		rot.apply_ciphertext(x_enc, evaluator, gk, result_enc);

		Plaintext result;
		decryptor.decrypt(result_enc, result);

		poly result_poly(result.data(), result.data() + slot_ring.N());
		poly expected;
		poly_add(expected, slot_ring.from_slot_value({ 1, 1 }, 67), slot_ring.R().scalar_mod);
		poly_add(expected, slot_ring.from_slot_value({ 2, 0, 3 }, 2), slot_ring.R().scalar_mod);
		poly_add(expected, slot_ring.from_slot_value({ 0, 1 }, 4), slot_ring.R().scalar_mod);

		assert(rot(a) == expected);
		assert(result_poly == expected);
	}
	std::cout << "test_rotate_noncyclic(): success" << std::endl;
}

size_t SlotRing::Rotation::effective_block_size() const
{
	if (block_size == slot_ring.n && !slot_ring.slot_group_cyclic) {
		return slot_ring.n / 2;
	}
	else {
		return block_size;
	}
}

poly SlotRing::Rotation::operator()(const poly& x) const
{
	if (s == 0) {
		return x;
	}
	poly a = poly_mul_mod(x, fmask, slot_ring.mod_pe, slot_ring.mod_2N_cyclotomic);
	poly b = poly_mul_mod(x, bmask, slot_ring.mod_pe, slot_ring.mod_2N_cyclotomic);
	assert(b.size() == slot_ring.N());
	slot_ring.g1_automorphism(a, s % effective_block_size());
	slot_ring.g1_automorphism(b, slot_ring.g1_ord() + s % effective_block_size() - effective_block_size());
	if (block_size == slot_ring.n && !slot_ring.slot_group_cyclic) {
		slot_ring.g2_automorphism(b);
	}
	poly_add(a, b, slot_ring.mod_pe);
	if (s >= effective_block_size()) {
		slot_ring.g2_automorphism(a);
	}
	return a;
}

void SlotRing::Rotation::apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const
{
	if (&in == &result) {
		throw std::invalid_argument("This function does not accept the same reference for in and result.");
	}
	if (s == 0) {
		result = in;
		return;
	}
	if (s == effective_block_size()) {
		result = in;
		eval.apply_galois_inplace(result, static_cast<uint32_t>(std::get<0>(galois_elements())), gk);
		return;
	}
	if (fmask_plain == nullptr) {
		fmask_plain = std::make_unique<seal::Plaintext>(gsl::span(fmask));
		bmask_plain = std::make_unique<seal::Plaintext>(gsl::span(bmask));
	}
	seal::Ciphertext backward;
	eval.multiply_plain(in, *fmask_plain, result);
	eval.multiply_plain(in, *bmask_plain, backward);
	eval.apply_galois_inplace(result, static_cast<uint32_t>(std::get<0>(galois_elements())), gk); log_galois();
	eval.apply_galois_inplace(backward, static_cast<uint32_t>(std::get<1>(galois_elements())), gk); log_galois();
	eval.add_inplace(result, backward);
}

const std::tuple<uint64_t, uint64_t> SlotRing::Rotation::get_forward_g1_g2_decomp() const
{
	if (block_size == slot_ring.n && !slot_ring.slot_group_cyclic) {
		if (s >= effective_block_size()) {
			return std::make_tuple(s % effective_block_size(), 1);
		}
		else {
			return std::make_tuple(s % effective_block_size(), 0);
		}
	}
	else {
		return std::make_tuple(s % effective_block_size(), 0);
	}
}

const std::tuple<uint64_t, uint64_t> SlotRing::Rotation::get_backward_g1_g2_decomp() const
{
	if (block_size == slot_ring.n && !slot_ring.slot_group_cyclic) {
		if (s >= effective_block_size()) {
			return std::make_tuple(slot_ring.g1_ord() + s % effective_block_size() - effective_block_size(), 0);
		}
		else {
			return std::make_tuple(slot_ring.g1_ord() + s % effective_block_size() - effective_block_size(), 1);
		}
	}
	else {
		return std::make_tuple(slot_ring.g1_ord() + s - effective_block_size(), 0);
	}
}

std::tuple<uint32_t, uint32_t> SlotRing::Rotation::galois_elements() const
{
	std::tuple<uint64_t, uint64_t> forward_g1_g2 = get_forward_g1_g2_decomp();
	std::tuple<uint64_t, uint64_t> backward_g1_g2 = get_backward_g1_g2_decomp();
	return std::make_tuple(
		static_cast<uint32_t>(seal::util::multiply_uint_mod(
			seal::util::exponentiate_uint_mod(slot_ring.g1, std::get<0>(forward_g1_g2), slot_ring.index_mod()),
			seal::util::exponentiate_uint_mod(slot_ring.g2, std::get<1>(forward_g1_g2), slot_ring.index_mod()),
			slot_ring.index_mod()
		)),
		static_cast<uint32_t>(seal::util::multiply_uint_mod(
			seal::util::exponentiate_uint_mod(slot_ring.g1, std::get<0>(backward_g1_g2), slot_ring.index_mod()),
			seal::util::exponentiate_uint_mod(slot_ring.g2, std::get<1>(backward_g1_g2), slot_ring.index_mod()),
			slot_ring.index_mod()
		))
	);
}

bool SlotRing::Rotation::is_identity() const
{
	return s == 0;
}

poly SlotRing::Frobenius::operator()(poly x) const
{
	if (iters == 0) {
		return x;
	}
	assert(x.size() <= slot_ring.N());
	x.resize(slot_ring.N());
	slot_ring.apply_frobenius(x, iters);
	return x;
}

void SlotRing::Frobenius::apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const
{
	if (iters == 0) {
		result = in;
		return;
	}
	eval.apply_galois(in, static_cast<uint32_t>(galois_element()), gk, result); log_galois();
}

const std::tuple<uint64_t, uint64_t> SlotRing::Frobenius::get_g1_g2_decomp() const
{
	return std::make_tuple(
		seal::util::multiply_uint_mod(std::get<0>(slot_ring.p_log), iters, slot_ring.index_mod()),
		seal::util::multiply_uint_mod(std::get<1>(slot_ring.p_log), iters, slot_ring.index_mod())
	);
}

uint32_t SlotRing::Frobenius::galois_element() const
{
	std::tuple<uint64_t, uint64_t> g1_g2_decomp = get_g1_g2_decomp();
	return static_cast<uint32_t>(seal::util::multiply_uint_mod(
		seal::util::exponentiate_uint_mod(slot_ring.g1, std::get<0>(g1_g2_decomp), slot_ring.index_mod()),
		seal::util::exponentiate_uint_mod(slot_ring.g2, std::get<1>(g1_g2_decomp), slot_ring.index_mod()),
		slot_ring.index_mod()
	));
}

bool SlotRing::Frobenius::is_identity() const
{
	return iters == 0;
}

poly SlotRing::RawAuto::operator()(poly x) const
{
	assert(x.size() <= slot_ring.N());
	x.resize(slot_ring.N());
	slot_ring.g1_automorphism(x, std::get<0>(g1_g2_decomp));
	if (std::get<1>(g1_g2_decomp) % 2 != 0) {
		slot_ring.g2_automorphism(x);
	}
	return x;
}

void SlotRing::RawAuto::apply_ciphertext(const seal::Ciphertext& in, const seal::Evaluator& eval, const seal::GaloisKeys& gk, seal::Ciphertext& result) const
{
	if (g1_g2_decomp == std::make_tuple(0, 0)) {
		result = in;
		return;
	}
	eval.apply_galois(in, galois_element(), gk, result); log_galois();
}

uint32_t SlotRing::RawAuto::galois_element() const
{
	uint64_t result = seal::util::exponentiate_uint_mod(slot_ring.g1, std::get<0>(g1_g2_decomp), slot_ring.index_mod());
	if (std::get<1>(g1_g2_decomp) % 2 != 0) {
		result = seal::util::multiply_uint_mod(result, slot_ring.g2, slot_ring.index_mod());
	}
	return static_cast<uint32_t>(result);
}

bool SlotRing::RawAuto::is_identity() const
{
	return galois_element() == 1;
}
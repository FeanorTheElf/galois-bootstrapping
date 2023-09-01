#include "polyarith.h"
#include "karatsuba.h"
#include <cstdlib>

struct uint64_mod_t {
	const seal::Modulus& modulus;
	uint64_t val;
};

const uint64_mod_t& operator-=(uint64_mod_t& lhs, const uint64_mod_t& rhs) {
	lhs.val = seal::util::sub_uint_mod(lhs.val, rhs.val, lhs.modulus);
	return lhs;
}

const uint64_mod_t& operator+=(uint64_mod_t& lhs, const uint64_mod_t& rhs) {
	lhs.val = seal::util::add_uint_mod(lhs.val, rhs.val, lhs.modulus);
	return lhs;
}

uint64_mod_t operator*(const uint64_mod_t& lhs, const uint64_mod_t& rhs) {
	uint64_mod_t result = { lhs.modulus, 0 };
	result.val = seal::util::multiply_uint_mod(lhs.val, rhs.val, lhs.modulus);
	return result;
}

void poly_normalize(poly& p) {
	int64_t i = p.size() - 1;
	while (i >= 0 && p[i] == 0) {
		--i;
	}
	p.resize(i + 1);
	assert(p.size() == 0 || p[p.size() - 1] != 0);
}

void poly_add(poly& lhs, const poly& rhs, const seal::Modulus& mod, uint64_t scale, size_t power_x_scale) {
	lhs.resize(std::max(lhs.size(), rhs.size() + power_x_scale));
	seal::util::MultiplyUIntModOperand scale_op;
	scale_op.set(scale, mod);
	for (size_t i = 0; i < rhs.size(); ++i) {
		lhs[i + power_x_scale] = seal::util::add_uint_mod(lhs[i + power_x_scale], seal::util::multiply_uint_mod(rhs[i], scale_op, mod), mod);
	}
}

void poly_scale(poly& p, uint64_t scale, const seal::Modulus& mod)
{
	seal::util::MultiplyUIntModOperand scale_op;
	scale_op.set(scale, mod);
	for (size_t i = 0; i < p.size(); ++i) {
		p[i] = seal::util::multiply_uint_mod(p[i], scale_op, mod);
	}
}

poly poly_mul_mod(const poly& lhs, const poly& rhs, const seal::Modulus& mod, const PolyModulus& pmod) {
	if (lhs.size() == 0 || rhs.size() == 0) {
		poly result;
		result.resize(pmod.n);
		return result;
	}
#ifdef CONTRACT_TEST
	poly check;
	poly_add_mul(check, lhs, rhs, mod);
#endif
	poly result;
	result.resize(lhs.size() + rhs.size());
	assert(lhs.size() <= pmod.n);
	assert(rhs.size() <= pmod.n);
	karatsuba::karatsuba<4>(
		&result[0], result.size(), 
		&lhs[0], lhs.size(), 
		&rhs[0], rhs.size(),
		[&mod](uint64_t a, uint64_t b) { return seal::util::add_uint_mod(a, b, mod); },
		[&mod](uint64_t a, uint64_t b) { return seal::util::sub_uint_mod(a, b, mod); },
		[&mod](uint64_t a, uint64_t b) { return seal::util::multiply_uint_mod(a, b, mod); }
	);
	poly_reduce_mod(result, mod, pmod);
#ifdef CONTRACT_TEST
	poly_reduce_mod(check, mod, pmod);
	assert(check.size() <= result.size());
	check.resize(result.size());
	assert(check == result);
#endif
	return result;
}

void poly_reduce_mod(poly& lhs, const seal::Modulus& mod, const PolyModulus& pmod)
{
	if (lhs.size() < pmod.n) {
		return;
	}
	for (int64_t i = lhs.size() - 1; i >= static_cast<int64_t>(pmod.n); --i) {
		poly_add(lhs, pmod.x_power_n, mod, lhs[i], i - pmod.n);
	}
	lhs.resize(pmod.n);
}

poly poly_exponentiate_mod(const poly& basis, size_t exp, const seal::Modulus& mod, const PolyModulus& pmod)
{
	bool changed = false;
	poly result = { 1 };
	for (int64_t i = std::numeric_limits<size_t>().digits - 1; i >= 0; --i) {
		if (((exp >> i) & 1) == 1) {
			result = poly_mul_mod(basis, poly_mul_mod(result, result, mod, pmod), mod, pmod);
			changed = true;
		}
		else if (changed) { // a performance optimization, no sense in squaring 1
			result = poly_mul_mod(result, result, mod, pmod);
		}
	}
	return result;
}

void poly_sub_mul(poly& dst, const poly& lhs, const poly& rhs, const seal::Modulus& mod)
{
	if (lhs.size() == 0 || rhs.size() == 0) {
		return;
	}
	dst.resize(std::max(dst.size(), lhs.size() + rhs.size() - 1));
	for (size_t i = 0; i < lhs.size(); ++i) {
		for (size_t j = 0; j < rhs.size(); ++j) {
			dst[i + j] = seal::util::sub_uint_mod(dst[i + j], seal::util::multiply_uint_mod(lhs[i], rhs[j], mod), mod);
		}
	}
}

void poly_add_mul(poly& dst, const poly& lhs, const poly& rhs, const seal::Modulus& mod)
{
	if (lhs.size() == 0 || rhs.size() == 0) {
		return;
	}
	dst.resize(std::max(dst.size(), lhs.size() + rhs.size() - 1));
	for (size_t i = 0; i < lhs.size(); ++i) {
		for (size_t j = 0; j < rhs.size(); ++j) {
			dst[i + j] = seal::util::add_uint_mod(dst[i + j], seal::util::multiply_uint_mod(lhs[i], rhs[j], mod), mod);
		}
	}
}

std::tuple<int64_t, int64_t> eea(uint64_t a, uint64_t b) {
	int64_t sa = 1;
	int64_t ta = 0;
	int64_t sb = 0;
	int64_t tb = 1;

	while (b != 0) {
		uint64_t quo = a / b;
		uint64_t rem = a % b;
		ta = ta - quo * tb;
		sa = sa - quo * sb;
		a = rem;
		std::swap(a, b);
		std::swap(sa, sb);
		std::swap(ta, tb);
	}
	return std::make_tuple(sa, ta);
}

uint64_t poly_eval(const poly& p, uint64_t x, const seal::Modulus& mod) {
	uint64_t result = p[p.size() - 1];
	for (int64_t i = p.size() - 2; i >= 0; --i) {
		result = seal::util::multiply_uint_mod(result, x, mod);
		result = seal::util::add_uint_mod(result, p[i], mod);
	}
	return result;
}

poly poly_div_rem(poly& lhs, const poly& rhs, const seal::Modulus& mod)
{
	assert(rhs.size() > 0);
	assert(rhs[rhs.size() - 1] == 1);
	if (rhs.size() > lhs.size()) {
		return {};
	}
#ifdef CONTRACT_TEST
	poly initial = lhs;
#endif
	poly result;
	result.resize(lhs.size() - rhs.size() + 1);
	for (int64_t i = lhs.size() - 1; i >= static_cast<int64_t>(rhs.size() - 1); --i) {
		result[i - rhs.size() + 1] = lhs[i];
		poly_sub(lhs, rhs, mod, lhs[i], i - rhs.size() + 1);
	}
	poly_normalize(result);
	poly_normalize(lhs);
#ifdef CONTRACT_TEST
	poly check = lhs;
	poly_add_mul(check, rhs, result, mod);
	assert(initial == check);
#endif
	return result;
}

std::tuple<poly, poly, poly> poly_eea(poly a, poly b, const seal::Modulus& mod)
{
	poly_normalize(a);
	poly_normalize(b);
	poly sa = { 1 };
	poly ta = {};
	poly sb = {};
	poly tb = { 1 };

#ifdef CONTRACT_TEST
	poly initial_a = a;
	poly initial_b = b;
#endif

	while (b.size() > 0) {
		if (b[b.size() - 1] != 1) {
			int64_t lc_inv = inv_mod(b[b.size() - 1], mod);
			poly_scale(b, lc_inv, mod);
			poly_scale(sb, lc_inv, mod);
			poly_scale(tb, lc_inv, mod);
		}
		poly quo = poly_div_rem(a, b, mod);
		poly_sub_mul(ta, quo, tb, mod);
		poly_sub_mul(sa, quo, sb, mod);
		std::swap(a, b);
		std::swap(sa, sb);
		std::swap(ta, tb);
	}
#ifdef CONTRACT_TEST
	poly check;
	poly_add_mul(check, sa, initial_a, mod);
	poly_add_mul(check, ta, initial_b, mod);
	poly_normalize(check);
	assert(check == a);
#endif
	return std::make_tuple(sa, ta, a);
}

uint64_t inv_mod(uint64_t x, const seal::Modulus& mod)
{
	int64_t result = std::get<0>(eea(x, mod.value()));
	if (result < 0) {
		result += mod.value();
	}
	assert(result >= 0);
	assert(static_cast<uint64_t>(result) < mod.value());
	assert(seal::util::multiply_uint_mod(x, result, mod) == 1);
	return static_cast<uint64_t>(result);
}

std::ostream& print_poly(std::ostream& os, const poly& f)
{
	if (f.size() == 0) {
		return os << "0";
	}
	for (size_t i = f.size() - 1; i >= 1; --i) {
		os << std::dec << f[i] << std::dec << " * x^" << i << "+";
	}
	return os << std::dec << f[0] << std::dec;
}

bool is_hex_digit(char c) {
	return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
}

poly parse_poly(std::string_view in)
{
	enum class ParseState {
		read_coeff, expect_power, read_power
	};

	poly result;
	std::string current;
	uint64_t coeff;
	uint64_t power;
	ParseState state = ParseState::read_coeff;
	auto finalize = [&]() {
		if (result.size() < power) {
			result.resize(power + 1);
		}
		result[power] = coeff;
	};
	for (char c : in) {
		if (c == ' ') {
			continue;
		}
		switch (state) {
		case ParseState::read_coeff:
			if (is_hex_digit(c)) {
				current.push_back(c);
			}
			else if (c == 'x') {
				coeff = strtoull(current.c_str(), nullptr, 16);
				state = ParseState::expect_power;
			}
			else if (c == '+') {
				coeff = strtoull(current.c_str(), nullptr, 16);
				power = 0;
				finalize();
				state = ParseState::read_coeff;
				current.clear();
			}
			else {
				throw "invalid format";
			}
			break;
		case ParseState::expect_power:
			if (c == '^') {
				state = ParseState::read_power;
				current.clear();
			}
			else if (c == '+') {
				power = 1;
				finalize();
				state = ParseState::read_coeff;
				current.clear();
			}
			else {
				throw "invalid format";
			}
			break;
		case ParseState::read_power:
			if (is_hex_digit(c)) {
				current.push_back(c);
			}
			else if (c == '+') {
				power = strtoull(current.c_str(), nullptr, 10);
				finalize();
				state = ParseState::read_coeff;
				current.clear();
			}
			else {
				throw "invalid format";
			}
			break;
		}
	}
	switch (state) {
	case ParseState::read_coeff:
		coeff = strtoull(current.c_str(), nullptr, 16);
		power = 0;
		finalize();
		break;
	case ParseState::expect_power:
		power = 1;
		finalize();
		break;
	case ParseState::read_power:
		power = strtoull(current.c_str(), nullptr, 10);
		finalize();
		state = ParseState::read_coeff;
		current.clear();
		break;
	}
	return result;
}

bool is_zero(const poly& x)
{
	for (uint64_t a : x) {
		if (a != 0) {
			return false;
		}
	}
	return true;
}

bool is_irreducible(const poly& x, const seal::Modulus& prime_modulus)
{
	assert(x.size() > 0);
	if (x.size() <= 2) {
		return true;
	}
	assert(x[x.size() - 1] == 1);
	PolyModulus mod_x;
	mod_x.n = x.size() - 1;
	poly mod_poly = x;
	mod_poly[x.size() - 1] = 0;
	poly_normalize(mod_poly);
	poly_scale(mod_poly, seal::util::negate_uint_mod(1, prime_modulus), prime_modulus);
	mod_x.x_power_n = std::move(mod_poly);

	const size_t p = prime_modulus.value();
	const size_t d = x.size() - 1;
	poly tmp = poly_exponentiate_mod({ 0, 1 }, ipow(p, d), prime_modulus, mod_x);
	poly_sub(tmp, { 0, 1 }, prime_modulus);
	poly_reduce_mod(tmp, prime_modulus, mod_x);
	poly_normalize(tmp);
	if (tmp.size() > 0) {
		return false;
	}

	std::vector<uint64_t> primes;
	for (uint64_t n = 2; n <= d; ++n) {
		if (std::any_of(primes.cbegin(), primes.cend(), [n](uint64_t p) { return n % p == 0; })) {
			continue;
		}
		primes.push_back(n);

		poly tmp = poly_exponentiate_mod({ 0, 1 }, ipow(p, d / n), prime_modulus, mod_x);
		poly_sub(tmp, { 0, 1 }, prime_modulus);
		poly_reduce_mod(tmp, prime_modulus, mod_x);
		poly_normalize(tmp);
		if (tmp.size() == 0) {
			return false;
		}
	}

	return true;
}

void test_is_irreducible() {
	{
		seal::Modulus mod(2);
		assert(is_irreducible({ 0, 1 }, mod));
		assert(is_irreducible({ 1, 1, 1 }, mod));
		assert(!is_irreducible({ 1, 0, 1 }, mod));
		assert(!is_irreducible({ 0, 0, 1 }, mod));
		assert(is_irreducible({ 1, 0, 1, 1 }, mod));
	}
	{
		seal::Modulus mod(127);
		assert(is_irreducible({ 115, 22, 36, 1 }, mod));
		assert(!is_irreducible({ 1, 0, 0, 1 }, mod));
		assert(!is_irreducible({ 1, 1, 0, 1 }, mod));
	}
	std::cout << "test_is_irreducible(): success" << std::endl;
}
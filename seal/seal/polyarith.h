#pragma once
#include <vector>
#include <assert.h>
#include <ostream>
#include "seal/util/uintarithsmallmod.h"

/**
 * Somewhat fast arithmetic on polynomials, used to compute
 * in the plaintext space of the scheme.
*/

inline uint64_t ipow(uint64_t base, uint64_t exp) {
	uint64_t result = 1;
	for (int64_t i = std::numeric_limits<size_t>().digits - 1; i >= 0; --i) {
		if (((exp >> i) & 1) == 1) {
			result = result * result * base;
		}
		else {
			result = result * result;
		}
	}
	return result;
}

typedef std::vector<uint64_t> poly;

struct PolyModulus {
	poly x_power_n;
	size_t n;
};

void poly_normalize(poly& p);
void poly_add(poly& lhs, const poly& rhs, const seal::Modulus& mod, uint64_t scale = 1, size_t power_x_scale = 0);
poly poly_exponentiate_mod(const poly& basis, size_t exp, const seal::Modulus& mod, const PolyModulus& pmod);
void poly_scale(poly& p, uint64_t scale, const seal::Modulus& mod);
poly poly_mul_mod(const poly& lhs, const poly& rhs, const seal::Modulus& mod, const PolyModulus& pmod);
void poly_reduce_mod(poly& lhs, const seal::Modulus& mod, const PolyModulus& pmod);
void poly_sub_mul(poly& dst, const poly& lhs, const poly& rhs, const seal::Modulus& mod);
void poly_add_mul(poly& dst, const poly& lhs, const poly& rhs, const seal::Modulus& mod);
std::tuple<int64_t, int64_t> eea(uint64_t a, uint64_t b);
uint64_t poly_eval(const poly& p, uint64_t x, const seal::Modulus& mod);
poly poly_div_rem(poly& lhs, const poly& rhs, const seal::Modulus& mod);
std::tuple<poly, poly, poly> poly_eea(poly lhs, poly rhs, const seal::Modulus& mod);
uint64_t inv_mod(uint64_t x, const seal::Modulus& mod);
std::ostream& print_poly(std::ostream& os, const poly& f);
poly parse_poly(std::string_view in);
bool is_zero(const poly& x);

// prime_modulus must be a prime!
bool is_irreducible(const poly& x, const seal::Modulus& prime_modulus);

inline void poly_sub(poly& lhs, const poly& rhs, const seal::Modulus& mod, uint64_t scale = 1, size_t power_x_scale = 0) {
	poly_add(lhs, rhs, mod, seal::util::negate_uint_mod(scale, mod), power_x_scale);
}

void test_is_irreducible();
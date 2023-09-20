#pragma once
#include <vector>
#include <assert.h>

/**
 * For plaintext arithmetic, we use karatsuba multiplication
*/
namespace karatsuba {

	inline size_t leading_zeros(uint64_t mask) {
		unsigned long result = 64;
		while (mask != 0) {
			mask = mask >> 1;
			result -= 1;
		}
		return result;
	}

	template<typename T, typename F>
	void elementwise_assign(T* dst, const T* src, size_t n, const F& f);

	template<bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	void naive_assign_mul(T* dst, const T* lhs, const T* rhs, size_t in_size_log2, const Add& add, const Sub& sub, const Mul& mul);

	template<unsigned in_size_log2, unsigned threshold_log2, bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	void karatsuba_assign_mul_impl(T* dst, const T* lhs, const T* rhs, T* memory, const Add& add, const Sub& sub, const Mul& mul);

	template<unsigned threshold_log2, bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	void dispatch_karatsuba_assign_mul_impl(unsigned in_size_log2, T* dst, const T* lhs, const T* rhs, T* memory, const Add& add, const Sub& sub, const Mul& mul);

	template<unsigned threshold_log2, typename T, typename Add, typename Sub, typename Mul>
	void karatsuba(T* dst, size_t dst_len, const T* lhs, size_t lhs_len, const T* rhs, size_t rhs_len, const Add& add, const Sub& sub, const Mul& mul);

	template<typename T>
	size_t karatsuba_mem_element_count(size_t in_size_log2, size_t threshold_log2);

	template<typename T, typename F>
	inline void elementwise_assign(T* dst, const T* src, size_t n, const F& f)
	{
		for (size_t i = 0; i < n; ++i) {
			dst[i] = f(dst[i], src[i]);
		}
	}

	template<bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	inline void naive_assign_mul(T* dst, const T* lhs, const T* rhs, size_t in_size_log2, const Add& add, const Sub& sub, const Mul& mul)
	{
		const size_t n = (size_t)1 << in_size_log2;
		for (size_t i = 0; i < 2 * n; ++i) {
			T acc = 0;
			const size_t start = static_cast<size_t>(std::max<int64_t>(static_cast<int64_t>(i) - static_cast<int64_t>(n) + 1, 0));
			for (size_t j = start; j < std::min(n, i + 1); ++j) {
				acc = add(acc, mul(lhs[i - j], rhs[j]));
			}
			if constexpr (add_assign) {
				dst[i] = add(dst[i], acc);
			}
			else {
				dst[i] = acc;
			}
		}
	}

	template<unsigned in_size_log2, unsigned threshold_log2, bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	void karatsuba_assign_mul_impl(T* dst, const T* lhs, const T* rhs, T* memory, const Add& add, const Sub& sub, const Mul& mul)
	{
		if constexpr (in_size_log2 <= threshold_log2) {
			naive_assign_mul<add_assign>(dst, lhs, rhs, in_size_log2, add, sub, mul);
		}
		else {
			// n is half the length of lhs resp. rhs
			constexpr size_t n = (size_t)1 << (in_size_log2 - 1);
			// write lhs = a0 + X^n a1 and same for rhs
			T* lower = memory;
			T* rest = memory + (2 * n);

			// set lower to a0 * b0
			karatsuba_assign_mul_impl<in_size_log2 - 1, threshold_log2, false>(lower, lhs, rhs, rest, add, sub, mul);
			if constexpr (add_assign) {
				elementwise_assign(dst, lower, 2 * n, add);
			}
			else {
				std::memcpy(dst, lower, 2 * n * sizeof(T));
				std::memset(dst + 2 * n, 0, 2 * n * sizeof(T));
			}
			elementwise_assign(dst + n, lower, 2 * n, sub);

			// set upper to a1 * b1
			T* upper = lower; // reuse memory
			karatsuba_assign_mul_impl<in_size_log2 - 1, threshold_log2, false>(upper, lhs + n, rhs + n, rest, add, sub, mul);
			elementwise_assign(dst + 2 * n, upper, 2 * n, add);
			elementwise_assign(dst + n, upper, 2 * n, sub);

			// and now we come to the mid term X^n (a0 * b1 + a1 * b0)
			T* lhs_combined = upper; // reuse memory
			T* rhs_combined = upper + n;
			std::memcpy(lhs_combined, lhs, n * sizeof(T));
			std::memcpy(rhs_combined, rhs, n * sizeof(T));
			elementwise_assign(lhs_combined, lhs + n, n, add);
			elementwise_assign(rhs_combined, rhs + n, n, add);
			karatsuba_assign_mul_impl<in_size_log2 - 1, threshold_log2, true>(dst + n, lhs_combined, rhs_combined, rest, add, sub, mul);
		}
	}

	template<unsigned threshold_log2, bool add_assign, typename T, typename Add, typename Sub, typename Mul>
	void dispatch_karatsuba_assign_mul_impl(unsigned in_size_log2, T* dst, const T* lhs, const T* rhs, T* memory, const Add& add, const Sub& sub, const Mul& mul)
	{
		if (in_size_log2 > 16) {
			throw std::invalid_argument("Size must be between 1 and 2^16");
		}
		switch (in_size_log2) {
		case 0:
			karatsuba_assign_mul_impl<0, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 1:
			karatsuba_assign_mul_impl<1, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 2:
			karatsuba_assign_mul_impl<2, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 3:
			karatsuba_assign_mul_impl<3, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 4:
			karatsuba_assign_mul_impl<4, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 5:
			karatsuba_assign_mul_impl<5, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 6:
			karatsuba_assign_mul_impl<6, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 7:
			karatsuba_assign_mul_impl<7, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 8:
			karatsuba_assign_mul_impl<8, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 9:
			karatsuba_assign_mul_impl<9, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 10:
			karatsuba_assign_mul_impl<10, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 11:
			karatsuba_assign_mul_impl<11, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 12:
			karatsuba_assign_mul_impl<12, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 13:
			karatsuba_assign_mul_impl<13, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 14:
			karatsuba_assign_mul_impl<14, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 15:
			karatsuba_assign_mul_impl<15, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		case 16:
			karatsuba_assign_mul_impl<16, threshold_log2, add_assign, T, Add, Sub, Mul>(dst, lhs, rhs, memory, add, sub, mul);
			break;
		}
	}

	template<unsigned threshold_log2, typename T, typename Add, typename Sub, typename Mul>
	void karatsuba(T* dst, size_t dst_len, const T* lhs, size_t lhs_len, const T* rhs, size_t rhs_len, const Add& add, const Sub& sub, const Mul& mul)
	{
		if (lhs_len == 0 || rhs_len == 0) {
			return;
		}
		assert(dst_len >= lhs_len + rhs_len);
		int block_size_log2 = (63 - static_cast<int>(std::max(leading_zeros(lhs_len), leading_zeros(rhs_len))));
		size_t n = (size_t)1 << block_size_log2;
		assert(lhs_len >= n);
		assert(rhs_len >= n);
		assert((lhs_len < 2 * n) || (rhs_len < 2 * n));
		const size_t mem_element_count = karatsuba_mem_element_count<T>(static_cast<unsigned int>(block_size_log2), threshold_log2);
		std::unique_ptr<T[] > memory(new T[mem_element_count]);

		for (size_t i = 0; i + n <= lhs_len; i += n) {
			for (size_t j = 0; j + n <= rhs_len; j += n) {
				dispatch_karatsuba_assign_mul_impl<threshold_log2, true>(static_cast<unsigned int>(block_size_log2), dst + i + j, lhs + i, rhs + j, memory.get(), add, sub, mul);
			}
		}

		// now make block sizes smaller and compute the remaining elements
		size_t lhs_rem = (lhs_len / n) * n;
		size_t rhs_rem = (rhs_len / n) * n;
		block_size_log2 -= 1;
		n = n / 2;
		while (block_size_log2 >= 0) {
			if (lhs_len >= lhs_rem + n) {
				for (size_t j = 0; j + n <= rhs_rem; j += n) {
					dispatch_karatsuba_assign_mul_impl<threshold_log2, true>(static_cast<unsigned int>(block_size_log2), dst + lhs_rem + j, lhs + lhs_rem, rhs + j, memory.get(), add, sub, mul);
				}
				lhs_rem += n;
			}
			if (rhs_len >= rhs_rem + n) {
				for (size_t i = 0; i + n <= lhs_len; i += n) {
					dispatch_karatsuba_assign_mul_impl<threshold_log2, true>(static_cast<unsigned int>(block_size_log2), dst + rhs_rem + i, lhs + i, rhs + rhs_rem, memory.get(), add, sub, mul);
				}
				rhs_rem += n;
			}
			n = n / 2;
			block_size_log2 -= 1;
		}
	}

	template<typename T>
	inline size_t karatsuba_mem_element_count(size_t in_size_log2, size_t threshold_log2)
	{
		if (in_size_log2 <= threshold_log2) {
			return 0;
		}
		const size_t elements = ((size_t)1 << (in_size_log2 + 1)) - ((size_t)1 << (threshold_log2 + 1));
		return elements;
	}

	inline void test_karatsuba() {
		auto add = [](int a, int b) { return a + b; };
		auto sub = [](int a, int b) { return a - b; };
		auto mul = [](int a, int b) { return a * b; };
		{
			std::vector<int> lhs = { 0, 1, 1 };
			std::vector<int> rhs = { 1, 0, 0, 1 };
			std::vector<int> dst;
			dst.resize(7);
			karatsuba<1>(&dst[0], dst.size(), &lhs[0], lhs.size(), &rhs[0], rhs.size(), add, sub, mul);
			std::vector<int> expected = { 0, 1, 1, 0, 1, 1, 0 };
			assert(dst == expected);
		}
		{
			std::vector<int> lhs = { 1, 1, 0, 1, 0, 2, 2, 2 };
			std::vector<int> rhs = { 2, 0, 0, -1, 1, 1, 1 };
			std::vector<int> dst;
			dst.resize(15);
			karatsuba<1>(&dst[0], dst.size(), & lhs[0], lhs.size(), &rhs[0], rhs.size(), add, sub, mul);
			std::vector<int> expected = { 2, 2, 0, 1, 0, 6, 5, 6, -1, 1, 2, 6, 4, 2, 0 };
			assert(dst == expected);
		}
		{
			std::vector<int> lhs = { 0, 1, 0, 2, 2 };
			std::vector<int> rhs = { 2, 0, 0, -1, 1, 1, 1 };
			std::vector<int> dst;
			dst.resize(12);
			karatsuba<2>(&dst[0], dst.size(), &lhs[0], lhs.size(), &rhs[0], rhs.size(), add, sub, mul);
			std::vector<int> expected = { 0, 2, 0, 4, 3, 1, -1, 1, 4, 4, 2, 0 };
			assert(dst == expected);
		}
	}
}
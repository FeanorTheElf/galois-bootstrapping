#include <iostream>
#include <vector>
#include <assert.h>
#include <fstream>
#include <chrono>
#include "slots.h"
#include "transform.h"
#include "seal/seal.h"
#include "karatsuba.h"
#include "bootstrapping.h"

using namespace seal;
using namespace std;

int main()
{
	// set up context
	EncryptionParameters parms(scheme_type::bfv);
	std::shared_ptr<SlotRing> slot_ring = p_257_test_parameters();
	SlotRing basic_slot_ring = slot_ring->change_exponent(1);
	parms.set_poly_modulus_degree(slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring->N()));
	parms.set_plain_modulus(slot_ring->prime());
	SEALContext context(parms);
	std::cout << "Q bit count is " << context.key_context_data()->total_coeff_modulus_bit_count() << std::endl;

	// set up keys
	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	RelinKeys rk;
	keygen.create_relin_keys(rk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);

	// create an example plaintext
	poly data;
	poly_add(data, basic_slot_ring.from_slot_value({ 1 }, 0), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 8 }, 1), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 63 }, 2), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 17 }, 3), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 12 }, 4), basic_slot_ring.R().scalar_mod);
	Plaintext x_plain{ gsl::span<const uint64_t>(data) };

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);
	std::cout << "encrypted message" << std::endl;
	std::cout << "noise budget is " << decryptor.invariant_noise_budget(x_enc) << " bits" << std::endl;

	for (size_t i = 0; i < 24; ++i) {
		evaluator.square_inplace(x_enc);
		evaluator.relinearize_inplace(x_enc, rk);
	}

	std::cout << "performed computations" << std::endl;
	std::cout << "noise budget is " << decryptor.invariant_noise_budget(x_enc) << " bits" << std::endl;

	std::cout << "setup context" << std::endl;

	Bootstrapper bootstrapper(context, slot_ring, p_257_test_parameters_digit_extractor(*slot_ring));
	bootstrapper.initialize();

	std::cout << "initialized bootstrapper" << std::endl;

	BootstrappingKey bk;
	bootstrapper.create_bootstrapping_key(sk, bk);

	std::cout << "created bootstrapping key" << std::endl;

	auto start = std::chrono::steady_clock::now();

	Ciphertext in_coeffs;
	bootstrapper.slots_to_coeffs(x_enc, bk, in_coeffs, MemoryManager::GetPool() DEBUG_PASS(sk));

	auto end = std::chrono::steady_clock::now();
	std::cout << "moved to coefficients in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;
	kswitch_counter = 0;
	start = end;

	Ciphertext noisy_dec;
	bootstrapper.homomorphic_noisy_decrypt(in_coeffs, bk, noisy_dec, MemoryManager::GetPool() DEBUG_PASS(sk));

	end = std::chrono::steady_clock::now();
	std::cout << "homomorphically decrypted in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;
	start = end;
	kswitch_counter = 0;

	Ciphertext in_slots;
	bootstrapper.coeffs_to_slots(noisy_dec, bk, in_slots, MemoryManager::GetPool() DEBUG_PASS(sk));

	end = std::chrono::steady_clock::now();
	std::cout << "moved to slots in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;
	start = end;
	kswitch_counter = 0;

	Ciphertext digit_extracted;
	bootstrapper.slotwise_digit_extract(in_slots, bk, digit_extracted, MemoryManager::GetPool() DEBUG_PASS(sk));

	end = std::chrono::steady_clock::now();
	std::cout << "bootstrapped in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;

	Plaintext result;
	decryptor.decrypt(digit_extracted, result);

	std::cout << "noise budget is " << decryptor.invariant_noise_budget(digit_extracted) << " bits" << std::endl;

	poly result_poly(result.data(), result.data() + result.coeff_count());
	result_poly.resize(basic_slot_ring.N());

	for (size_t i = 0; i < basic_slot_ring.slot_group_len(); ++i) {
		std::cout << basic_slot_ring.extract_slot_value(result_poly, i)[0] << std::endl;
	}

	return 0;
	test_compile_slot_basis();
	test_rotate_noncyclic();
	test_is_irreducible();
	test_g_automorphisms();
	test_block_rotate();
	test_apply_ciphertext();
	//test_apply_ciphertext_subring();
	//test_compile_pair_transformation();
	test_slotwise_digit_extract();
	test_first_coeffs_to_scalar_slots();
	test_coeffs_to_slots();
	test_galois_poly_evaluator();
	test_homomorphic_noisy_decrypt();
	return 0;
}
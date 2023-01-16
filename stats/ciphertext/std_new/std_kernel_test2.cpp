#include "openfhe.h"
#include <random>
#include <chrono>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
        int n = 15; //The number of data
        double a = 1./(double)n;

	CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

	for (int d = 1; d < 2; d++) {
		cout << "Test dimension " << d << endl;
		int nd = max(n, d);
		//Rotation keys
		vector<int> ro;
		for (int i = 0; i < nd; i++) {
			ro.push_back(i+1);
			ro.push_back(-(i+1));

		}
        	cc -> EvalRotateKeyGen(keys.secretKey, ro);

        
		//Inputs & Encoding & Encryption
                std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));

     		for (int z = 0; z < n; z++) {
                     std::vector<double> x;
                     int count = 0;
                     while (count < d) {
			     double data = dist(gen);
			     x.push_back(data);
			     count += 1;                                
		     }
		     dataset_original[z] = x;
		}

	
		std::vector<Plaintext> plaintext;

	        for (int i = 0; i < n; i++) {
        	        Plaintext ptxt = cc -> MakeCKKSPackedPlaintext(dataset_original[i]);
                	plaintext.push_back(ptxt);
        	}

	        std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dataset;

        	for (int i = 0; i < n; i++) {
                	auto c = cc -> Encrypt(keys.publicKey, plaintext[i]);
	                dataset.push_back(c);
        	}

		//Generating zero vector; its dimension is d
		vector<double> zero_d_original(d, 0.0);
		Plaintext ptxt_zero_d = cc -> MakeCKKSPackedPlaintext(zero_d_original);
		auto ct_zero_d = cc -> Encrypt(keys.publicKey, ptxt_zero_d);

		//Generating zero vector; its dimension is n
		vector<double> zero_n_original(n, 0.0);
		Plaintext ptxt_zero_n = cc -> MakeCKKSPackedPlaintext(zero_n_original);
		auto ct_zero_n = cc -> Encrypt(keys.publicKey, ptxt_zero_n);

		//Generating one vectors; its dimension is d
		vector<Plaintext> one_d_plaintext;
		for (int i = 0; i < d; i++) {
			vector<double> one_d = zero_d_original;
			one_d[i] = 1.0;
			Plaintext ptxt_one_d = cc -> MakeCKKSPackedPlaintext(one_d);
			one_d_plaintext.push_back(ptxt_one_d);
		}

                //Generating one vectors; its dimension is n
                vector<Plaintext> one_n_plaintext;
                for (int i = 0; i < n; i++) {
                        vector<double> one_n = zero_n_original;
                        one_n[i] = 1.0;
                        Plaintext ptxt_one_n = cc -> MakeCKKSPackedPlaintext(one_n);
                        one_n_plaintext.push_back(ptxt_one_n);
                }


		//Making kernel; Scalar multiplication kernel
		vector<vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel_zero(n ,vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > (n, ct_zero_n));
		vector<vector<vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > > kernel_vector(d, kernel_zero); //The vector of kernels; the number of elements is d

		for (int i = 0; i < n; i++) {//Select data_i
			auto data_i = dataset[i];
			for (int j = 0; j < n; j++) {//Select data_j
				auto data_j = dataset[j];
				auto increment = cc -> EvalMult(data_i, data_j); //For each coordinate k, its element is the increment of i,j-element of kernel_k; For example, increment[k] = kernel_k(i, j)

				for (int k = 0; k < d; k++) {//Select dimension
					auto kernel_k = kernel_vector[k];
					auto kernel_k_i_j = kernel_k[i][j]; //i,j-element of kernel_k
					auto increment_k = cc -> EvalMult(increment, one_d_plaintext[k]); //Delete except k-th element

					if (k == 0) {
						kernel_k_i_j = cc -> EvalAdd(kernel_k_i_j, increment_k);
					} else {
						increment_k = cc -> EvalRotate(increment_k, k);
						kernel_k_i_j = cc -> EvalAdd(kernel_k_i_j, increment_k);
					}
					kernel_k[i][j] = kernel_k_i_j; //Update i-th row of kernel_k
					kernel_vector[k] = kernel_k; //Update kernel_k
				}
			}
		}
		
		//Calculate std
	
		system_clock::time_point start_time = system_clock::now();

		//First term
		auto first_term = ct_zero_n; //Initial value
		for (int i = 0; i < d; i++) {//Select dimension
			auto kernel_i = kernel_vector[i];
			auto sum_diag_kernel = ct_zero_n; //Initial value
			for (int j = 0; j < n; j++) {//Select row & column
				sum_diag_kernel = cc -> EvalAdd(sum_diag_kernel, kernel_i[j][j]);	
			} //The first element of sum_diag_kernel is the summation of diagonal elements
			sum_diag_kernel = cc -> EvalMult(sum_diag_kernel, one_n_plaintext[0]); //Eliminate except the first element
			if (i == 0) {
				;
			} else {
				sum_diag_kernel = cc -> EvalRotate(sum_diag_kernel, -i);
			} //The i-th element of sum_diag_kernel is the increment
			first_term = cc -> EvalAdd(first_term, sum_diag_kernel); //Update i-coordinate of first_term
		}
		first_term = cc -> EvalMult(first_term, a); //Multiply 1/n

		//Second term
		auto second_term = ct_zero_n;
		for (int i = 0; i < d; i++) {//Select dimension
			auto kernel_i = kernel_vector[i];
			auto sum_all_kernel = ct_zero_n; //Initial value
			for (int j = 0; j < n; j++) {//Select row
				for (int k = 0; k < n; k++) {//Select column
					sum_all_kernel = cc -> EvalAdd(sum_all_kernel, kernel_i[j][k]);
				}
			} //The first element of sum_all_kernel is summation of all elements	
			sum_all_kernel = cc -> EvalMult(sum_all_kernel, one_n_plaintext[0]); //Eliminate except the first element
			sum_all_kernel = cc -> EvalMult(sum_all_kernel, a*a); //Multiply 1/n^2
			if (i == 0) {
				;
			} else {
				sum_all_kernel = cc -> EvalRotate(sum_all_kernel, -i);
			} //The i-th element of sum_all_kernel is the increment
			second_term = cc -> EvalAdd(second_term, sum_all_kernel); //Update i-coordinate of second_term
		}
		second_term = cc -> EvalMult(second_term, a*a); //Multiply 1/n^2

		//Variance & std
		auto variance = cc -> EvalSub(first_term, second_term);
		auto std = sqrt(variance, cc);

		system_clock::time_point end_time = system_clock::now();
		microseconds micro_eval = duration_cast<microseconds>(end_time - start_time);

		//Decryption & Decoding
		Plaintext result;
		cc -> Decrypt(keys.secretKey, std, &result);
		result -> SetLength(d);
		
		cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;
		
	}

}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
        ctxs = cc -> EvalMult(ctxs, 0.0001);
        auto a = ctxs;
        auto b = cc->EvalSub(ctxs, 1.0);

        for (int i = 0; i<6; i++) {
                auto b_half = cc -> EvalMult(b, -0.5); // -b/2
                a = cc->EvalMult(a, cc->EvalAdd(b_half, 1.0)); // a = a*(-b/2 + 1)
                auto tmp = cc -> EvalSub(b, 3.0); // b-3
                auto b_quater = cc -> EvalMult(tmp, 0.25); // (b-3)/4
                auto b_square = cc -> EvalMult(b, b); // b**2
                b = cc -> EvalMult(b_square, b_quater); // b = (b**2)*((b-3)/4)
        }

        float constant = sqrt(10000);

        a = cc -> EvalMult(a, constant);

        return a;
}


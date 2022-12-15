#include "openfhe.h"
#include <time.h>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 50;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
	int n = 2; //The number of data

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);
        std::cout << "CKKS scheme is using ring dimension " << cc -> GetRingDimension() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using scaling mod size " << parameters.GetScalingModSize() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using security level " << parameters.GetSecurityLevel() <<std::endl << std::endl;

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

	std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

        //test 10 case
        for (int d = 14; d < 30; d++) { //dimension
                for (int testcase = 0; testcase < 1; testcase++) { //Test 10 case for each dimension
                        cout << "Test dimension " << d << " and case " << testcase + 1 << endl;

                        //Generate rotation key
                        int nd = max(n, d);
                        vector<int> ro;
                                for (int i = 1; i < nd + 1; i++) {
                                        ro.push_back(i);
                                        ro.push_back(-i);
                                }
                        cc -> EvalRotateKeyGen(keys.secretKey, ro);

                        //Generating dataset
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

			//Generate one vectors; {1, 0, 0, ... }, {0, 1, 0, ... }
                        std::vector<Plaintext> ptxt_one;
                        std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > vector_one;

                        for (int i = 0; i < d; i++) {
                                std::vector<double> one(d, 0.0);
                                one[i] = 1.0;
                                Plaintext one_tmp = cc -> MakeCKKSPackedPlaintext(one);
                                ptxt_one.push_back(one_tmp);
                                auto ct_one = cc -> Encrypt(keys.publicKey, one_tmp);
                                vector_one.push_back(ct_one);
                        }

			//Generate n-dimensional zero vector
                        vector<double> zero_n(n, 0.0);
                        Plaintext ptxt_zero_n = cc -> MakeCKKSPackedPlaintext(zero_n);
                        auto ct_zero_n = cc -> Encrypt(keys.publicKey, ptxt_zero_n);

		        //Encoding & Encryption

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

			//Make kernel
	
			system_clock::time_point start_time_kernel = system_clock::now(); //Measure the time of making kernel

			std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel(n, ct_zero_n); //Kernel matrix

                        //Diagonal elements
                        for (int i = 0; i < n; i++) {
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i]; //ith row of kernel
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square = cc -> EvalMult(dataset[i], dataset[i]); //(d1*d1, d2*d2, ...)
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = data_square;
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotate = data_square;

                                for (int j = 1; j < d; j++) {
                                        data_square_rotate = cc -> EvalRotate(data_square_rotate, 1);
                                        inner_product = cc -> EvalAdd(inner_product, data_square_rotate);
                                } //The first element of inner_product is the inner product between data

                                inner_product = cc -> EvalMult(inner_product, vector_one[0]); //Eliminate dummy elements

                                if (i == 0) {
                                        ith_row = cc -> EvalAdd(ith_row, inner_product);
                                } else {
                                        inner_product = cc -> EvalRotate(inner_product, -i);
                                        ith_row = cc -> EvalAdd(ith_row, inner_product);
                                } //Update ith row

                                kernel[i] = ith_row;

                        }

                        //Non-diagonal elements
                        for (int i = 0; i < n-1; i++) { //select row
                                for (int j = i+1; j < n; j++) { //select column
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i]; //ith row of kernel
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[j]);
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;

                                        for (int k = 1; k < d; k++) {
                                                d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
                                                inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
                                        } //The first element of inner_product is the inner product between data_i and data_j

                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element = cc -> EvalMult(inner_product, ptxt_one[0]); //Eliminate dummy elements
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element_rotate = cc -> EvalRotate(i_j_element, -j);
                                        ith_row = cc -> EvalAdd(ith_row, i_j_element_rotate); //Updata ith row
                                        kernel[i] = ith_row;

                                        //By symmetry, we can also update jth row
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > jth_row = kernel[j];
                                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > j_i_element_rotate;
                                        if (i == 0) {
                                                j_i_element_rotate = i_j_element;
                                        } else {
                                                j_i_element_rotate = cc -> EvalRotate(i_j_element, -i);
                                        }
                                        jth_row = cc -> EvalAdd(jth_row, j_i_element_rotate); //Update jth row
                                        kernel[j] = jth_row;
                                }
                        }

			system_clock::time_point end_time_kernel = system_clock::now();
			microseconds micro_kernel = duration_cast<microseconds>(end_time_kernel - start_time_kernel);
			std::cout << "Time of making kernel: " << micro_kernel.count()/1000000. << " s" << std::endl;
			
			//Calculate fraction
	
			system_clock::time_point start_time_eval = system_clock::now();

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > kernel0 = kernel[0];
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > kernel1 = kernel[1];

			kernel1 = cc -> EvalRotate(kernel1, 1); //The first element is kernel[1][1]

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > fraction_square = cc -> EvalMult(kernel0, kernel1); //The first element is kernel[0][0]*kernel[1][1]

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > fraction = sqrt(fraction_square, cc); //The first element is fraction

			//Inverse of fraction
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > fraction_inverse = inverse(fraction, cc); //The first element is the inverse of fraction

			//Calculate denominator
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > denominator = cc -> EvalRotate(kernel0, 1); //Denominator; The first element is kernel[0][1]

			//Calculate similarity
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > similarity = cc -> EvalMult(fraction_inverse, denominator);

			system_clock::time_point end_time_eval = system_clock::now();
	                microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);

        		cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, similarity, &result1);
			result1 -> SetLength(1);
		}
	}
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
        ctxs = cc -> EvalMult(ctxs, 0.001);
	auto a = ctxs;
        auto b = cc->EvalSub(ctxs, 1.0);

        for (int i = 0; i<10; i++) {
                auto b_half = cc -> EvalMult(b, -0.5); // -b/2
                a = cc->EvalMult(a, cc->EvalAdd(b_half, 1.0)); // a = a*(-b/2 + 1)
                auto tmp = cc -> EvalSub(b, 3.0); // b-3
                auto b_quater = cc -> EvalMult(tmp, 0.25); // (b-3)/4
                auto b_square = cc -> EvalMult(b, b); // b**2
                b = cc -> EvalMult(b_square, b_quater); // b = (b**2)*((b-3)/4)
        }

	float constant = sqrt(1000);
	a = cc -> EvalMult(a, constant);
        return a;
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc) {
	ct = cc -> EvalMult(ct, 0.001);
	auto ct1 = cc->EvalSub(1.0,ct);
	auto ct2 = cc->EvalAdd(ct1,1.0);

	for(int i=0; i<10; i++){
        	ct1 = cc->EvalMult(ct1,ct1);
	        ct2 = cc->EvalMult(ct2,cc->EvalAdd(ct1,1.0));
    	}
	
	ct2 = cc -> EvalMult(ct2, 0.001);

	return ct2;
}


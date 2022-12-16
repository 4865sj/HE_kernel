#include "openfhe.h"
#include <time.h>
#include <random>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace std;
using namespace lbcrypto;
using namespace chrono;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
	int d = 1; //Dimension

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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

        ofstream fout;
        fout.open("variance_kernel.txt"); //Save the results of time

        //test 10 case
        for (int n = 1; n < 16; n++) { //The number of data
                for (int testcase = 0; testcase < 1; testcase++) { //Test 1 case for each dimension
                        cout << "Test the number of data " << n << " and case " << testcase + 1 << endl;
			double a = 1./(double)n;

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

                        for (int i = 0; i < n; i++) {
                                std::vector<double> one(n, 0.0);
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

			system_clock::time_point start_time_kernel = system_clock::now(); //Measure the time of evaluation

			std::vector<std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel(n, std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > (n, ct_zero_n)); //Kernel matrix

			//Diagonal elements
			for (int i = 0; i < n; i++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_i = dataset[i];
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > element_ii = cc -> EvalMult(data_i, data_i);
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotation = element_ii;
				for (int j = 1; j < d; j++) {
					data_square_rotation = cc -> EvalRotate(data_square_rotation, 1);
					element_ii = cc -> EvalAdd(element_ii, data_square_rotation);
				} //The first element of element_ii is the inner product between data
				kernel[i][i] = element_ii;
			}

			//Non-diagonal elements
			for (int i = 0; i < n-1; i ++) {
				for (int j = i+1; j < n; j++) { //i < j; upper element
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > element_ij = cc -> EvalMult(dataset[i], dataset[j]);
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotation = element_ij;
					for (int k = 1; k < d; k++) {
						data_square_rotation = cc -> EvalRotate(data_square_rotation, 1);
						element_ij = cc -> EvalAdd(element_ij, data_square_rotation);
					} //The first element of element_ij is the inner product between data
					kernel[i][j] = element_ij;
					kernel[j][i] = element_ij; //By symmetry
				}
			}

			system_clock::time_point end_time_kernel = system_clock::now();
                        microseconds micro_kernel = duration_cast<microseconds>(end_time_kernel - start_time_kernel);

                        cout << "Time of making kernel: " << micro_kernel.count()/1000000. << " s" << endl;

			//Calculate the first term: (K(x0, x0) + K(x1, x1) + ...)/n

			system_clock::time_point start_time_eval = system_clock::now(); //Measure the time of evaluation
				
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > first_term = kernel[0][0]; //Initial value
			for (int i = 1; i < n; i++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = kernel[i][i];
				first_term = cc -> EvalAdd(first_term, increment);
			} //The first element of first_term is K(x0, x0) + K(x1, x1) + ...

			first_term = cc -> EvalMult(first_term, a); // The first element of first_term is (K(x0, x0) + K(x1, x1) + ...)/n

			//Calculate the second term: (K(x0, x0) + K(x0, x1) + ... + K(x1, x0) + K(x1, x1) + ...)/(n**2), in other words, the numerator is summation of all element of kernel
	
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > second_term = ct_zero_n; //Initial value

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = kernel[i][j];
					second_term = cc -> EvalAdd(second_term, increment);
				}
			} //THe first element of second_term is K(x0, x0) + K(x0, x1) + ... + K(x1, x0) + K(x1, x1) + ...

			second_term = cc -> EvalMult(second_term, a*a); //The first element of second_term is (K(x0, x0) + K(x0, x1) + ... + K(x1, x0) + K(x1, x1) + ...)/(n**2)

			//Calculate variance
	
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > variance = cc -> EvalSub(first_term, second_term);

                        //Calculate standard deviation
                        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > standard_deviation = sqrt(variance, cc);

			system_clock::time_point end_time_eval = system_clock::now();
			microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);
			std::cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << std::endl;

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, standard_deviation, &result1);
			result1 -> SetLength(1);

			fout << micro_kernel.count()/1000000. << endl;
			fout << micro_eval.count()/1000000. << endl;
		}
	}
	fout.close();
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
        ctxs = cc -> EvalMult(ctxs, 0.001);
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

        float constant = sqrt(1000);

        a = cc -> EvalMult(a, constant);

        return a;
}



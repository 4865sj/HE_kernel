#include "openfhe.h"
#include <time.h>
#include <random>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace std;
using namespace lbcrypto;
using namespace chrono;

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
	int n = 4; //The number of data
	double a = 1./(double)n;

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);
        std::cout << "CKKS scheme is using ring dimension " << cc -> GetRingDimension() <<std::endl << std::endl;

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

        ofstream fout;
        fout.open("variance_kernel.txt"); //Save the results of time

        //test 10 case
        for (int d = 15; d < 30; d++) { //dimension
                for (int testcase = 0; testcase < 1; testcase++) { //Test 1 case for each dimension
                        cout << "Test dimension " << d << " and case " << testcase + 1 << endl;
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

			int rows = n; //The number of rows of kernel
			int cols = n; //The number of columns of kernel
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > val = dataset[0]; //Temprory data
			std::vector< std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel(rows, std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > (cols, val)); //Kernel matrix

			for (int i = 0; i < n; i++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square = cc -> EvalMult(dataset[i], dataset[i]); //(d1*d1, d2*d2, ...)
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = data_square;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotate = data_square;

				for (int j = 1; j < d; j++) {
					data_square_rotate = cc -> EvalRotate(data_square_rotate, 1);
					inner_product = cc -> EvalAdd(inner_product, data_square_rotate);
				} //The first element of inner_product is the inner product between data

				kernel[i][i] = inner_product; //Diagonal element
			} //Calculate the diagonal elements of kernel

			for (int i = 0; i < n-1; i++) {
				for (int j = i+1; j < n; j++) {
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[j]);
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;
					for (int k = 1; k < d; k++) {
						d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
						inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
					} //The first element of inner_product is the inner product between data_i and data_j

					kernel[i][j] = inner_product;
					kernel[j][i] = inner_product; //Non-diagonal element
				}
			}

			system_clock::time_point end_time_kernel = system_clock::now();
			microseconds micro_kernel = duration_cast<microseconds>(end_time_kernel - start_time_kernel);
			std::cout << "Time of making kernel: " << micro_kernel.count()/1000000. << " s" << std::endl;

			//Calculate the first term: (K(x1, x1) + K(x2, x2) + ...)/n

			system_clock::time_point start_time_eval = system_clock::now(); //Measure the time of evaluation
				
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > first_term = kernel[0][0];
			for (int i = 1; i < n; i++) {
				first_term = cc -> EvalAdd(first_term, kernel[i][i]);
			} //The first element of first_term is K(x1, x1) + K(x2, x2) + ...

			first_term = cc -> EvalMult(first_term, a); // The first element of first_term is (K(x1, x1) + K(x2, x2) + ...)/n

			//Calculate the second term: (K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...)/(n**2), in other words, the numerator is summation of all element of kernel
	
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > second_term = cc -> EvalSub(kernel[0][0], kernel[0][0]); //Temporary value

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					second_term = cc -> EvalAdd(second_term, kernel[i][j]);
				}
			} //The first element of second term is K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...

			second_term = cc -> EvalMult(second_term, pow(a, 2)); //The first element of second_term is (K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...)/(n**2)

			//Calculate variance
	
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > variance = cc -> EvalSub(first_term, second_term);

			system_clock::time_point end_time_eval = system_clock::now();
			microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);
			std::cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << std::endl;

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, variance, &result1);
			result1 -> SetLength(1);

			fout << micro_kernel.count()/1000000. << endl;
			fout << micro_eval.count()/1000000. << endl;
		}
	}
	fout.close();
}




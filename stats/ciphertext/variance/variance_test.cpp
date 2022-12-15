#include "openfhe.h"
#include <time.h>
#include <random>
#include <fstream>
#include <cmath>
#include <chrono>
#include <ctime>

using namespace std;
using namespace lbcrypto;
using namespace chrono;

int main() {
	srand(time(NULL));

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
	
        std::cout << "CKKS scheme is using ring dimension " << cc -> GetRingDimension() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using scaling mod size " << parameters.GetScalingModSize() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using security level " << parameters.GetSecurityLevel() <<std::endl << std::endl;
	
        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

	std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

	//test 15 case
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

			//Calculate mean

			system_clock::time_point start_time_eval = system_clock::now(); //Measure the time of evaluation

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > mean = dataset[0];

			for (int i = 1; i < n; i++) {
				mean = cc -> EvalAdd(mean, dataset[i]);
			}

			mean = cc -> EvalMult(mean, a); //Mean vector

			//Calculate variance, (x1-m)**2 + (x2-m)**2 + ...
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > variance = cc -> EvalSub(dataset[0], dataset[0]);

			for (int i = 0; i < n; i++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta = cc -> EvalSub(dataset[i], mean);
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square = cc -> EvalMult(delta, delta);
			        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > distance_square = delta_square;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square_rotate = delta_square;

				for (int j = 1; j < d; j++) {
					delta_square_rotate = cc -> EvalRotate(delta_square_rotate, 1);
					distance_square = cc -> EvalAdd(distance_square, delta_square_rotate);
				} //The first element of distance_square is (x_i - m)**2
				variance = cc -> EvalAdd(variance, distance_square);
			} //The first element of variance is the summation of (x_i - m)**2

			variance = cc -> EvalMult(variance, a);

			system_clock::time_point end_time_eval = system_clock::now();
	                microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);

			cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, variance, &result1);
			result1 -> SetLength(1);

		}
	}
}



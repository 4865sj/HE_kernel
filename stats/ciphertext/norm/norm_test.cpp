#include "openfhe.h"
#include <time.h>
#include <random>
#include <fstream>
#include <chrono>
#include <ctime>

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
	int n = 1; //The number of data

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);
        cout << "CKKS scheme is using ring dimension " << cc -> GetRingDimension() <<endl;
        std::cout << "CKKS scheme is using scaling mod size " << parameters.GetScalingModSize() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using security level " << parameters.GetSecurityLevel() <<std::endl << std::endl;

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
	cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

	std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

	ofstream fout;
	fout.open("norm.txt"); //Save the results of time

	//test 10 case
	for (int d = 14; d < 31; d++) { //dimension
		for (int testcase = 0; testcase < 1; testcase++) { //Test 10 case for each dimension
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

			//Norm
			system_clock::time_point start_time_eval = system_clock::now(); //Measure the time of evaluation

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square;

			x_square = cc -> EvalMult(dataset[0], dataset[0]);

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_sum = x_square;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_rotate = x_square;
		
			for (int i = 1; i < d; i++) {
				x_square_rotate = cc -> EvalRotate(x_square_rotate, 1);
				x_square_sum = cc -> EvalAdd(x_square_sum, x_square_rotate);
			} //The first element of x_square_sum is (x1)**2 + (x2)**2 + ...

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_norm = sqrt(x_square_sum, cc); //The first element is the norm of x

			system_clock::time_point end_time_eval = system_clock::now();
                        microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);

			cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, x_norm, &result1);
			result1 -> SetLength(1);

	
		}	
	}
	fout.close();
}



Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
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

        return a;
}

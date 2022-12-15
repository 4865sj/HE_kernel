#include "openfhe.h"
#include <time.h>
#include <cmath>
#include <fstream>
#include <random>
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

        ofstream fout;
        fout.open("distance.txt"); //Save the results of time

        //test 10 case
        for (int d = 14; d < 30; d++) { //dimension
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

			//Calculate x-y

			system_clock::time_point start_time_eval = system_clock::now(); //Measure the time of evaluation

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta = cc -> EvalSub(dataset[0], dataset[1]);

			//Calculate (x1-y1)**2 + (x2-y2)**2 + ...
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square = cc -> EvalMult(delta, delta);
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > distance_square = delta_square;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square_rotate = delta_square;

			for (int i = 1; i < d; i++) {
				delta_square_rotate = cc -> EvalRotate(delta_square_rotate, 1);
				distance_square = cc -> EvalAdd(distance_square, delta_square_rotate);
			} //The first element of distance_square is (x1-y1)**2 + (x2-y2)**2 + ...

			//Calculate distance
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > distance = sqrt(distance_square, cc); //Distance

                        system_clock::time_point end_time_eval = system_clock::now();
                        microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);

                        cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;
	
			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, distance, &result1);
			result1 -> SetLength(1);

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

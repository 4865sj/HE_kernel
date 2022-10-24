#include "openfhe.h"
#include <time.h>
#include <random>
#include <fstream>
#include <cmath>

using namespace std;
using namespace lbcrypto;

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
        fout.open("variance.txt"); //Save the results of time

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

			float time1 = -clock(); //Measure time of encoding and encryption
	
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

			time1 += clock(); //Measure time of encoding and encryption
			time1 = time1/CLOCKS_PER_SEC;

			//Calculate mean
	
			float time2 = -clock();

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

			time2 += clock();
			time2 = time2/CLOCKS_PER_SEC;

			//Decrypting and Decoding
			Plaintext result1;

			float time3 = -clock(); //Measure time of decoding

			cc -> Decrypt(keys.secretKey, variance, &result1);
			result1 -> SetLength(1);

			time3 += clock(); //Measure time of decoding
			time3 = time3/CLOCKS_PER_SEC;

			fout << time1 << endl;
			fout << time2 << endl;
			fout << time3 << endl;
		}
	}
	fout.close();
}



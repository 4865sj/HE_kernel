#include "openfhe.h"
#include <time.h>

using namespace lbcrypto;

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 8;
	int d = 4; //dimension
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

        //Inputs & Encoding & Encryption
        std::vector<double> x1 = {2.0, 2.1, 2.2, 2.3};
        std::vector<double> x2 = {3.1, 3.2, 3.3, 3.4};
	std::vector<double> x3 = {4.1, 4.2, 4.3, 4.4};
	std::vector<double> x4 = {5.1, 5.2, 5.3, 5.4};

	std::vector<std::vector<double>> dataset_original(4, std::vector<double> (4, 0));
	dataset_original = {x1, x2, x3, x4};

	float time1 = -clock(); //Measure time of encoding and encryption
	
	std::vector<Plaintext> plaintext;

	for (int i = 0; i < n; i++) {
		Plaintext ptxt = cc -> MakeCKKSPackedPlaintext(dataset_original[i]);
		plaintext.push_back(ptxt);
		std::cout << "Input x" << i << ": " << ptxt << std::endl;
	}

        std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dataset;

	for (int i = 0; i < n; i++) {
		auto c = cc -> Encrypt(keys.publicKey, plaintext[i]);
		dataset.push_back(c);
	}

	time1 += clock(); //Measure time of encoding and encryption
	time1 = time1/CLOCKS_PER_SEC;
	std::cout << "Time of encoding and encryption: " << time1 << " s" << std::endl;

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
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square_rotate;

		for (int j = 1; j < d; j++) {
			delta_square_rotate = cc -> EvalRotate(delta_square, j);
			distance_square = cc -> EvalAdd(distance_square, delta_square_rotate);
		} //The first element of distance_square is the square of distance of (x_i - m)
		variance = cc -> EvalAdd(variance, distance_square);
	}

	variance = cc -> EvalMult(variance, a);

	time2 += clock();
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time2 << " s" << std::endl;

	//Decrypting and Decoding
	Plaintext result1;

	float time3 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, variance, &result1);
	result1 -> SetLength(1);

	time3 += clock(); //Measure time of decoding
	time3 = time3/CLOCKS_PER_SEC;
	std::cout << "Time of decoding: " << time3 << " s" << std::endl;

	std::cout << result1 <<  std::endl; //Print variance over encrypted data

	//Calculate original value for testing
	double variance_original = 0.0;
	std::vector<double> mean_original = {0.0, 0.0, 0.0, 0.0};
	for (int i = 0; i< n; i++) {
		mean_original[i] = (x1[i] + x2[i] + x3[i] + x4[i])/4;
	}
	for (int i = 0; i < n; i++) {
		std::vector<double> difference = {0.0, 0.0, 0.0, 0.0};
		for (int j = 0; j < d; j++) {
			difference[j] = dataset_original[i][j] - mean_original[j];
		}
		double difference_square = 0.0;
		for (int j = 0; j < d; j++) {
			difference_square += pow(difference[j], 2);
		}
		variance_original += difference_square;
	}
	variance_original = variance_original/n;
	std::cout << variance_original << std::endl;

}



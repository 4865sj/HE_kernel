#include "openfhe.h"
#include <time.h>

using namespace lbcrypto;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 8;
	int d = 4; //dimension
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

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

        //Inputs & Encoding & Encryption
        std::vector<double> x = {0.3, 0.6, 0.5, 0.7};
        std::vector<double> y = {0.7, 0.4, 0.6, 0.3};

	std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));
	dataset_original = {x, y};

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

	//Calculate x-y
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta = cc -> EvalSub(dataset[0], dataset[1]);

	//Calculate (x1-y1)**2 + (x2-y2)**2 + ...
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square = cc -> EvalMult(delta, delta);
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > distance_square = delta_square;
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta_square_rotate;

	for (int i = 1; i < d; i++) {
		delta_square_rotate = cc -> EvalRotate(delta_square, i);
		distance_square = cc -> EvalAdd(distance_square, delta_square_rotate);
	} //The first element of distance_square is (x1-y1)**2 + (x2-y2)**2 + ...

	//Calculate distance
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > distance = sqrt(distance_square, cc); //Distance

	//Decrypting and Decoding
	Plaintext result1;

	float time3 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, distance, &result1);
	result1 -> SetLength(1);

	time3 += clock(); //Measure time of decoding
	time3 = time3/CLOCKS_PER_SEC;
	std::cout << "Time of decoding: " << time3 << " s" << std::endl;

	std::cout << result1 <<  std::endl; //Print similarity over encrypted data

	//Calculate original value for testing
	double sum = 0.0;
	for (int i = 0; i < d; i++) {
		sum += pow(x[i] - y[i], 2);
	}
	double distance_original = sqrt(sum);
	std::cout << distance_original << std::endl;

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

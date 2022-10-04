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

	float time1 = -clock(); //Measure time of encoding and encryption

        Plaintext ptxt1 = cc -> MakeCKKSPackedPlaintext(x);

        std::cout << "Input x: " << ptxt1 << std::endl;

        auto c1 = cc -> Encrypt(keys.publicKey, ptxt1);

	time1 += clock(); //Measure time of encoding and encryption
	time1 = time1/CLOCKS_PER_SEC;
	std::cout << "Time of encoding and encryption: " << time1 << " s" << std::endl;

	//Norm
	float time2 = -clock(); //Measure time of evaluation

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square;

	x_square = cc -> EvalMult(c1, c1);

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_sum = x_square;
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_rotate;

	for (int i = 1; i < d; i++) {
		x_square_rotate = cc -> EvalRotate(x_square, i);
		x_square_sum = cc -> EvalAdd(x_square_sum, x_square_rotate);
	} //The first element of x_square_sum is (x1)**2 + (x2)**2 + ...

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_norm = sqrt(x_square_sum, cc); //The first element is the norm of x

	time2 += clock(); //Measure time of evaluation
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time2 << " s" << std::endl;

	//Decrypting and Decoding
	Plaintext result1;

	float time3 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, x_norm, &result1);
	result1 -> SetLength(1);

	time3 += clock(); //Measure time of decoding
	time3 = time3/CLOCKS_PER_SEC;
	std::cout << "Time of decoding: " << time3 << " s" << std::endl;

	std::cout << result1 <<  std::endl; //Print x_norm over encrypted data

	//Calculate original value for testing
	double x_norm_square = 0.0;
	for (int i = 0; i < d; i++) {
		x_norm_square += x[i]*x[i];
	}
	double x_norm_original = sqrt(x_norm_square);
	std::cout << x_norm_original << std::endl;

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

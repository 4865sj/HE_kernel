#include "openfhe.h"
#include <time.h>

using namespace lbcrypto;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 8;
	int d = 2; //dimension

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
        cc -> EvalRotateKeyGen(keys.secretKey, {1, -2});

        //Inputs & Encoding & Encryption
        std::vector<double> x = {0.3, 0.6};
        std::vector<double> y = {0.7, 0.4};

	float time1 = -clock(); //Measure time of encoding and encryption

        Plaintext ptxt1 = cc -> MakeCKKSPackedPlaintext(x);
        Plaintext ptxt2 = cc -> MakeCKKSPackedPlaintext(y);

        std::cout << "Input x: " << ptxt1 << std::endl;
        std::cout << "Input y: " << ptxt2 << std::endl;

        auto c1 = cc -> Encrypt(keys.publicKey, ptxt1);
        auto c2 = cc -> Encrypt(keys.publicKey, ptxt2);

	time1 += clock(); //Measure time of encoding and encryption
	time1 = time1/CLOCKS_PER_SEC;
	std::cout << "Time of encoding and encryption: " << time1 << std::endl;

        //Norm of x
	
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

	//Norm of y
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square;

        y_square = cc -> EvalMult(c2, c2);

        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square_sum = y_square;
        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square_rotate;

        for (int i = 1; i < d; i++) {
                y_square_rotate = cc -> EvalRotate(y_square, i);
                y_square_sum = cc -> EvalAdd(y_square_sum, y_square_rotate);
        } //The first element of y_square_sum is (y1)**2 + (y2)**2 + ...

        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_norm = sqrt(y_square_sum, cc); //The first element is the norm of y

	//Inverse of x_norm and y_norm
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_norm_inverse = inverse(x_norm, cc); //The first element is the inverse of x_norm
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_norm_inverse = inverse(y_norm, cc); //The first element is the inverse of y_norm

	//Inner product of x and y
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_dot_y = cc -> EvalMult(c1, c2); //(x1*y1, x2*y2, ...)
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = x_dot_y;
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_dot_y_rotate;

	for (int i = 1; i < d; i++) {
		x_dot_y_rotate = cc -> EvalRotate(x_dot_y, i);
		inner_product = cc -> EvalAdd(inner_product, x_dot_y_rotate);
	} //The first element of inner_product is the inner product between x and y

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product_divided_by_x_norm = cc -> EvalMult(inner_product, x_norm_inverse); //The first element is (x dot y)/|x|
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > similarity = cc -> EvalMult(inner_product_divided_by_x_norm, y_norm_inverse); //The first element is (x dot y)/|x||y|

	time2 += clock(); //Measure time of evaluation
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time2 << std::endl;

	//Decrypting and Decoding
	Plaintext result1;

	float time3 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, similarity, &result1);
	result1 -> SetLength(1);

	time3 += clock(); //Measure time of decoding
	time3 = time3/CLOCKS_PER_SEC;
	std::cout << "Time of decoding: " << time3 << std::endl;

	std::cout << result1 <<  std::endl; //Print similarity over encrypted data

	//Calculate original value for testing
	double dot = 0.0, x_norm_square = 0.0, y_norm_square = 0.0;
	for (int i = 0; i < d; i++) {
		dot += x[i]*y[i];
		x_norm_square += x[i]*x[i];
		y_norm_square += y[i]*y[i];
	}
	double similarity_original = dot/(sqrt(x_norm_square)*sqrt(y_norm_square));
	std::cout << similarity_original << std::endl;

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

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc) {
    auto ct1 = cc->EvalSub(1.0,ct);
    auto ct2 = cc->EvalAdd(ct1,1.0);

    for(int i=0; i<5; i++){
        ct1 = cc->EvalMult(ct1,ct1);
        ct2 = cc->EvalMult(ct2,cc->EvalAdd(ct1,1.0));
    }

    return ct2;
}


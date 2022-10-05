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

	//Make kernel
	
	float time2 = -clock(); //Measure the time of making kernel
	
	int rows = 2; //The number of rows of kernel
	int cols = 2; //The number of columns of kernel
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > val = dataset[0]; //Temprory data
	std::vector< std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel(rows, std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > (cols, val)); //Kernel matrix

	for (int i = 0; i < n; i++) {
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square = cc -> EvalMult(dataset[i], dataset[i]); //(d1*d1, d2*d2, ...)
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = data_square;
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotate;

		for (int j = 1; j < d; j++) {
			data_square_rotate = cc -> EvalRotate(data_square, j);
			inner_product = cc -> EvalAdd(inner_product, data_square_rotate);
		} //The first element of inner_product is the inner product between data

		kernel[i][i] = inner_product; //Diagonal element
	} //Calculate the diagonal elements of kernel
	
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[j]);
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate;
			for (int k = 1; k < d; k++) {
				d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j, k);
				inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
			} //The first element of inner_product is the inner product between data_i and data_j

			kernel[i][j] = inner_product;
			kernel[j][i] = inner_product; //Non-diagonal element
		}
	}

	time2 += clock();
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of making kernel: " << time2 << " s" << std::endl;

	//Calculate denominator
	
	float time3 = -clock(); //Measure time of evaluation

	float time4 = -clock();

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > denominator_square = cc -> EvalMult(kernel[0][0], kernel[1][1]); //The first element is the square of fraction
	time4 += clock();
	time4 = time4/CLOCKS_PER_SEC;
	std::cout << "The time of multiplication once time: " << time4 << " s" << std::endl;

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > denominator = sqrt(denominator_square, cc); //The first element is fraction

	//Inverse of denominator
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > denominator_inverse = inverse(denominator, cc); //The first element is the inverse of fraction

	//Calculate numerator
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > numerator = kernel[0][1]; //numerator

	//Calculate similarity
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > similarity = cc -> EvalMult(denominator_inverse, numerator);

	time3 += clock(); //Measure time of evaluation
	time3 = time3/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time3 << " s" << std::endl;

	//Decrypting and Decoding
	Plaintext result1;

	float time5 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, similarity, &result1);
	result1 -> SetLength(1);

	time5 += clock(); //Measure time of decoding
	time5 = time5/CLOCKS_PER_SEC;
	std::cout << "Time of decoding: " << time5 << " s" << std::endl;

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


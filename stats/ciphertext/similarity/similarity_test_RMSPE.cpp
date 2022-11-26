#include "openfhe.h"
#include <fstream>
#include <cmath>

using namespace lbcrypto;
using namespace std;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 50;
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

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, {1, 2, 3, 4});

        //test 10 case
        for (int d = 14; d < 30; d++) { //dimension
                for (int testcase = 0; testcase < 1; testcase++) { //Test 1 case for each dimension
                        cout << "Test dimension " << d << " and case " << testcase + 1 << endl;
                        //Generating dataset
                        std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));
			
			for (int z = 0; z < n; z++) {
                                   std::vector<double> x;
                                   double data = (double)z;
                                   int count = 0;
                                   while (count < d) {
                                           x.push_back(data);
                                           data += 1;
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
        
			//Norm of x	
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square;

			x_square = cc -> EvalMult(dataset[0], dataset[0]);

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_sum = x_square;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_square_rotate = x_square;

			for (int i = 1; i < d; i++) {
				x_square_rotate = cc -> EvalRotate(x_square_rotate, 1);
				x_square_sum = cc -> EvalAdd(x_square_sum, x_square_rotate);
			} //The first element of x_square_sum is (x1)**2 + (x2)**2 + ...

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_norm = sqrt(x_square_sum, cc); //The first element is the norm of x

			//Norm of y
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square;

		        y_square = cc -> EvalMult(dataset[1], dataset[1]);

		        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square_sum = y_square;
		        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_square_rotate = y_square;

		        for (int i = 1; i < d; i++) {
		                y_square_rotate = cc -> EvalRotate(y_square_rotate, 1);
		                y_square_sum = cc -> EvalAdd(y_square_sum, y_square_rotate);
		        } //The first element of y_square_sum is (y1)**2 + (y2)**2 + ...
 			
		        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_norm = sqrt(y_square_sum, cc); //The first element is the norm of y

			//Inverse of x_norm and y_norm
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_norm_inverse = inverse(x_norm, cc); //The first element is the inverse of x_norm
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > y_norm_inverse = inverse(y_norm, cc); //The first element is the inverse of y_norm
			
			//Inner product of x and y
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_dot_y = cc -> EvalMult(dataset[0], dataset[1]); //(x1*y1, x2*y2, ...)
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = x_dot_y;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > x_dot_y_rotate = x_dot_y;
			
			for (int i = 1; i < d; i++) {
				x_dot_y_rotate = cc -> EvalRotate(x_dot_y_rotate, 1);
				inner_product = cc -> EvalAdd(inner_product, x_dot_y_rotate);
			} //The first element of inner_product is the inner product between x and y

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product_divided_by_x_norm = cc -> EvalMult(inner_product, x_norm_inverse); //The first element is (x dot y)/|x|
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > similarity = cc -> EvalMult(inner_product_divided_by_x_norm, y_norm_inverse); //The first element is (x dot y)/|x||y|

			//Decrypting and Decoding
			Plaintext result1;

			cc -> Decrypt(keys.secretKey, similarity, &result1);
			result1 -> SetLength(1);

			cout << "similarity: " << result1 << endl;

		}
	}
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
        ctxs = cc -> EvalMult(ctxs, 0.0001);
	auto a = ctxs;
        auto b = cc->EvalSub(ctxs, 1.0);

        for (int i = 0; i<10; i++) {
                auto b_half = cc -> EvalMult(b, -0.5); // -b/2
                a = cc->EvalMult(a, cc->EvalAdd(b_half, 1.0)); // a = a*(-b/2 + 1)
                auto tmp = cc -> EvalSub(b, 3.0); // b-3
                auto b_quater = cc -> EvalMult(tmp, 0.25); // (b-3)/4
                auto b_square = cc -> EvalMult(b, b); // b**2
                b = cc -> EvalMult(b_square, b_quater); // b = (b**2)*((b-3)/4)
        }
	
	float constant = sqrt(10000);

	a = cc -> EvalMult(a, constant);

        return a;
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc) {
	ct = cc -> EvalMult(ct, 0.001);
	auto ct1 = cc->EvalSub(1.0,ct);
	auto ct2 = cc->EvalAdd(ct1,1.0);

    
	for(int i=0; i<10; i++){
		ct1 = cc->EvalMult(ct1,ct1);
		ct2 = cc->EvalMult(ct2,cc->EvalAdd(ct1,1.0));
    	}

	ct2 = cc -> EvalMult(ct2, 0.001);

	return ct2;
}


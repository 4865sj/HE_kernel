#include "openfhe.h"
#include <time.h>

using namespace lbcrypto;
using namespace std;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 512;
	int n = 1; //The number of data

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
        std::uniform_real_distribution<> dist(10, 20); //Random distribution

        ofstream fout;
        fout.open("norm_kernel.txt"); //Save the results of time

        //test 100 case
        for (int d = 4; d < 401; d = d + 4) { //dimension
                for (int testcase = 0; testcase < 10; testcase++) { //Test 10 case for each dimension
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

			//Make kernel

			float time2 = -clock(); //Measure the time of making kernel

			int rows = n; //The number of rows of kernel
			int cols = n; //The number of columns of kernel
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > val = dataset[0]; //Temprory data
			std::vector< std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel(rows, std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > (cols, val)); //Kernel matrix

			for (int i = 0; i < n; i++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square = cc -> EvalMult(dataset[i], dataset[i]); //(d1*d1, d2*d2, ...)
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = data_square;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotate = data_square;

				for (int j = 1; j < d; j++) {
					data_square_rotate = cc -> EvalRotate(data_square_rotate, 1);
					inner_product = cc -> EvalAdd(inner_product, data_square_rotate);
				} //The first element of inner_product is the inner product between data

				kernel[i][i] = inner_product; //Diagonal element
			} //Calculate the diagonal elements of kernel

			for (int i = 0; i < n; i++) {
				for (int j = i+1; j < n; j++) {
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[j]);
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;
					for (int k = 1; k < d; k++) {
						d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
						inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
					} //The first element of inner_product is the inner product between data_i and data_j

					kernel[i][j] = inner_product;
					kernel[j][i] = inner_product; //Non-diagonal element
				}
			}

			time2 += clock();
			time2 = time2/CLOCKS_PER_SEC;

  			//Measure the time of multiplication once time;parallel
		        float time3 = -clock();
		        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > test = cc -> EvalMult(dataset[0], dataset[0]);
			time3 += clock();
			time3 = time3/CLOCKS_PER_SEC;

			//Norm
			float time4 = -clock();

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > norm = sqrt(kernel[0][0], cc);

			time4 += clock();
			time4 = time4/CLOCKS_PER_SEC;

			//Decrypting and Decoding
			Plaintext result1;

			float time5 = -clock(); //Measure time of decoding

			cc -> Decrypt(keys.secretKey, norm, &result1);
			result1 -> SetLength(1);

			time5 += clock(); //Measure time of decoding
			time5 = time5/CLOCKS_PER_SEC;

			fout << time1 << endl;
			fout << time2 << endl;
			fout << time3 << endl;
			fout << time4 << endl;
			fout << time5 << endl;
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

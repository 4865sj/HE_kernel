#include "openfhe.h"

using namespace lbcrypto;
using namespace std;

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 8;
        int d = 4; //dimension
        int n = 4; //The number of data
	int nd = max(n, d);
        double a = 1./(double)n;

	CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);

        auto keys = cc -> KeyGen();

        cc -> EvalMultKeyGen(keys.secretKey);
	
	//Rotation keys
	vector<int> ro;
	for (int i = 0; i < nd; i++) {
		ro.push_back(i+1);
		ro.push_back(-(i+1));

	}
        cc -> EvalRotateKeyGen(keys.secretKey, ro);

        //Inputs & Encoding & Encryption
        std::vector<double> x1 = {2.0, 2.1, 2.2, 2.3};
        std::vector<double> x2 = {3.1, 3.2, 3.3, 3.4};
        std::vector<double> x3 = {4.1, 4.2, 4.3, 4.4};
        std::vector<double> x4 = {5.1, 5.2, 5.3, 5.4};

	std::vector<std::vector<double>> dataset_original(4, std::vector<double> (4, 0));
        dataset_original = {x1, x2, x3, x4};

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

	//Generating zero vector; its dimension is d
	vector<double> zero_d_original(d, 0.0);
	Plaintext ptxt_zero_d = cc -> MakeCKKSPackedPlaintext(zero_d_original);
	auto ct_zero_d = cc -> Encrypt(keys.publicKey, ptxt_zero_d);

	//Generating zero vector; its dimension is n
	vector<double> zero_n_original(n, 0.0);
	Plaintext ptxt_zero_n = cc -> MakeCKKSPackedPlaintext(zero_n_original);
	auto ct_zero_n = cc -> Encrypt(keys.publicKey, ptxt_zero_n);

	//Generating one vectors; its dimension is d
	vector<Plaintext> one_d_plaintext;
	for (int i = 0; i < d; i++) {
		vector<double> one_d = zero_d_original;
		one_d[i] = 1.0;
		Plaintext ptxt_one_d = cc -> MakeCKKSPackedPlaintext(one_d);
		one_d_plaintext.push_back(ptxt_one_d);
	}

	//Making kernel; Scalar multiplication kernel
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_zero(n, ct_zero_n);
	vector<vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel_vector(d, kernel_zero); //The vector of kernels; the number of elements is d

	for (int i = 0; i < n; i++) {//Select data_i
		auto data_i = dataset[i];
		for (int j = 0; j < n; j++) {//Select data_j
			auto data_j = dataset[j];
			auto increment = cc -> EvalMult(data_i, data_j); //For each dimension d, its element is the increment of i,j-element of kernel_d; For example, increment[k] = kernel_k(i, j)

			for (int k = 0; k < d; k++) {//Select dimension
				auto kernel_k = kernel_vector[k];
				auto kernel_k_i = kernel_k[i]; //i-th row of kernel_k
				auto increment_k = cc -> EvalMult(increment, one_d_plaintext[k]); //Delete except k-th element

				if (k == j) {
					kernel_k_i = cc -> EvalAdd(kernel_k_i, increment_k);
				} else {
					increment_k = cc -> EvalRotate(increment_k, k-j);
					kernel_k_i = cc -> EvalAdd(kernel_k_i, increment_k);
				}
				kernel_k[i] = kernel_k_i; //Update i-th row of kernel_k
				kernel_vector[k] = kernel_k; //Update kernel_k
			}
		}
	}

	//Calculate mean
	
	auto sum_all_kernel = ct_zero_n; //Initial value
	for (int i = 0; i < n; i++) {//Select kernel
		auto kernel_i = kernel_vector[i];
		for (int j = 0; j < n; j++) {//Select row
			sum_all_kernel = cc -> EvalAdd(sum_all_kernel, kernel_i[j]);
		}
	}


	auto mean_square = sum_all_kernel; //Initial value
	for (int i = 1; i < n; i++) {
		sum_all_kernel = cc -> EvalRotate(sum_all_kernel, 1);
		mean_square = cc -> EvalAdd(mean_square, sum_all_kernel);
	} //The first element of mean_square is the summation of all elements of kernels

	mean_square = cc -> EvalMult(mean_square, a*a); //Divide by n^2
	
	//Decryption & Decoding
	Plaintext result;
	cc -> Decrypt(keys.secretKey, mean_square, &result);
	result -> SetLength(1);

	cout << "The square of mean is " << result << endl;

	//Original value
	vector<double> mean_vector_original(d, 0.); //Initial value
	for (int i = 0; i < n; i++) {//Select data
		vector<double> data_i_original = dataset_original[i];
		for (int j = 0; j < d; j++) {//Select dimension
			mean_vector_original[j] = mean_vector_original[j] + data_i_original[j]*a; //a is 1/n
		}
	}
	
	double mean_square_original = 0.; //Initial value

	for (int i = 0; i < d; i++) {
		mean_square_original = mean_square_original + pow(mean_vector_original[i], 2);
	}

	cout << "The original value is " << mean_square_original << endl;
		


}



			



					




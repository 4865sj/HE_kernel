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

	//Make kernel

	float time2 = -clock(); //Measure the time of making kernel

	int rows = n; //The number of rows of kernel
	int cols = n; //The number of columns of kernel
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

        //Measure the time of multiplication once time;parallel
        float time3 = -clock();
        Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > test = cc -> EvalMult(dataset[0], dataset[0]);
	time3 += clock();
	time3 = time2/CLOCKS_PER_SEC;
	std::cout << "The time of multiplication once time: " << time3 << " s" << std::endl;

	//Calculate the first term: (K(x1, x1) + K(x2, x2) + ...)/n
	
	float time4 = -clock(); //Measure the time of evaluation

	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > first_term = kernel[0][0];
	for (int i = 1; i < n; i++) {
		first_term = cc -> EvalAdd(first_term, kernel[i][i]);
	} //The first element of first_term is K(x1, x1) + K(x2, x2) + ...

	first_term = cc -> EvalMult(first_term, a); // The first element of first_term is (K(x1, x1) + K(x2, x2) + ...)/n

	//Calculate the second term: (K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...)/(n**2), in other words, the numerator is summation of all element of kernel
	
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > second_term = cc -> EvalSub(kernel[0][0], kernel[0][0]); //Temporary value

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			second_term = cc -> EvalAdd(second_term, kernel[i][j]);
		}
	} //The first element of second term is K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...

	second_term = cc -> EvalMult(second_term, pow(a, 2)); //The first element of second_term is (K(x1, x1) + K(x1, x2) + ... + K(x2, x1) + K(x2, x2) + ...)/(n**2)

	//Calculate variance
	
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > variance = cc -> EvalSub(first_term, second_term);

	time4 += clock();
	time4 = time4/CLOCKS_PER_SEC;
	std::cout << "Time of evaluation: " << time4 << " s" << std::endl;

	//Decrypting and Decoding
	Plaintext result1;

	float time5 = -clock(); //Measure time of decoding

	cc -> Decrypt(keys.secretKey, variance, &result1);
	result1 -> SetLength(1);

	time5 += clock(); //Measure time of decoding
	time5 = time5/CLOCKS_PER_SEC;
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




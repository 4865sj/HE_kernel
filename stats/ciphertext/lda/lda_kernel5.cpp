#include "openfhe.h"
#include <time.h>
#include <cmath>

using namespace lbcrypto;
using namespace std;

double inner_product(std::vector<double> a, std::vector<double> b);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
diagonalize(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > 
transpose(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
matrix_vector_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxt, int n, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
add_many(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, CryptoContext<DCRTPoly> cc);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > 
matrix_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, int n, CryptoContext<DCRTPoly> cc);

vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q);

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >>
matrix_inverse(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, double dv, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
power_method(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> c, int n, CryptoContext<DCRTPoly> cc);

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > 
inverse(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ct, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 70;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
        int d = 4; //The number of x; dimension
        int n = 5; //The number of data
	int nd = max(n, d);
        int c = 3; //The number of classes
        std::vector<int> classes = {2, 2, 1}; //The number of data per class; In this case, The number of data in first class is 2

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

	vector<int> ro;
	for (int i = 1; i < nd + 1; i++) {
		ro.push_back(i);
		ro.push_back(-i);
	}

        cc -> EvalMultKeyGen(keys.secretKey);
        cc -> EvalRotateKeyGen(keys.secretKey, ro);

        //Inputs & Encoding & Encryption
        std::vector<double> x1 = {0.2, 0.3, 0.4, 0.5};
        std::vector<double> x2 = {0.6, 0.5, 0.4, 0.3};
        std::vector<double> x3 = {0.2, 0.5, 0.3, 0.1};
        std::vector<double> x4 = {0.2, 0.4, 0.6, 0.8};
        std::vector<double> x5 = {0.1, 0.4, 0.7, 0.9};
        std::vector<double> zero(d, 0.0); //Zero vector
	std::vector<double> zero_n(n, 0.0);

	//Generate one vectors; {1, 0, 0, ... }, {0, 1, 0, ... }
	std::vector<Plaintext> ptxt_one;
	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > vector_one;

	for (int i = 0; i < n; i++) {
		std::vector<double> one(n, 0.0);
		one[i] = 1.0;
		Plaintext one_tmp = cc -> MakeCKKSPackedPlaintext(one);
		ptxt_one.push_back(one_tmp);
		auto ct_one = cc -> Encrypt(keys.publicKey, one_tmp);
		vector_one.push_back(ct_one);
	}

	//Generate various zero vectors
	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > vector_zeros;
	for (int i = 0; i < c; i++) {
		std::vector<double> zeros(classes[i], 0.0);
		Plaintext zeros_tmp = cc -> MakeCKKSPackedPlaintext(zeros);
		auto ct_zeros = cc -> Encrypt(keys.publicKey, zeros_tmp);
		vector_zeros.push_back(ct_zeros);
	}

	//Generate all_ones vector
	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > vector_all_ones;
	for (int i = 0; i < c; i++) {
		std::vector<double> all_ones(classes[i], 1.0);	
		Plaintext all_ones_tmp = cc -> MakeCKKSPackedPlaintext(all_ones);
		auto ct_all_ones = cc -> Encrypt(keys.publicKey, all_ones_tmp);
		vector_all_ones.push_back(ct_all_ones);
	}

	//Generate all_n_i vector
	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > vector_all_n_i;
	for (int i = 0; i < c; i++) {
		vector<double> all_n_i(classes[i], 1.0/(double)classes[i]);
		Plaintext all_n_i_tmp = cc -> MakeCKKSPackedPlaintext(all_n_i);
		auto ct_all_n_i = cc -> Encrypt(keys.publicKey, all_n_i_tmp);
		vector_all_n_i.push_back(ct_all_n_i);
	}

        Plaintext ptxt_zero = cc -> MakeCKKSPackedPlaintext(zero); //Encoding zero vector
        auto ct_zero = cc -> Encrypt(keys.publicKey, ptxt_zero); //Encrypting zero vector

	Plaintext ptxt_zero_n = cc -> MakeCKKSPackedPlaintext(zero_n);
	auto ct_zero_n = cc -> Encrypt(keys.publicKey, ptxt_zero_n);

        std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));
        dataset_original = {x1, x2, x3, x4, x5};

        float time1 = -clock();

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

	std::vector<std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > > kernel_vector;

	int count = 0;

	for (int r = 0; r < c; r++) {
		int number = classes[r];

		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel(n, vector_zeros[r]); //Kernel matrix

		//First, calculate for duplicate data
		//Diagonal element of duplicate data; kernel[n_i][0], kernel[n_i + 1][1], ...
		for (int i = count; i < count + number; i++) {
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i]; //ith row of kernel
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square = cc -> EvalMult(dataset[i], dataset[i]); //(d1*d1, d2*d2, ...)
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = data_square;
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > data_square_rotate = data_square;

			for (int j = 1; j < d; j++) {
				data_square_rotate = cc -> EvalRotate(data_square_rotate, 1);
				inner_product = cc -> EvalAdd(inner_product, data_square_rotate);
			} //The first element of inner_product is the inner product between data

			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element = cc -> EvalMult(inner_product, ptxt_one[0]); //Reset data except the first element
			if (i == count) {
				ith_row = cc -> EvalAdd(ith_row, i_j_element);
			} else {
				i_j_element = cc -> EvalRotate(i_j_element, -(i - count));	
				ith_row = cc -> EvalAdd(ith_row, i_j_element);
			}
			kernel[i] = ith_row; //Update kernel[i][i - count]
		} //Calculate the diagonal elements of kernel of duplicate data
		
		//Non-diagonal duplicate data; kernel[n_i][1], kernel[n_i][2], ... , kernel[n_i + 1][2], kernel[n_i + 1][3], ...
		for (int i = count; i < count + number - 1; i++) {
			for (int j = i + 1; j < count + number; j++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i];
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[j]);
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;
				for (int k = 1; k < d; k++) {
					d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
					inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
				} //The first element of inner_product is the inner product between data_i and data_j
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element = cc -> EvalMult(inner_product, ptxt_one[0]);
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element_rotate = cc -> EvalRotate(i_j_element, -(j - count));
				
				ith_row = cc -> EvalAdd(ith_row, i_j_element_rotate);
				kernel[i] = ith_row; //Update kernel[i][j - count]

				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > jth_row = kernel[j]; //By using symmetry, we can update kernel[j][i - count]
				
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > j_i_element;

				if (i == count) {
					j_i_element = i_j_element;
				} else {
					j_i_element = cc -> EvalRotate(i_j_element, -(i - count));
				}

				jth_row = cc -> EvalAdd(jth_row, j_i_element);
				kernel[j] = jth_row; //Update kernel[j][i - count]
			}
		} //Calculate the non-diagonal elements of kernel of duplicate data
		
		//Next, calculate remainders
		//Upper remainders
		for (int i = 0; i < count; i++) {
			for (int j = 0; j < number; j++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i];
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[count + j]);
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;
				for (int k = 1; k < d; k++) {
					d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
					inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
				} //The first element of inner_product is the inner product between data_i and data_j
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element = cc -> EvalMult(inner_product, ptxt_one[0]);
				
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element_rotate;

				if (j == 0) {
					i_j_element_rotate = i_j_element;
				} else {
					i_j_element_rotate = cc -> EvalRotate(i_j_element, -j);
				}
			
                                ith_row = cc -> EvalAdd(ith_row, i_j_element_rotate);
                                kernel[i] = ith_row; //Update kernel[i][j]
			}
		}
		
		//Lower remainders
		for (int i = count + number; i < n; i++) {
			for (int j = 0; j < number; j++) {
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ith_row = kernel[i];
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j = cc -> EvalMult(dataset[i], dataset[count + j]);
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > inner_product = d_i_dot_d_j;
                                Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > d_i_dot_d_j_rotate = d_i_dot_d_j;
                                for (int k = 1; k < d; k++) {
                                        d_i_dot_d_j_rotate = cc -> EvalRotate(d_i_dot_d_j_rotate, 1);
                                        inner_product = cc -> EvalAdd(inner_product, d_i_dot_d_j_rotate);
                                } //The first element of inner_product is the inner product between data_i and data_j
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element = cc -> EvalMult(inner_product, ptxt_one[0]);
                                
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > i_j_element_rotate;

				if (j == 0) {
					i_j_element_rotate = i_j_element;
				} else {
					i_j_element_rotate = cc -> EvalRotate(i_j_element, -j);
				}
	
                                ith_row = cc -> EvalAdd(ith_row, i_j_element_rotate);
                                kernel[i] = ith_row; //Update kernel[i][j]
                        }
                }
		kernel_vector.push_back(kernel);
		count = count + classes[r];
	}

	time2 += clock();
	time2 = time2/CLOCKS_PER_SEC;
	std::cout << "Time of making kernel: " << time2 << " s" << std::endl;

	//Evaluation
	//Caculate m_i
	
	float time3 = -clock();
	
	float time_m_i = -clock();

	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > m_vector;

	for (int i = 0; i < c; i++) { //Select class
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i = kernel_vector[i]; //n x n_i matrix
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i_trans = transpose(kernel_i, ptxt_one, n, cc); //n_i x n matrix
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i_diag = diagonalize(kernel_i_trans, ptxt_one, n, cc);
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > m_i = matrix_vector_mult(kernel_i_diag, vector_all_ones[i], n, cc);
		m_i = cc -> EvalMult(m_i, 1.0/(double)classes[i]);
		m_vector.push_back(m_i);
	}
		
	cout << "m_i complete" << endl;
	time_m_i += clock();
	time_m_i = time_m_i/CLOCKS_PER_SEC;
	cout << "Time of m_i part: " << time_m_i << " s" << endl;

	//Calculate M
	
	float time_M = -clock();

	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > M(n, ct_zero_n);
	for (int i = 0; i < c - 1; i++) { //Select m_i
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > m_i = m_vector[i];
		for (int j = i + 1; j < c; j++) { //Select m_j
			Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > m_j = m_vector[j];
			//First, calculate diagonal elements
			for (int k = 0; k < n; k++) { //Select row
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta = m_i - m_j;
				Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = cc -> EvalMult(delta, delta); //kth element is increment of M[k][k]
				increment = cc -> EvalMult(increment, ptxt_one[k]); //Reset data except kth element
				M[k] = cc -> EvalAdd(M[k], increment); //Update M[k][k]
			} //Diagonal elements

			//Next, calculate non-diagonal elements
			for (int k = 0; k < n - 1; k++) { //Select row
				for (int l = k + 1; l < n; l++) { //Select column
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta1 = m_i - m_j; //We need kth element
					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta1_rotate = cc -> EvalRotate(delta1, -(l - k)); //delta1_rotate[l] = delta1[k]

					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > delta2 = m_i - m_j; //We need lth element

					Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > increment = cc -> EvalMult(delta1_rotate, delta2); //lth element is increment of M[k][l]
					increment = cc -> EvalMult(increment, ptxt_one[l]); //Reset data except lth element
					M[k] = cc -> EvalAdd(M[k], increment); //Update M[k][l]
					
					increment = cc -> EvalRotate(increment, l - k); //kth element is increment of M[l][k]
					M[l] = cc -> EvalAdd(M[l], increment); // By symmetry, update M[l][k]
				}
			} //Non-diagonal elements
		}
	}
	cout << "M complete" << endl;
	time_M += clock();
	time_M = time_M/CLOCKS_PER_SEC;
	cout << "Time of M part: " << time_M << " s" << endl;

	//Calculate N
	
	float time_N = -clock();

	std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > N(n, ct_zero_n);

	for (int i = 0; i < c; i++) { //Select class
		//Calculate middle matrix
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > middle_matrix(classes[i], vector_zeros[i]);
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > all_n_i = vector_all_n_i[i];
		for (int j = 0; j < classes[i]; j++) {
			middle_matrix[j] = cc -> EvalSub(vector_one[j], all_n_i);
		}

		for (int j = classes[i]; j < n; j++) {
			middle_matrix.push_back(ct_zero_n);
		}

		//First multiplication
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i = kernel_vector[i]; //n x n_i matrix
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i_trans = transpose(kernel_i, ptxt_one, n, cc); //n_i x n matrix
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_i_diag = diagonalize(kernel_i_trans, ptxt_one, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > first_mult = matrix_mult(kernel_i_diag, middle_matrix, n, cc);
		//Second multiplication
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > first_mult_diag = diagonalize(first_mult, ptxt_one, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > increment = matrix_mult(first_mult_diag, kernel_i, n, cc);
		
		for (int j = 0; j < n; j++) {
			N[j] = cc -> EvalAdd(N[j], increment[j]);
		}
	
	} //The output is the tranpose of N		

	cout << "N complete" << endl;
	time_N += clock();
	time_N = time_N/CLOCKS_PER_SEC;
	cout << "Time of N part: " << time_N << " s" << endl;

	//Calculate dominant eigenvector of (N_inverse)*(M)
	N = transpose(N, ptxt_one, n, cc);
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> N_inverse = matrix_inverse(N, n, 20.0, cc); //The transpose of N_inverse
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> N_inverse_diag = diagonalize(N_inverse, ptxt_one, n, cc);
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> M_trans = transpose(M, ptxt_one, n, cc);
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix_trans = matrix_mult(N_inverse_diag, M_trans, n, cc);
	vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix = transpose(Final_matrix_trans, ptxt_one, n, cc);
	Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > eigenvector = power_method(Final_matrix, n, cc);

        time3 += clock();
        time3 = time3/CLOCKS_PER_SEC;
        cout << "Time of evaluation: " << time3 << " s" << endl;

	Plaintext result1;

	cc -> Decrypt(keys.secretKey, eigenvector, &result1);
	result1 -> SetLength(n);

	cout << result1 << endl;

	//LDA_kernel over plaintext for testing
	std::vector<std::vector<std::vector<double>>> kernel_vector_original;
	int count_original = 0;
	for (int i = 0; i < c; i++) {
		int number = classes[i];
		int row = n;
		int col = number;

		std::vector<std::vector<double>> kernel_original(row, std::vector<double> (col, 0.0));

		for (int j = count_original; j < count_original + number; j++) {
			kernel_original[j][j-count_original] = inner_product(dataset_original[j], dataset_original[j]);
		}

		for (int j = count_original; j < count_original + number - 1; j++) {
			for (int k = j + 1; k < count_original + number; k++) {
				kernel_original[j][k-count_original] = inner_product(dataset_original[j], dataset_original[k]);
				kernel_original[k][j-count_original] = inner_product(dataset_original[j], dataset_original[k]);
			}
		}
		
		for (int j = 0; j < count_original; j++) {
			for (int k = 0; k < number; k++) {
				kernel_original[j][k] = inner_product(dataset_original[j], dataset_original[count_original + k]);
			}
		}
		
		for (int j = count_original + number; j < n; j++) {
			for (int k = 0 ; k < number; k++) {
				kernel_original[j][k] = inner_product(dataset_original[j], dataset_original[count_original + k]);
			}
		}
		
		kernel_vector_original.push_back(kernel_original);
		count_original = count_original + number;
	}
	
	std::vector<std::vector<double>> m_vector_original;

	for (int i = 0; i < c; i++) {
		std::vector<double> m_i_original;
		std::vector<std::vector<double>> kernel_i_original = kernel_vector_original[i];
		for (int j = 0; j < n; j++) {
			double increment_original = 0.0;
			for (int k = 0; k < classes[i]; k++) {
				increment_original = increment_original + kernel_i_original[j][k];
			}
			increment_original = increment_original/(double)classes[i];
			m_i_original.push_back(increment_original);
		}
		m_vector_original.push_back(m_i_original);
	}

	std::vector<std::vector<double>> M_original(n, std::vector(n, 0.0));

	for (int i = 0; i < c - 1; i++) {
		std::vector<double> m_i_original = m_vector_original[i];
		for (int j = i + 1; j < c; j++) {
			std::vector<double> m_j_original = m_vector_original[j];
			for (int p = 0; p < n; p++) {
				M_original[p][p] = M_original[p][p] +  pow(m_i_original[p] - m_j_original[p], 2);
			}

			for (int p = 0; p < n - 1; p++) {
				for (int q = p + 1; q < n; q++) {
					M_original[p][q] = M_original[p][q] + (m_i_original[p] - m_j_original[p])*(m_i_original[q] - m_j_original[q]);
					M_original[q][p] = M_original[q][p] + (m_i_original[p] - m_j_original[p])*(m_i_original[q] - m_j_original[q]);
				}
			}
		}
	}

	std::vector<std::vector<double>> N_original(n, std::vector(n, 0.0));

	for (int i = 0; i < c; i++) {
		std::vector<std::vector<double>> kernel_i_original = kernel_vector_original[i];
		std::vector<std::vector<double>> middle_matrix_original(classes[i], std::vector(classes[i], -(1.0/(double)classes[i])));
		for (int j = 0; j < classes[i]; j++) {
			middle_matrix_original[j][j] = 1.0 - (1.0/(double)classes[i]);
		}

		vector<vector<double>> first_mult_original = pmatrix_mult(kernel_i_original, middle_matrix_original, n, classes[i], classes[i]);

		cout << "\n" << endl;
		vector<vector<double>> kernel_i_trans_original(classes[i], vector(n, 0.0));
		for (int p = 0; p < classes[i]; p++) {
			for (int q = 0; q < n; q++) {
				kernel_i_trans_original[p][q] = kernel_i_original[q][p];
			}
		}
		std::vector<std::vector<double>> increment_original = pmatrix_mult(first_mult_original, kernel_i_trans_original, n, classes[i], n);

		for (int p = 0; p < n; p++) {
			for (int q = 0; q < n; q++) {
				N_original[p][q] = N_original[p][q] + increment_original[p][q];
			}
		}
	}

	//test kernel
	cout << "kernel" << endl;
	for (int i = 0; i < c; i++) {
		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > kernel_cp = kernel_vector[i];
		std::vector<std::vector<double>> kernel_pt = kernel_vector_original[i];
		for (int p = 0; p < n; p++) {
			Plaintext result;
			cc -> Decrypt(keys.secretKey, kernel_cp[p], &result);
			result -> SetLength(classes[i]);
			cout << "kernel[" << p << "]: " << result << endl;
			for (int q = 0; q < classes[i]; q++) {
				cout << "kernel" << i << "[" << p << "][" << q << "] pt: " << kernel_pt[p][q] << endl;
			}
			cout << "\n" << endl;
		}
	}

	//test m_i
	cout << "m_i" << endl;
	for (int i = 0; i < c; i++) {
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > m_i_cp = m_vector[i];
		std::vector<double> m_i_pt = m_vector_original[i];
		
		Plaintext result;
		cc -> Decrypt(keys.secretKey, m_i_cp, &result);
		result -> SetLength(n);
		cout << "m_" << i << result << endl;
		for (int p = 0; p < n; p++) {
                        cout << "m_" << i << "[" << p << "] pt: " << m_i_pt[p] <<  endl;
		}
		cout << "\n" << endl;
	}

	//test M
	cout << "M" << endl;
	for (int i = 0; i < n; i++) {
		Plaintext result;
		cc -> Decrypt(keys.secretKey, M[i], &result);
		result -> SetLength(n);
		cout << "M[" << i << "]: " << result << endl;
		for (int j = 0; j < n; j++) {
			cout << "M[" << i << "][" << j << "] pt: " << M_original[i][j] << endl;
		}
		cout << "\n" << endl;
	}

	//test N
	cout << "N" << endl;
        for (int i = 0; i < n; i++) {
		Plaintext result;
		cc -> Decrypt(keys.secretKey, N[i], &result);
		result -> SetLength(n);
		cout << "N[" << i << "]: " << result << endl;
                for (int j = 0; j < n; j++) {
                        cout << "N[" << i << "][" << j << "] pt: " << N_original[i][j] << endl;
                }
		cout << "\n" << endl;
        }
}

double inner_product(std::vector<double> a, std::vector<double> b) {
	int n = a.size();
	int m = b.size();
	double result = 0.0;

	if (n != m) {
		std::cout << "Two input vectors have differenct size!" << std::endl;

		return 0.0;
	}
	else {
		for (int i = 0; i < n; i++) {
			result += a[i]*b[i];
		}

		return result;
	}
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
diagonalize(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > dev;
    for(int i=0; i<n; i++){
        auto tmp = cc->EvalMult(ctxst[i],diag[i]);
        tmpv.push_back(tmp);
    }
    auto tmp = add_many(tmpv,n,cc);
    dev.push_back(tmp);
    tmpv.clear();

    for(int i=1; i<n; i++){
        for(int j=0; j<n; j++){
            auto tmp = cc->EvalMult(ctxst[j],diag[(n+j-i)%n]);
            tmpv.push_back(tmp);
        }
        auto tmp = add_many(tmpv,n,cc);
        dev.push_back(tmp);
        tmpv.clear();
    }

    return dev;
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
transpose(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, vector<Plaintext> diag, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > ctxst;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            auto tmp = cc->EvalMult(ctxs[j],diag[i]);
            if(j < i){
                tmp = cc->EvalRotate(tmp,i-j);
            }
            else if(j > i){
                tmp = cc->EvalRotate(tmp,-j+i);
            }
            tmpv.push_back(tmp);
        }
        auto tmp = add_many(tmpv,n,cc);
        tmpv.clear();
        ctxst.push_back(tmp);
    }

    return ctxst;
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
matrix_vector_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxt, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    for(int j=0; j<n; j++){
        if(j==0){
            auto tmp = cc->EvalMult(dev[j],ctxt);
            tmpv.push_back(tmp);
        }
        else{
            auto tmp = cc->EvalRotate(ctxt,j);
            auto tmp1 = cc->EvalRotate(ctxt,-n+j);
            tmp1 = cc->EvalAdd(tmp,tmp1);
            tmp = cc->EvalMult(dev[j],tmp1);
            tmpv.push_back(tmp);
        }
    }
    auto matrix_vector_mult = add_many(tmpv,n,cc);
    tmpv.clear();

    return matrix_vector_mult;

}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
add_many(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, CryptoContext<DCRTPoly> cc) {

    auto res = cc->EvalAdd(ctxs[0],ctxs[1]);
    if (n > 2){
        for(int i=2; i<n; i++){
            res = cc->EvalAdd(res,ctxs[i]);
        }
    }
    return res;
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > >
matrix_mult(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> dev, vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxst, int n, CryptoContext<DCRTPoly> cc){

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;
    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > matrix_mult;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(j==0){
                auto tmp = cc->EvalMult(dev[j],ctxst[i]);
                tmpv.push_back(tmp);
            }
            else{
                auto tmp = cc->EvalRotate(ctxst[i],j);
                auto tmp1 = cc->EvalRotate(ctxst[i],-n+j);
                tmp1 = cc->EvalAdd(tmp,tmp1);
                tmp = cc->EvalMult(dev[j],tmp1);
                tmpv.push_back(tmp);
            }
        }
        auto tmp = add_many(tmpv,n,cc);
        tmpv.clear();
        matrix_mult.push_back(tmp);
    }

    return matrix_mult;
}

vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q){

    vector<vector<double>> m;
    vector<double> tmp;
    for(int p=0; p<n; p++){
        for(int i=0; i<q; i++){
            double zero = 0.;
            for(int j=0; j<d; j++){
                zero += m1[p][j] * m2[j][i];
            }
            tmp.push_back(zero);
        }
        m.push_back(tmp);
        tmp.clear();
    }

    return m;
}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
power_method(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> c, int n, CryptoContext<DCRTPoly> cc){

    vector<Plaintext> p;
    for(int i=0; i<n; i++){
        Plaintext ptxt;
        p.push_back(ptxt);
    }
    vector<Plaintext> diag;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(1.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag.push_back(p[i]);
        d.clear();
    }

    auto c_diag = diagonalize(c,diag,n,cc);

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    auto tmp = add_many(c,n,cc);
    tmpv.push_back(tmp);
    for(int i=0; i<4; i++){
        auto tmp = matrix_vector_mult(c_diag,tmpv[i],n,cc);
        tmpv.push_back(tmp);
    }

    return tmpv[4];
}

vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >>
matrix_inverse(vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> ctxs, int n, double dv, CryptoContext<DCRTPoly> cc){

    vector<Plaintext> p;
    for(int i=0; i<n; i++){
        Plaintext ptxt;
        p.push_back(ptxt);
    }
    vector<Plaintext> diag;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(1.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag.push_back(p[i]);
        d.clear();
    }

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > tmpv;

    auto ctxst = transpose(ctxs,diag,n,cc);

    auto c = matrix_mult(diagonalize(ctxst,diag,n,cc),ctxs,n,cc);
    auto c_diag = diagonalize(c,diag,n,cc);

    auto tmp = add_many(c,n,cc);
    tmpv.push_back(tmp);
    for(int i=0; i<4; i++){
        auto tmp = matrix_vector_mult(c_diag,tmpv[i],n,cc);
        tmpv.push_back(tmp);
    }

    auto binverse = cc->EvalMult(tmpv[4],1./dv);
    binverse = cc->EvalMult(binverse,diag[0]);
    auto alpha = cc->EvalMult(inverse(binverse,cc),tmpv[3]);
    alpha = cc->EvalMult(alpha,1./dv);
    alpha = cc->EvalMult(alpha,diag[0]);
    tmpv.clear();
    tmpv.push_back(alpha);
    for(int i=1; i<n; i++){
        auto tmp = cc->EvalRotate(alpha,-i);
        tmpv.push_back(tmp);
    }
    alpha = add_many(tmpv,n,cc);

    vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > b;
    for(int i=0; i<n; i++){
        auto tmp = cc->EvalMult(alpha,ctxs[i]);
        b.push_back(tmp);
    }

    vector<Plaintext> diag2;
    for(int i=0; i<n; i++){
        vector<double> d;
        for(int j=0; j<n; j++){
            if(i == j){
                d.push_back(2.);
            }
            else{
                d.push_back(0.);
            }
        }
        p[i] = cc->MakeCKKSPackedPlaintext(d);
        diag2.push_back(p[i]);
        d.clear();
    }
    for(int i=0; i<15; i++){
        auto ab = matrix_mult(diagonalize(ctxst,diag,n,cc),b,n,cc);
        for(int j=0; j<n; j++){
            ab[j] = cc->EvalSub(diag2[j],ab[j]);
        }
        b = matrix_mult(diagonalize(b,diag,n,cc),ab,n,cc);
    }

    return b;
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

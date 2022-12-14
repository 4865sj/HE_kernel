#include "openfhe.h"
#include <time.h>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>
#include <chrono>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

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
	srand(time(NULL));
	
        // Setting
        uint32_t multDepth = 65;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
        int n = 10; //The number of data
        int c = 3; //The number of classes
        std::vector<int> classes = {3, 3, 4}; //The number of data per class; In this case, The number of data in first class is 3

        CCParams<CryptoContextCKKSRNS> parameters;
        parameters.SetMultiplicativeDepth(multDepth);
        parameters.SetScalingModSize(scaleModSize);
        parameters.SetBatchSize(batchSize);

        CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

        cc -> Enable(PKE);
        cc -> Enable(KEYSWITCH);
        cc -> Enable(LEVELEDSHE);


        std::cout << "CKKS scheme is using ring dimension " << cc -> GetRingDimension() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using scaling mod size " << parameters.GetScalingModSize() <<std::endl << std::endl;
        std::cout << "CKKS scheme is using security level " << parameters.GetSecurityLevel() <<std::endl << std::endl;
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dist(0.3, 0.7);

	for (int d = 5; d < 6; d++) {
		cout << "d: " << d << " starts" << endl;
        	// Setting
		int nd = max(n, d);

		auto keys = cc -> KeyGen();

		vector<int> ro;
		for (int i = 1; i < nd + 1; i++) {
			ro.push_back(i);
			ro.push_back(-i);
		}

        	cc -> EvalMultKeyGen(keys.secretKey);
	        cc -> EvalRotateKeyGen(keys.secretKey, ro);

		//Generate dataset
		std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0));
		
		for (int z = 0; z < n; z++) {
			vector<double> x;
			int x_index = 0;
			while (x_index < d) {
				double data = dist(gen);
				x.push_back(data);
				x_index += 1;
			}
			dataset_original[z] = x;
		}
		
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
		
		//Generate epsilon
		std::vector<Plaintext> ptxt_epsilon;
		std::vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > > epsilon;
		for (int i = 0; i < n; i++) {
                        std::vector<double> epsilon_row(n, 0.0);
                        epsilon_row[i] = 0.001;
                        Plaintext epsilon_tmp = cc -> MakeCKKSPackedPlaintext(epsilon_row);
                        ptxt_epsilon.push_back(epsilon_tmp);
                        auto ct_epsilon = cc -> Encrypt(keys.publicKey, epsilon_tmp);
                        epsilon.push_back(ct_epsilon);
                }

		
		vector<double> zero(d, 0.0);
        	Plaintext ptxt_zero = cc -> MakeCKKSPackedPlaintext(zero); //Encoding zero vector
	        auto ct_zero = cc -> Encrypt(keys.publicKey, ptxt_zero); //Encrypting zero vector
	
		vector<double> zero_n(n, 0.0);
		Plaintext ptxt_zero_n = cc -> MakeCKKSPackedPlaintext(zero_n);
		auto ct_zero_n = cc -> Encrypt(keys.publicKey, ptxt_zero_n);

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
 
		//Make kernel

		system_clock::time_point start_time_kernel = system_clock::now(); //Measure the time of making kernel

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

		system_clock::time_point end_time_kernel = system_clock::now();
		microseconds micro_kernel = duration_cast<microseconds>(end_time_kernel - start_time_kernel);
		std::cout << "Time of making kernel: " << micro_kernel.count()/1000000. << " s" << std::endl;

		//Evaluation
		//Caculate m_i

		system_clock::time_point start_time_eval = system_clock::now();
	
		system_clock::time_point start_time_m_i = system_clock::now();

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
		system_clock::time_point end_time_m_i = system_clock::now();
		microseconds micro_m_i = duration_cast<microseconds>(end_time_m_i - start_time_m_i);
		cout << "Time of m_i part: " << micro_m_i.count()/1000000. << " s" << endl;

		//Calculate M
	
		system_clock::time_point start_time_M = system_clock::now();

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
		system_clock::time_point end_time_M = system_clock::now();
                microseconds micro_M = duration_cast<microseconds>(end_time_M - start_time_M);
		cout << "Time of M part: " << micro_M.count()/1000000. << " s" << endl;

		//Calculate N
	
		system_clock::time_point start_time_N = system_clock::now();

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
		system_clock::time_point end_time_N = system_clock::now();
                microseconds micro_N = duration_cast<microseconds>(end_time_N - start_time_N);
		cout << "Time of N part: " << micro_N.count()/1000000. << " s" << endl;

		//Calculate dominant eigenvector of (N_inverse)*(M)
		for (int i = 0; i < n; i++) {
			N[i] = cc -> EvalAdd(N[i], epsilon[i]);
		} //Since N is singular, we need to add epsilon
		N = transpose(N, ptxt_one, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> N_inverse = matrix_inverse(N, n, 20.0, cc); //The transpose of N_inverse
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> N_inverse_diag = diagonalize(N_inverse, ptxt_one, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> M_trans = transpose(M, ptxt_one, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix_trans = matrix_mult(N_inverse_diag, M_trans, n, cc);
		vector<Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >> Final_matrix = transpose(Final_matrix_trans, ptxt_one, n, cc);
		Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > eigenvector = power_method(Final_matrix, n, cc);
	       
 	        system_clock::time_point end_time_eval = system_clock::now();
                microseconds micro_eval = duration_cast<microseconds>(end_time_eval - start_time_eval);

        	cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;

		Plaintext result1;

		cc -> Decrypt(keys.secretKey, eigenvector, &result1);
		result1 -> SetLength(n);

		cout << result1 << endl;

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

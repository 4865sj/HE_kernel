#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <random>
#include <numeric>

using namespace std;
using namespace chrono;

vector<vector<double>> pmatrix_mult(vector<vector<double>> m1,vector<vector<double>> m2, int n, int d, int q);
vector<double> pmatrix_vector_mult(vector<vector<double>> m1,vector<double> m2, int n, int d);
vector<vector<double>> ptranspose(vector<vector<double>> m, int n, int d);
vector<vector<double>> pmatrix_inverse(vector<vector<double>> A);
vector<double> power_method(vector<vector<double>> A);
double inner_product(std::vector<double> a, std::vector<double> b);
double calculator_mean(std::vector<double> input);
double calculator_standard_deviation(std::vector<double> input);

int main() {
	//Setting
        int n = 10; //The number of data
        int c = 3; //The number of classes
        std::vector<int> classes = {3, 3, 4}; //The number of data per class; In this case, The number of data in first class is 3

	std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.3, 0.7); //Random distribution

        int trials = 20000;

        for (int d = 4; d < 16; d++) { //Dimension

                std::vector<double> result;
                cout << "d: " << d << " starts with trials " << trials << endl;
                for (int testcase = 0; testcase < trials; testcase++) {
                        std::vector<std::vector<double>> dataset_original(n, std::vector<double> (d, 0.0));

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

			system_clock::time_point start_time_eval = system_clock::now();

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
		                vector<vector<double>> middle_matrix(classes[i], vector<double> (classes[i], 0.));
				for (int j = 0; j < classes[i]; j++) {
					for (int k = 0; k < classes[i]; k++) {
						if (j == k) {
							middle_matrix[j][k] = 1.0 - 1./(double)classes[i];
						} else {
							middle_matrix[j][k] = 01./(double)classes[i];
						}
					}
				}
				vector<vector<double>> kernel_i_original_trans = ptranspose(kernel_i_original, n, classes[i]);
				vector<vector<double>> N_i = pmatrix_mult(kernel_i_original, middle_matrix, n, classes[i], classes[i]);
				N_i = pmatrix_mult(N_i, kernel_i_original_trans, n, classes[i], n);
		        }

			vector<vector<double>> N_original_inverse = pmatrix_inverse(N_original);
		        vector<vector<double>> Final_matrix_original = pmatrix_mult(N_original_inverse, M_original, n, n, n);

		        vector<double> eigenvector = power_method(Final_matrix_original);
			eigenvector = eigenvector;

			system_clock::time_point end_time_eval = system_clock::now();
                        nanoseconds nano = end_time_eval - start_time_eval;

                        result.push_back(nano.count());
		}
                double mean = calculator_mean(result);
                double standard_deviation = calculator_standard_deviation(result);
                cout << "The mean of time of evaluation: " << mean << " and the standard deviation of time: " << standard_deviation << endl;
	}
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

vector<double> pmatrix_vector_mult(vector<vector<double>> m1,vector<double> m2, int n, int d){

    vector<double> tmp;
    for(int p=0; p<n; p++){
        double zero = 0.;
        for(int i=0; i<d; i++){
           zero += m1[p][i] * m2[i];
        }
        tmp.push_back(zero);
    }

    return tmp;
}

vector<vector<double>> ptranspose(vector<vector<double>> m, int n, int d){

    vector<vector<double>> transpose(d,vector<double>(n,0.));
    for(int i=0; i<d; i++){
        for(int j=0; j<n; j++){
            transpose[i][j] = m[j][i];
        }
    }

    return transpose;
}

double pinverse(double ct){

    double c = 1. - ct;
    double v = 2. - ct;

    for(int i=0; i<5; i++){
        c = c * c;
        v = v * (1. + c);
    }

    return v;
}

vector<vector<double>> pmatrix_inverse(vector<vector<double>> A) {
	//Input: A is nxn matrix
	
	//Step 1) Setting & calculate C
	int n = A.size();
	vector<vector<double>> A_trans = ptranspose(A, n, n);
	vector<vector<double>> C = pmatrix_mult(A, A_trans, n, n, n);

	//Step 2) Iteration part for getting the dominant eigenvalue
	
	//The first procedure of iteration part; x_new[i] = C[i][0] + C[i][1] + ...
	vector<double> x_new(n, 0.0);
	for (int i = 0; i < n; i++) { //Select row
		double increment = 0.0;
		for (int j = 0; j < n; j++) { //Select column
			increment += C[i][j];
		}
		x_new[i] += increment; //Update x_new[i]
	}

	//Other procedures of interation part
	vector<double> x_old = x_new;
	for (int i = 0; i < 4; i++) { //Iteration with 4 times except first procedure
		x_old = x_new;
		x_new = pmatrix_vector_mult(C, x_old, n, n);
	}

	double alpha_denominator = alpha_denominator; //For eliminating error; this variable is not used
	alpha_denominator = pinverse(x_new[0]); //Because of domain of piverse function, we have some limitation for using this function directly. However, we should evaluate these value for comparing time between for ciphertext and plaintext
	double alpha = x_old[0]/x_new[0]; //alpha = x_old[0]/x_new[0] is the inverse of the dominant eigenvalue

	//Step 3) Iteration part for getting inverse matrix of A
	
	//Generating the necessary matrix for calculating
	vector<vector<double>> alpha_I(n, vector<double> (n, 0.0));
	for (int i = 0; i < n; i++) {
		alpha_I[i][i] = alpha; //All diagonal elements are alpha
	}//alpha_I = (alpha)*I

	vector<vector<double>> two_I(n, vector<double> (n, 0.0));
	for (int i = 0; i < n; i++) {
		two_I[i][i] = 2.0; //All diagonal elements are 2
	}//two_I = 2*I

	//Iteration part
	vector<vector<double>> B = pmatrix_mult(alpha_I, A_trans, n, n, n); //Initial value

	for (int i = 0; i < 15; i++) { //Iteration with 15 times
		//First, calculate latter matrix; 2*I - A*B
		vector<vector<double>> latter_matrix = two_I;
		vector<vector<double>> AB = pmatrix_mult(A, B, n, n, n);
		for (int j = 0; j < n; j++) { //Select row
			for (int k = 0; k < n; k++) { //Select column
				latter_matrix[j][k] = latter_matrix[j][k] - AB[j][k];
			}
		}

		//Update B
		B = pmatrix_mult(B, latter_matrix, n, n, n);
	}

	return B; //The output is the inverse of A
}

vector<double> power_method(vector<vector<double>> A) {
        //Input: A is nxn matrix

        //Step 1) Setting
        int n = A.size();

	//Step 2) Iteration part for getting the dominant eigenvalue

        //The first procedure of iteration part; x_new[i] = A[i][0] + A[i][1] + ...
        vector<double> x_new(n, 0.0);
        for (int i = 0; i < n; i++) { //Select row
                double increment = 0.0;
                for (int j = 0; j < n; j++) { //Select column
                        increment += A[i][j];
                }
                x_new[i] += increment; //Update x_new[i]
        }

        //Other procedures of interation part
        vector<double> x_old = x_new;
        for (int i = 0; i < 4; i++) { //Iteration with 4 times except first procedure
                x_old = x_new;
                x_new = pmatrix_vector_mult(A, x_old, n, n);
        }

	return x_new; //The output is the dominant eigenvector
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

double calculator_mean(std::vector<double> input) {
        int n = input.size();
        double sum = accumulate(input.begin(), input.end(), 0.);
        double mean = sum/double(n);
        return mean;
}

double calculator_standard_deviation(std::vector<double> input) {
        int n = input.size();
        double mean = calculator_mean(input);
        double sum = 0.;

        for (int i = 0; i < n; i++) {
                double increment = pow(input[i] - mean, 2);
                sum = sum + increment;
        }

        double standard_deviation = sqrt(sum/double(n));

        return standard_deviation;
}


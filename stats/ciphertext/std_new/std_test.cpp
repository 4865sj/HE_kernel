#include "openfhe.h"
#include <random>
#include <chrono>

using namespace lbcrypto;
using namespace std;
using namespace chrono;

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc);

int main() {
        // Setting
        uint32_t multDepth = 33;
        uint32_t scaleModSize = 50;
        uint32_t batchSize = 64;
        int n = 15; //The number of data
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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(0.5, 1.0); //Random distribution

	for (int d = 1; d < 2; d++) {
		cout << "Test dimension " << d << endl;
		int nd = max(n, d);
		//Rotation keys
		vector<int> ro;
		for (int i = 0; i < nd; i++) {
			ro.push_back(i+1);
			ro.push_back(-(i+1));

		}
        	cc -> EvalRotateKeyGen(keys.secretKey, ro);

        
		//Inputs & Encoding & Encryption
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

                //Generating one vectors; its dimension is n
                vector<Plaintext> one_n_plaintext;
                for (int i = 0; i < n; i++) {
                        vector<double> one_n = zero_n_original;
                        one_n[i] = 1.0;
                        Plaintext ptxt_one_n = cc -> MakeCKKSPackedPlaintext(one_n);
                        one_n_plaintext.push_back(ptxt_one_n);
                }

		//Calculate std
	
		system_clock::time_point start_time = system_clock::now();
		
		//Caculate mean vector
		auto mean_vector = ct_zero_d; //Initial value

		for (int i = 0; i < n; i++) {
			mean_vector = cc -> EvalAdd(mean_vector, dataset[i]);
		}

		mean_vector = cc -> EvalMult(mean_vector, a); //Multiply 1/n

		//Calculate variance vector
		auto variance_vector = ct_zero_d; //Initial value

		for (int i = 0; i < n; i++) {
			auto increment = cc -> EvalSub(mean_vector, dataset[i]);
			increment = cc -> EvalMult(increment, increment);
			variance_vector = cc -> EvalAdd(variance_vector, increment);
		}

		variance_vector = cc -> EvalMult(variance_vector, a*a);

		//Caculate std vector
		auto std_vector = sqrt(variance_vector, cc);

		system_clock::time_point end_time = system_clock::now();
		microseconds micro_eval = duration_cast<microseconds>(end_time - start_time);

		//Decryption & Decoding
		Plaintext result;
		cc -> Decrypt(keys.secretKey, std_vector, &result);
		result -> SetLength(d);
		
		cout << "Time of evaluation: " << micro_eval.count()/1000000. << " s" << endl;
		
	}

}

Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > >
sqrt(Ciphertext<lbcrypto::DCRTPolyImpl<bigintdyn::mubintvec<bigintdyn::ubint<unsigned int> > > > ctxs, CryptoContext<DCRTPoly> cc) {
        ctxs = cc -> EvalMult(ctxs, 0.0001);
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

        float constant = sqrt(10000);

        a = cc -> EvalMult(a, constant);

        return a;
}


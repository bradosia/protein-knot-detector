/*
 * Name        : Protein Knot Detector
 * Author      : Brad Lee
 * Version     : 1.00
 * License     : GNU LGPL v3
 * Description : Detects knots
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916–919 (2000) doi:10.1038/35022623
 */

#ifndef PKD_AMALGAMATED_H
#define PKD_AMALGAMATED_H

// c++17
#include <string>
#include <functional>
#include <iostream>
#include <optional>
#include <memory>

// c
#include <stdio.h>
#include <string.h>

/*
 * PKD = Protein Knot Detector
 */
namespace PKD {

class CommandLineOptions {
public:
// optional can be used as the return type of a factory that may fail
	static std::optional<bool> output_each_iteration(int argc, char **argv);
	static std::optional<std::string> output_type(int argc, char **argv);
	static std::optional<std::string> input_type(int argc, char **argv);
	static std::optional<std::string> input_file(int argc, char **argv);
};

/*
 * our s x 3 matrix
 * s = amino acid chain length
 */
class DoubleMatrix {
public:
	/*
	 * 1D array is used instead of a 2D array to guarantee contiguous
	 * memory is used, so the algorithm will run faster.
	 */
	double* m;
	std::size_t n;
	std::size_t s;
	DoubleMatrix(std::size_t size) {
		s = size;
		n = s * 3;
		m = new double[n];
	}
	DoubleMatrix(const DoubleMatrix&) = delete;
	DoubleMatrix& operator=(const DoubleMatrix&) = delete;
	~DoubleMatrix() {
		delete[] m;
	}
	void printMatrix() {
		for (size_t i = 0; i < n; i += 3) {
			printf("%f %f %f\n", m[i], m[i + 1], m[i + 2]);
		}
	}
};

/*
 * William R. Taylor Knot Detection Algorithm
 */
class TaylorKnotAlgorithm {
private:
	std::unique_ptr<DoubleMatrix> m;
public:
	std::unique_ptr<DoubleMatrix> getMatrix();
	void setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr);
	void smooth(unsigned int nRepeat);
	void smoothAuto();
	/* After about 50 iterations of smoothing,
	 * the knot now may be detected.
	 */
};

std::optional<bool> CommandLineOptions::output_each_iteration(int argc,
		char **argv) {
	bool returnValue = { };
	char * token;
	for (int i = 0; i < argc; i++) {
		token = strtok(argv[i], "=");
		if (strcmp("--output_each_iteration", token) == 0) {
			token = strtok(NULL, "=");
			if (strcmp("true", token) == 0) {
				returnValue = true;
			} else if (strcmp("false", token) == 0) {
				returnValue = false;
			} else {
				printf("Warning: option 'output_each_iteration' invalid");
			}
			break;
		}
	}
	return returnValue;
}

std::optional<std::string> CommandLineOptions::output_type(int argc,
		char **argv) {
	std::string returnValue = { };
	char * token;
	for (int i = 0; i < argc; i++) {
		token = strtok(argv[i], "=");
		if (strcmp("--output_type", token) == 0) {
			token = strtok(NULL, "=");
			if (strcmp("pdb", token) == 0) {
				returnValue = "pdb";
			} else {
				returnValue = "pdb";
				printf("Warning: option 'output_type' invalid");
			}
			break;
		}
	}
	return returnValue;
}

std::optional<std::string> CommandLineOptions::input_type(int argc,
		char **argv) {
	std::string returnValue = { };
	char * token;
	for (int i = 0; i < argc; i++) {
		token = strtok(argv[i], "=");
		if (strcmp("--input_type", token) == 0) {
			token = strtok(NULL, "=");
			if (strcmp("pdb", token) == 0) {
				returnValue = "pdb";
			} else {
				returnValue = "pdb";
				printf("Warning: option 'input_type' invalid");
			}
			break;
		}
	}
	return returnValue;
}

std::optional<std::string> CommandLineOptions::input_file(int argc,
		char **argv) {
	std::string returnValue = { };
	char * token;
	for (int i = 0; i < argc; i++) {
		token = strtok(argv[i], "=");
		if (strcmp("--input_file", token) == 0) {
			printf("WHYYY");
			token = strtok(NULL, "=");
			if (strcmp("", token) != 0) {
				returnValue = token;
			}
			break;
		}
	}
	return returnValue;
}

std::unique_ptr<DoubleMatrix> TaylorKnotAlgorithm::getMatrix() {
	return std::move(m);
}
void TaylorKnotAlgorithm::setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr) {
	m = std::move(matrixPtr);
}
void TaylorKnotAlgorithm::smooth(unsigned int nRepeat = 1) {
	int n = m->n - 3;
	for (unsigned int j = 0; j < nRepeat; j++) {
		for (int i = 3; i < n; i++) {
			m->m[i] = ((m->m[i - 3] + m->m[i + 3]) / 2 + m->m[i]) / 2;
		}
	}
}
void TaylorKnotAlgorithm::smoothAuto() {

}

} // namespace PKD

#endif

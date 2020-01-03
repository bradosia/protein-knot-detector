/*
 * Name        : Protein Knot Detector
 * Author      : Brad Lee
 * Version     : 1.00
 * License     : GNU LGPL v3
 * Description : Detects knots
 *
 * Works Cited:
 * Taylor, W. A deeply knotted protein structure and how it might fold.
 * Nature 406, 916�919 (2000) doi:10.1038/35022623
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

// Copyright (C) 2016 by Doug Baldwin.
// This work is licensed under a Creative Commons Attribution-ShareAlike 4.0 International
// License (http://creativecommons.org/licenses/by-sa/4.0/).
// High and low bounds on t values that are considered to represent a ray intersecting
// a triangle's plane.
#define tFar   10000.0f
#define tNear  0.0000001f
/* CROSS, DOT, and SUB3 Macros for 3-component vectors
 * used in original Moeller and Trumbore algorithm.
 *
 * Citation for macros from Moeller and Trumbore's original code
 * Created on: 1997
 * 		Author: Tomas Moeller -- Chalmers University of Technology
 * 		& Ben Trumbore -- Cornell University
 * 		Title:  Fast, Minimum Storage Ray /Triangle Intersection
 * 		Source: http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
 */

#define CROSS(res, v1, v2)\
		res[0]=v1[1]*v2[2]-v1[2]*v2[1];\
		res[1]=v1[2]*v2[0]-v1[0]*v2[2];\
		res[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB3(res,v1,v2)\
		res[0]=v1[0]-v2[0];\
		res[1]=v1[1]-v2[1];\
		res[2]=v1[2]-v2[2];
/*
 * CRITICAL ALGORITHM
 * Function calls would create too much overhead so we use #define
 * for maximum computational efficiency
 */
void TaylorKnotAlgorithm::smooth(unsigned int nRepeat = 1) {
	double *x = m->m; // x is an alias for the vertex matrix
	double *v0, *v1, *v2, *v0a, *v1a, *v2a, *rayOrigin, *rayDirection;
	double v1p[3];
	/* v# are the operated vertexes
	 * v#a are the committed vertexes
	 * v#p are the prime vertexes (vertex after move)
	 */
	int n = m->n - 3;
	for (unsigned int j = 0; j < nRepeat; j++) {
		for (int i = 3; i < n; i += 3) {
			v0a = x + i - 3;
			v1a = x + i;
			v2a = x + i + 3;
			v1p[0] = ((v0a[0] + v2a[0]) / 2 + v1a[0]) / 2;
			v1p[1] = ((v0a[1] + v2a[1]) / 2 + v1a[1]) / 2;
			v1p[2] = ((v0a[2] + v2a[2]) / 2 + v1a[2]) / 2;
			/* check that the triangles {i'-1,i,i'} and {i;i';i+1}
			 * did not intersect any line segment {j'-1;j'}(j<i) before the move point
			 * or any line {j;j+1}(j>i) following.
			 * implemented with the M�ller�Trumbore intersection algorithm
			 *
			 * This code is repeated 4 times to check intersection for each of two triangles twice
			 * Creating a function for this repeated code would create too much of an overhead
			 * so they code just copy and pasted 4 times.
			 */
			// triangle {i'-1,i,i'} and line {j'-1;j'}(j<i)
			/*v0 = v0a;
			v1 = v1a;
			v2 = v1p;
			for (int k = 3; k < i; k += 3) {
				rayOrigin = x + k - 3; // setup ray
				rayDirection = x + k;
				double edge1[3]; // Find vectors for two edges sharing vertex 0
				SUB3(edge1, v1, v0);
				double edge2[3];
				SUB3(edge2, v2, v0);
				double pvec[3]; // Begin calculating determinant;
				CROSS(pvec, rayDirection, edge2); // also used to calculate U parameter
				const double det = DOT(edge1, pvec); // If determinant is near zero, ray lies in plane of triangle
				if (det > -0.000001f && det < 0.000001f) // No backface culling in this experiment, determinant within "epsilon" as
					continue; // defined in M&T paper is considered 0
				const double inv_det = 1.0f / det;
				double tvec[3]; // Calculate vector from vertex to ray origin
				SUB3(tvec, rayOrigin, v0);
				const double u = DOT( tvec, pvec) * inv_det; // Calculate U parameter and test bounds
				if (u < 0.0f || u > 1.0f)
					continue;
				double qvec[3]; // Prepare to test V parameter
				CROSS(qvec, tvec, edge1);
				const double v = DOT( rayDirection, qvec ) * inv_det; // Calculate V parameter and test bounds
				if (v < 0.0f || u + v >= 1.0f)
					continue;
				const double t = DOT( edge2, qvec ) * inv_det; // Calculate t, final check to see if ray intersects triangle. Test to
				if (t <= tNear || t >= tFar) // see if t > tFar added for consistency with other algorithms in experiment.
					continue;
				goto intersect;
				// intersection found, don't move vertex
				break;// we won't reach this break but it's here anyway
			}
			// triangle {i;i';i+1} and line {j'-1;j'}(j<i)
			v0 = v0a;
			v1 = v1a;
			v2 = v1p;
			for (int k = 3; k < i; k += 3) {
				rayOrigin = x + k - 3;
				rayDirection = x + k;
				double edge1[3];
				SUB3(edge1, v1, v0);
				double edge2[3];
				SUB3(edge2, v2, v0);
				double pvec[3];
				CROSS(pvec, rayDirection, edge2);
				const double det = DOT(edge1, pvec);
				if (det > -0.000001f && det < 0.000001f)
					continue;
				const double inv_det = 1.0f / det;
				double tvec[3];
				SUB3(tvec, rayOrigin, v0);
				const double u = DOT( tvec, pvec) * inv_det;
				if (u < 0.0f || u > 1.0f)
					continue;
				double qvec[3];
				CROSS(qvec, tvec, edge1);
				const double v = DOT( rayDirection, qvec ) * inv_det;
				if (v < 0.0f || u + v >= 1.0f)
					continue;
				const double t = DOT( edge2, qvec ) * inv_det;
				if (t <= tNear || t >= tFar)
					continue;
				goto intersect;
				break;
			}*/
			// both triangles don't intersect, commit vertex move
			v1a[0] = v1p[0];
			v1a[1] = v1p[1];
			v1a[2] = v1p[2];
			/* for algorithm efficiency we use goto instead of if else statements
			 * if and else statements and unnecessary comparisons
			 */
			//intersect: ;
		}
	}
}
void TaylorKnotAlgorithm::smoothAuto() {

}

} // namespace PKD

#endif

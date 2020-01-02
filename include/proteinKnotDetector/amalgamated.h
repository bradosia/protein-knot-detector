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

// c++17
#include <string>
#include <functional>
#include <iostream>
#include <optional>
#include <memory>

// c
#include <stdio.h>
#include <string.h>

/* mmdb 1.25.6.1
 * MMDB is a macromolecular coordinate library,
 * written by Eugene Krissinel primarily for use by CCP4 group.
 * The Coordinate Library is designed to assist CCP4 developers in working with coordinate files.
 * The Library features work with the primary file formats of the Protein Data Bank (PDB),
 * the PDB file format and the mmCIF file format
 * License: GNU LGPL v3
 * Documentation: https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_object.html
 */
#include "mmdb_manager.h"
// MMDB conficts with OCCT so we must rename the definition
#define Abs Absx

/* openCascade (OCCT) 7.4.0
 * OCCT library is designed to be truly modular and extensible, providing C++ classes for:
 * -Basic data structures (geometric modeling, visualization, interactive selection and
 * application specific services);
 * -Modeling algorithms;
 * -Working with mesh (faceted) data;
 * -Data interoperability with neutral formats (IGES, STEP);
 */
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Edge.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>

using namespace std;

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
 * Mediates extraction of data between the MMDB Manager and the Alpha Carbon Matrix.
 * The MMDB Manager handles PDB, CIF, and MMDBF file formats.
 */
class MMDBAndCarbonAlphaMatrix {
private:
	std::unique_ptr<CMMDBManager> MMDBPtr_;
	std::unique_ptr<DoubleMatrix> matrix_;
	int modelId_;
	cpstr chainId_;
public:
	void setMMDBModel(std::unique_ptr<CMMDBManager> MMDBPtr, int modelId,
			cpstr chainId);
	void setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr);
	std::unique_ptr<CMMDBManager> getMMDBModel();
	std::unique_ptr<DoubleMatrix> getMatrix();
	std::unique_ptr<DoubleMatrix> toMatrix();
	std::unique_ptr<CMMDBManager> toMMDB();
};
/*
 * temp
 */
class OCCTManager {
public:
	int writeSTEP(char* path);
};
/*
 * The PDB, CIF, and MMDBF formats are not good for visualizing knots
 * because they do not encode bond data.
 * Instead openCascade is used to display paths between carbon
 * atoms and export as a STEP file to be visualized by a CAD program.
 */
class CarbonAlphaMatrixAndOCCT {
	std::unique_ptr<OCCTManager> OCCTPtr_;
	std::unique_ptr<DoubleMatrix> matrixPtr_;
public:
	void setMMDBModel(std::unique_ptr<CMMDBManager> MMDBPtr, int modelId,
			cpstr chainId);
	void setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr);
	std::unique_ptr<CMMDBManager> getOCCTManager();
	std::unique_ptr<DoubleMatrix> getMatrix();
	std::unique_ptr<DoubleMatrix> toMatrix();
	std::unique_ptr<CMMDBManager> toMMDB();
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
	void smooth();
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

void MMDBAndCarbonAlphaMatrix::setMMDBModel(
		std::unique_ptr<CMMDBManager> MMDBPtr, int modelId, cpstr chainId) {
	MMDBPtr_ = std::move(MMDBPtr);
	modelId_ = modelId;
	chainId_ = chainId;
}

void MMDBAndCarbonAlphaMatrix::setMatrix(
		std::unique_ptr<DoubleMatrix> matrixPtr) {
	matrix_ = std::move(matrixPtr);
}
std::unique_ptr<CMMDBManager> MMDBAndCarbonAlphaMatrix::getMMDBModel() {
	return std::move(MMDBPtr_);
}
std::unique_ptr<DoubleMatrix> MMDBAndCarbonAlphaMatrix::getMatrix() {
	return std::move(matrix_);
}
std::unique_ptr<DoubleMatrix> MMDBAndCarbonAlphaMatrix::toMatrix() {
	int ir, ia, iCA;
	int nResidues, nAtoms;
	std::size_t nCA;
	std::unique_ptr<DoubleMatrix> matrix;
	CChain *chain;
	CResidue **residueTable;
	CAtom **atomTable;

	//  get chain
	printf("Getting Chain...\n");
	chain = MMDBPtr_->GetChain(modelId_, chainId_);
	nCA = 0;
	// get residue table for current chain:
	chain->GetResidueTable(residueTable, nResidues);
	printf("Total Residues In Chain#%s Model#%d: %d\n", chainId_, modelId_,
			nResidues);
	// loop over all residues in current chain:
	for (ir = 0; ir < nResidues; ir++) {
		if (residueTable[ir]) {
			residueTable[ir]->GetAtomTable(atomTable, nAtoms);
			printf("Total Atoms In Residue#%d Chain#%s Model#%d: %d\n", ir,
					chainId_, modelId_, nAtoms);
			for (ia = 0; ia < nAtoms; ia++) {
				// atom exists
				if (atomTable[ia]) {
					if (strcmp((const char*) atomTable[ia]->name, " CA ")
							== 0) {
						nCA++;
					}
				}
			}
		}
	}

// create empty matrix with the carbon atom size
	matrix = std::make_unique<DoubleMatrix>(nCA);
// now get the coordinates
	iCA = 0;
	// get residue table for current chain:
	chain->GetResidueTable(residueTable, nResidues);
	// loop over all residues in current chain:
	for (ir = 0; ir < nResidues; ir++) {
		if (residueTable[ir]) {
			residueTable[ir]->GetAtomTable(atomTable, nAtoms);
			for (ia = 0; ia < nAtoms; ia++) {
				// atom exists
				if (atomTable[ia]) {
					if (strcmp((const char*) atomTable[ia]->name, " CA ")
							== 0) {
						matrix->m[iCA] = atomTable[ia]->x;
						matrix->m[iCA + 1] = atomTable[ia]->y;
						matrix->m[iCA + 2] = atomTable[ia]->z;
						/*printf("%f %f %f\n", atomTable[ia]->x,
						 atomTable[ia]->y, atomTable[ia]->z);*/
						iCA += 3;
					}
				}
			}
		}
	}

	return matrix;
}

std::unique_ptr<CMMDBManager> MMDBAndCarbonAlphaMatrix::toMMDB() {
	int RC, iResidue, modelId, bondReturn1, bondReturn2;
	cpstr chainId;
	CAtom *atom, *atomLast;
	std::unique_ptr<CMMDBManager> MMDB = std::make_unique<CMMDBManager>();
	//printf("Importing template PDB file %s...\n", "template.pdb");
	//RC = MMDB->ReadPDBASCII("template.pdb");
	RC = 0;
	MMDB->PutPDBString("COMPND    UNNAMED");
	MMDB->PutPDBString("AUTHOR    GENERATED BY PROTEIN KNOT DETECTOR 1.00");
	if (RC) {
		printf("Error #%i: %s\n", RC, GetErrorDescription(RC));
	} else {
		printf("Import Template SUCCESS\n");
		CModel *model = new CModel();
		MMDB->AddModel(model);
		modelId = model->GetSerNum();
		CChain *chain = new CChain();
		chainId = "A";
		chain->SetChainID(chainId);
		model->AddChain(chain);
		atomLast = nullptr;
		iResidue = 1; // Count Residues
		for (size_t i = 0; i < matrix_->n; i += 3) {
			/* the residue IS NOT associated with MMDB */
			CResidue *residue = new CResidue();
			residue->SetResID("ALA", iResidue, "");
			atom = new CAtom();
			atom->SetAtomName(" CA ");  // it has to be a PDB name!
			atom->SetElementName("C");
			/* In the PDB file format,
			 * each atom is given not only X, Y, and Z Cartesian coordinates,
			 * but two additional values immediately following
			 * called occupancy and temperature value
			 */
			atom->SetCoordinates(matrix_->m[i], matrix_->m[i + 1],
					matrix_->m[i + 2], 1.0, 1.0);
			/*printf("%f %f %f\n", matrix_->m[i], matrix_->m[i + 1],
			 matrix_->m[i + 2]);*/
			RC = residue->AddAtom(atom);
			if (RC <= 0) {
				// this may happen only if you try to add the same atom twice.
				printf("Error adding atom #%I64d\n", i);
			}
			chain->AddResidue(residue);
			/*
			 * bond this atom to the last
			 */
			if (atomLast != nullptr) {
				/*
				 * returns a positive (>0) number of total bonds for the atom
				 * in the case of successful completion
				 */
				/*RC = atom->AddBond(atomLast, 1, 1);
				 if (RC) {
				 printf("SUCCESS. BONDS: %i\n", RC);
				 }*/
			}
			atomLast = atom;
			iResidue++;
		}
		MMDB->PDBCleanup(
		PDBCLEAN_SERIAL | PDBCLEAN_INDEX);
		// MMDB->FinishStructEdit();
		/*
		 * Now we add bonds. We must do this after establishing the coordinate hierarchy
		 * because pointers obtained from CMMDBManager must be used.
		 */
		int ir;
		int nResidues;
		CResidue **residueTable;
		atomLast = nullptr;
		//  get chain
		printf("Getting Chain using modelId#%d chainId#%s\n", modelId,
				(char *) chainId);
		chain = MMDB->GetChain(modelId, chainId);
		// get residue table for current chain:
		chain->GetResidueTable(residueTable, nResidues);
		// loop over all residues in current chain:
		for (ir = 0; ir < nResidues; ir++) {
			if (residueTable[ir]) {
				atom = residueTable[ir]->GetAtom(0);
				if (atom) {
					if (atomLast != nullptr) {
						bondReturn1 = atom->AddBond(atomLast, 1);
						bondReturn2 = atomLast->AddBond(atom, 1);
						if (bondReturn1 && bondReturn2) {
							printf("OK #%d %f\n", ir, atom->x);
						} else {
							printf("ERROR #%d %f\n", ir, atom->x);
						}
					}
					atomLast = atom;
				}
			}
		}
	}
	return std::move(MMDB);
}

int OCCTManager::writeSTEP(char* path) {
	//CASCADE:Geom_Line = STEP:line
	STEPControl_Writer writer;
	if (!Interface_Static::SetIVal("write.precision.mode", 1)) {
		//error
	}
	if (!Interface_Static::SetIVal("write.step.assembly", 0)) {
		//error
	}
	if (!Interface_Static::SetCVal("write.step.schema", "DIS")) {
		//error
	}
	//STEP214Control_StepModelType mode = STEP214Control_ManifoldSolidBrep;
	gp_Pnt aPnt1(0, 0, 0);
	gp_Pnt aPnt2(10, 5, 2);
	gp_Pnt aPnt3(20, 3, 8);
	TopoDS_Edge aEdge1 = BRepBuilderAPI_MakeEdge(aPnt1, aPnt2);
	TopoDS_Edge aEdge2 = BRepBuilderAPI_MakeEdge(aPnt2, aPnt3);

	// Building the Resulting Compound
	TopoDS_Compound shape;
	BRep_Builder aBuilder;
	aBuilder.MakeCompound(shape);
	aBuilder.Add(shape, aEdge1);
	aBuilder.Add(shape, aEdge2);

	// Write file
	STEPControl_StepModelType mode = STEPControl_AsIs ;
	IFSelect_ReturnStatus stat1 = writer.Transfer(shape, mode);
	IFSelect_ReturnStatus stat2 = writer.Write("filename.stp");
	return 0;
}

std::unique_ptr<DoubleMatrix> TaylorKnotAlgorithm::getMatrix() {
	return std::move(m);
}
void TaylorKnotAlgorithm::setMatrix(std::unique_ptr<DoubleMatrix> matrixPtr) {
	m = std::move(matrixPtr);
}
void TaylorKnotAlgorithm::smooth() {
	int n = m->n - 3;
	for (int i = 3; i < n; i++) {
		m->m[i] = ((m->m[i - 3] + m->m[i + 3]) / 2 + m->m[i]) / 2;
	}
}
void TaylorKnotAlgorithm::smoothAuto() {

}

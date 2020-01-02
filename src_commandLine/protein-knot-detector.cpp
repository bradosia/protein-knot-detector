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

// proteinKnotDetector
#include "proteinKnotDetector/amalgamated.h"

using namespace std;

int main(int argc, char **argv) {
	int RC, im, ir;
	int nModels, nChains;
	CModel **modelTable;
	CChain **chainTable;
	unique_ptr<CMMDBManager> MMDB;
	unique_ptr<CMMDBManager> MMDBExport;
	std::unique_ptr<DoubleMatrix> carbonAlphaMatrix;

	MMDB = std::make_unique<CMMDBManager>();

	bool outputEachIteration = CommandLineOptions::output_each_iteration(argc,
			argv).value_or(false);
	std::string outputType =
			CommandLineOptions::output_type(argc, argv).value_or("pdb");
	std::string inputType = CommandLineOptions::input_type(argc, argv).value_or(
			"pdb");
	std::string inputFile = CommandLineOptions::input_file(argc, argv).value_or(
			"test.pdb");
	inputType = "pdb";
	inputFile = "test.pdb";

	MMDB->SetFlag(
			MMDBF_PrintCIFWarnings | MMDBF_FixSpaceGroup
					| MMDBF_IgnoreDuplSeqNum | MMDBF_IgnoreHash);

	if (inputType == "pdb") {
		printf("Reading PDB file %s\n", inputFile.c_str());
		RC = MMDB->ReadPDBASCII(inputFile.c_str());
	} else if (inputType == "cif") {
		printf("Reading CIF file %s\n", inputFile.c_str());
		RC = MMDB->ReadCIFASCII(inputFile.c_str());
	} else if (inputType == "mmdb") {
		printf("Reading MMDB binary file %s\n", inputFile.c_str());
		RC = MMDB->ReadMMDBF(inputFile.c_str());
	} else if (inputType == "crd") {
		printf("Reading Coordinate file %s\n", inputFile.c_str());

	}

//    4.3 Check for possible errors:
	if (RC) {
		//  An error was encountered. MMDB provides an error messenger
		//  function for easy error message printing.
		printf(" ***** ERROR #%i READ:\n\n %s\n\n", RC,
				GetErrorDescription(RC));
	} else {
		printf("Read file OK\n");

		int atomTotalNumber = MMDB->GetNumberOfAtoms();
		int modelTotalNumber = MMDB->GetNumberOfModels();

		printf("Total Atoms: %d\n"
				"Total Models: %d\n", atomTotalNumber, modelTotalNumber);

		/*
		 * select the first model
		 */
		printf("Selecting First Chain of Model...\n");
		int modelId = -1;
		cpstr chainId;
		MMDB->GetModelTable(modelTable, nModels);
		for (im = 0; im < nModels; im++) {
			if (modelTable[im]) {
				modelId = modelTable[im]->GetSerNum();
				modelTable[im]->GetChainTable(chainTable, nChains);
				for (ir = 0; ir < nChains; ir++) {
					chainId = chainTable[ir]->GetChainID();
					break;
				}
			}
		}

		if (chainId) {
			printf("Using Model SerNum#%d ChainId#%s\n", modelId, chainId);
			MMDBAndCarbonAlphaMatrix converter;
			printf("Setting Converter...\n");
			converter.setMMDBModel(std::move(MMDB), modelId, chainId);
			printf("Generating Alpha Carbon Matrix...\n");
			carbonAlphaMatrix = converter.toMatrix();
			//printf("Alpha Carbon Matrix:\n");
			//carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix1.txt");
			printf("Convert matrix to MMDB:\n");
			MMDBAndCarbonAlphaMatrix converter1;
			converter1.setMatrix(std::move(carbonAlphaMatrix));
			MMDBExport = converter1.toMMDB();
			carbonAlphaMatrix = converter1.getMatrix();
			printf("Exporting PDB\n");
			OCCTManager OCCT;
			OCCT.writeSTEP("out1.stp");
			//RC = MMDBExport->WritePDBASCII("out1.pdb");
			//RC = MMDBExport->WriteCIFASCII("out1.cif");
			//RC = MMDBExport->WriteMMDBF("out1.bin");
			printf("Running Taylor Knot Algorithm...\n");
			TaylorKnotAlgorithm taylorAlgorithm;
			taylorAlgorithm.setMatrix(std::move(carbonAlphaMatrix));
			taylorAlgorithm.smooth();
			carbonAlphaMatrix = taylorAlgorithm.getMatrix();
			//printf("Alpha Carbon Matrix:\n");
			//carbonAlphaMatrix->printMatrix();
			//carbonAlphaMatrix->writetoFileMatrix("matrix1.txt");
		}

	}

	system("pause");
	return 0;
}

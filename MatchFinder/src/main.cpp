// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>

#include "keyPoint.h"
#include "imageHelper.h"
#include "scannedScene.h"
#include "globalAppState.h"

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(1042);

	try {
		const std::string fileNameDescGlobalApp = "zParametersDefault.txt";
		GAS::loadGlobalAppState(fileNameDescGlobalApp);
		std::cout << std::endl;
		 
		const std::string srcPath = GAS::get().s_srcPath;
		const std::string outPath = GAS::get().s_outPath;
		Directory rootDir(srcPath);
		std::cout << "found " << rootDir.getDirectories().size() << " scenes " << std::endl;

		for (size_t dirIdx = 0; dirIdx < rootDir.getDirectories().size(); dirIdx++) {
			const std::string& s = rootDir.getDirectories()[dirIdx];
			if (s == "archive") continue;

			std::cout << "Loading Scene: " << s << std::endl;
			const std::string path = srcPath + "/" + s;
			 
			ScannedScene ss(path, s);			

			ss.findKeyPoints();
			ss.matchKeyPoints();
			ss.negativeKeyPoints();

			bool useTorchOutput = true;
			std::cout << "writing out matches to " << s << "_matches.txt | _negativ.txt" << std::endl;
			ss.saveMatches(outPath + "/" + s + "/" + "matches.txt", ss.getMatches(), useTorchOutput);
			ss.saveMatches(outPath + "/" + s + "/" + "negatives.txt", ss.getNegatives(), useTorchOutput);

			//const size_t numPairs = 10;
			//const size_t minMatches = 5;
			//ss.visulizeMatches(numPairs, minMatches);

			ss.saveImages(outPath + "/" + s + "/images/");

			if (GAS::get().s_maxNumScenes > 0 && dirIdx + 1 >= GAS::get().s_maxNumScenes) break;
		}

	}
	catch (const std::exception& e)
	{
		MessageBoxA(NULL, e.what(), "Exception caught", MB_ICONERROR);
		exit(EXIT_FAILURE);
	}
	catch (...)
	{
		MessageBoxA(NULL, "UNKNOWN EXCEPTION", "Exception caught", MB_ICONERROR);
		exit(EXIT_FAILURE);
	}
	
	std::cout << "<press key to continue>" << std::endl;
	getchar();
	return 0;
}


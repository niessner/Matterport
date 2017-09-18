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


	try {
		const std::string fileNameDescGlobalApp = "zParametersDefault.txt";
		GAS::loadGlobalAppState(fileNameDescGlobalApp);
		std::cout << std::endl;

		const std::string srcPath = GAS::get().s_srcPath;
		const std::string outPath = GAS::get().s_outPath;
		Directory rootDir(srcPath);
		std::vector<std::string> scenes;
		//scenes = rootDir.getDirectories();
		//scenes = { "zsNo4HB9uLZ" };

		std::vector<std::string> fileLists = {
			"../../scanner-ipad/Tasks/benchmark/checked_train.txt",
			"../../scanner-ipad/Tasks/benchmark/checked_test.txt" };
		for (const auto filelist : fileLists) {
			std::ifstream s(filelist); std::string line;
			while (std::getline(s, line)) scenes.push_back(line);
		}

		std::cout << "found " << scenes.size() << " scenes " << std::endl;

		//for (size_t dirIdx = 0; dirIdx < 50; dirIdx++) {
		for (size_t dirIdx = 0; dirIdx < scenes.size(); dirIdx++) {
			const std::string& s = scenes[dirIdx];
			if (s == "archive") continue;

			if (util::fileExists(outPath + "/" + s + "/matches1.txt")) {
			//if (util::directoryExists(outPath + "/" + s)) {
				std::cout << (outPath + "/" + s) << " already exists, skippping" << std::endl;
				continue;
			}

			std::cout << "Loading Scene: " << s << std::endl;
			const std::string path = srcPath + "/" + s;

			//ScannedScene::debug();

			ScannedScene ss(path, s);
			//ss.debugMatch();
			//getchar();

			//ss.computeNormals(ScannedScene::MESH_NORMALS);
			//ss.saveNormalImages(outPath + "/" + s + "/normals_mesh/");

			//if (!util::directoryExists(outPath + "/" + s + "/normals_depth/")) {
			//	ss.computeNormals(ScannedScene::DEPTH_NORMALS);
			//	ss.saveNormalImages(outPath + "/" + s + "/normals_depth/");
			//}
			//else {
			//	std::cout << "found depth normals [" << s << "]" << std::endl;
			//}

			ss.findKeyPoints();
			ss.matchKeyPoints();
			ss.negativeKeyPoints();

			if (!ss.getMatches().empty()) {
				bool useTorchOutput = true;
				std::cout << "writing out matches to " << outPath + "/" + s << "_matches.txt | _negativ.txt" << std::endl;

				if (!util::directoryExists(outPath + "/" + s)) util::makeDirectory(outPath + "/" + s);
				ss.saveMatches(outPath + "/" + s + "/" + "matches1.txt", ss.getMatches(), useTorchOutput);
				ss.saveMatches(outPath + "/" + s + "/" + "negatives1.txt", ss.getNegatives(), useTorchOutput);

				//const size_t numPairs = 10;
				//const size_t minMatches = 5;
				//ss.visulizeMatches(numPairs, minMatches);

				ss.saveImages(outPath + "/" + s + "/images/");
			}
			else {
				std::cout << "warning: no matches found for scene " << s << std::endl;
			}
			if (GAS::get().s_maxNumScenes > 0 && dirIdx + 1 >= GAS::get().s_maxNumScenes) break;

			std::cout << std::endl;
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


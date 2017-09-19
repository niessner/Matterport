// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>

#include "keyPoint.h"
#include "imageHelper.h"
#include "images.h"
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

		//check valid feature type
		std::unordered_set<std::string> validFeatureTypes = { "SIFT", "ORB", "SURF" };
		GAS::get().s_featureType = util::toUpper(GAS::get().s_featureType);
		if (validFeatureTypes.find(GAS::get().s_featureType) == validFeatureTypes.end())
			throw MLIB_EXCEPTION("invalid feature type: " + GAS::get().s_featureType);
		 
		std::string dataPath = GAS::get().s_dataPath;
		if (!util::directoryExists(dataPath)) throw MLIB_EXCEPTION("data path (" + dataPath + ") does not exist");
		if (!(dataPath.back() == '/' || dataPath.back() == '\\'))
			dataPath.push_back('/');
		const std::string sceneFileList = GAS::get().s_sceneFileList;
		if (!util::fileExists(sceneFileList)) throw MLIB_EXCEPTION("scene file list (" + sceneFileList + ") does not exist");
		std::vector<std::string> scenes;
		{
			std::ifstream s(sceneFileList);
			if (!s.is_open()) throw MLIB_EXCEPTION("failed to open " + sceneFileList + " for read");
			std::string line;
			while (std::getline(s, line))
				scenes.push_back(line);
		}
		////debugging
		//scenes = std::vector<std::string> {"17DRP5sb8fy"};
		//scenes.resize(2);
		////debugging
		std::cout << "found " << scenes.size() << " scenes " << std::endl;

		const std::string logFile = GAS::get().s_outputFile;
		if (util::fileExists(logFile)) {
			std::cout << "warning: log file " << logFile << " already exists, press key to delete and continue" << std::endl;
			getchar();
			util::deleteFile(logFile);
		}
		//const std::string statsLogFile = "stats.csv";
		//if (util::fileExists(statsLogFile)) {
		//	std::cout << "warning: stats log file " << statsLogFile << " already exists, press key to delete and continue" << std::endl;
		//	getchar();
		//	util::deleteFile(statsLogFile);
		//}

		Timer t;
		const bool bDumpTemp = scenes.size() > 1;
		const std::string tmpDir = "tmp/";
		if (util::directoryExists(tmpDir)) {
			std::cout << "warning: tmp dir " << tmpDir << " already exists, press key to delete and continue" << std::endl;
			getchar();
			util::deleteDirectory(tmpDir);
		}
		for (size_t dirIdx = 0; dirIdx < scenes.size(); dirIdx++) {
			const std::string& s = scenes[dirIdx];
			if (s == "archive") continue;

			const std::string path = dataPath + s;
			const std::string filenamePos = path + "/matches1.txt";
			const std::string filenameNeg = path + "/negatives1.txt";
			if (!util::fileExists(filenamePos) || !util::fileExists(filenameNeg)) continue;

			std::cout << "Loading Scene images: " << s << std::endl;
			 
			Images images(path, s);	

			images.loadGTMatches(filenamePos, filenameNeg);

			//debugging
			//images.computeKeyMatchStatistics(statsLogFile);
			//std::cout << "DONE" << std::endl;
			//getchar();

			images.matchKeyPoints(logFile, bDumpTemp, tmpDir);
		}
		if (bDumpTemp) { //aggregate
			Images::computePrecisionRecall(logFile, tmpDir);
		}
		std::cout << "total time: " << t.getElapsedTime() << " seconds" << std::endl;
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


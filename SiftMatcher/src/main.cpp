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

	if (false) { //count #train matches
		std::vector<std::string> scenes;
		{
			std::ifstream s("//tirion/share/datasets/Matterport/Matching/scenes_trainval.txt");
			MLIB_ASSERT(s.is_open());
			std::string line;
			while (std::getline(s, line))
				scenes.push_back(line);
		}
		std::vector<size_t> numMatchesPerScene; const unsigned int padding = 128 / 2; const unsigned int matchSkip = 4;
		std::vector<std::vector<vec2f>> keyposPos, keyposNeg;
		for (const auto& scene : scenes) {
			const std::string matchesFile = "//tirion/share/datasets/Matterport/Matching/" + scene + "/matches.txt";
			const std::string negativFile = "//tirion/share/datasets/Matterport/Matching/" + scene + "/negatives.txt";
			if (util::fileExists(matchesFile)) {
				keyposPos.push_back(std::vector<vec2f>());
				std::ifstream s(matchesFile); std::string line; size_t lineCount = 0;
				while (std::getline(s, line)) {
					if (lineCount >= 3 && ((lineCount - 3) % matchSkip == 0 || (lineCount - 3) % matchSkip == 1)) {
						const auto parts = util::split(line, '\t');
						const auto partsPixelPos = util::split(parts[3], ' ');
						vec2f pixPos(util::convertTo<float>(partsPixelPos[0]), util::convertTo<float>(partsPixelPos[1]));
						pixPos *= 0.5f;
						keyposPos.back().push_back(pixPos);
					}
					lineCount++;
				}
			}
			if (util::fileExists(negativFile)) {
				keyposNeg.push_back(std::vector<vec2f>());
				std::ifstream s(negativFile); std::string line; size_t lineCount = 0;
				while (std::getline(s, line)) {
					if (lineCount >= 3 && ((lineCount - 3) % matchSkip == 0 || (lineCount - 3) % matchSkip == 1)) {
						const auto parts = util::split(line, '\t');
						const auto partsPixelPos = util::split(parts[3], ' ');
						vec2f pixPos(util::convertTo<float>(partsPixelPos[0]), util::convertTo<float>(partsPixelPos[1]));
						pixPos *= 0.5f;
						keyposNeg.back().push_back(pixPos);
					}
					lineCount++;
				}
			}
		}
		for (unsigned int p = 0; p < keyposPos.size(); p++) {
			const auto& keysPos = keyposPos[p];
			const auto& keysNeg = keyposNeg[p];
			size_t count = 0;
			for (unsigned int i = 0; i < keysPos.size(); i += 2) {
				const vec2f pixPosA = keysPos[i];
				const vec2f pixPosB = keysPos[i + 1];
				const vec2f pixNegA = keysNeg[i];
				const vec2f pixNegB = keysNeg[i + 1];
				if (pixPosA.x >= padding && pixPosA.y >= padding && pixPosA.x + padding <= 640 && pixPosA.y + padding <= 480 &&
					pixPosB.x >= padding && pixPosB.y >= padding && pixPosB.x + padding <= 640 && pixPosB.y + padding <= 480 &&
					pixNegA.x >= padding && pixNegA.y >= padding && pixNegA.x + padding <= 640 && pixNegA.y + padding <= 480 &&
					pixNegB.x >= padding && pixNegB.y >= padding && pixNegB.x + padding <= 640 && pixNegB.y + padding <= 480) {
					count++;
				}
			}
			numMatchesPerScene.push_back(count);
		}
		std::sort(numMatchesPerScene.begin(), numMatchesPerScene.end());
		size_t total = std::accumulate(numMatchesPerScene.begin(), numMatchesPerScene.end(), 0);
		float avg = (float)total / (float)numMatchesPerScene.size();
		std::cout << numMatchesPerScene << std::endl;
		std::cout << "total: " << total << std::endl;
		std::cout << "avg: " << avg << std::endl;
		getchar();
		return 0;
	}
	 
	try {
		const std::string fileNameDescGlobalApp = "zParametersDefault.txt";
		GAS::loadGlobalAppState(fileNameDescGlobalApp);
		std::cout << std::endl;
		 
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
		//debugging
		scenes = std::vector<std::string> {"17DRP5sb8fy"};
		//debugging
		std::cout << "found " << scenes.size() << " scenes " << std::endl;

		const std::string logFile = "test.csv";
		if (util::fileExists(logFile)) {
			std::cout << "warning: log file " << logFile << " already exists, press key to delete and continue" << std::endl;
			getchar();
			util::deleteFile(logFile);
		}

		for (size_t dirIdx = 0; dirIdx < scenes.size(); dirIdx++) {
			const std::string& s = scenes[dirIdx];
			if (s == "archive") continue;
			  
			std::cout << "Loading Scene images: " << s << std::endl;
			const std::string path = dataPath + s;
			 
			Images images(path, s);	

			const std::string filenamePos = path + "/matches.txt";
			const std::string filenameNeg = path + "/negatives.txt";
			images.loadGTMatches(filenamePos, filenameNeg);
			images.matchKeyPoints(logFile);
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


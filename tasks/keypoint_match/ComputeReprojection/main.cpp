// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "scannedScene.h"
#include <iomanip>

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(7545);
	try {
		std::srand(0);

		const float maxDepth = 6.0f;
		const size_t numFramePairsToSample = 100;
#ifndef FIND_FAR_CAMERAS
		const std::string dataPath = "W:/data/matterport/v1_converted/";
		const std::vector<std::string> sceneFileLists = {
			"//tirion/share/datasets/Matterport/Matching1/scenes_trainval.txt",
			"//tirion/share/datasets/Matterport/Matching1/scenes_test.txt"
		};
		const size_t maxNumSampleTries = 10000;
#else
		const std::string dataPath = "W:/data/scan-net/scans/checked/";
		const std::vector<std::string> sceneFileLists = {
			"E:/Work/scanner-ipad/Tasks/benchmark/checked_train-0.txt",
			"E:/Work/scanner-ipad/Tasks/benchmark/checked_train-1.txt",
			"E:/Work/scanner-ipad/Tasks/benchmark/checked_val.txt",
			"E:/Work/scanner-ipad/Tasks/benchmark/checked_test.txt"
		};
		const size_t maxNumSampleTries = 100000;
#endif
		const std::string logFile = "log.csv";
		unsigned int maxNumScenes = 0;

		if (!util::directoryExists(dataPath)) throw MLIB_EXCEPTION("data path (" + dataPath + ") does not exist!");
		if (util::fileExists(logFile)) {
			std::cout << "warning: log file (" << logFile << ") already exists, press key to delete and continue" << std::endl;
			getchar();
			util::deleteFile(logFile);
		}
		std::vector<std::string> scenes; //read scenes
		for (const auto& file : sceneFileLists) {
			std::ifstream s(file); std::string line;
			if (!s.is_open()) throw MLIB_EXCEPTION("scene file list (" + file + ") does not exist!");
			while (std::getline(s, line)) 
				scenes.push_back(line);
			s.close();
		}
		//{
		//	std::vector<std::string> tmp;
		//	for (unsigned int i = 0; i < 380; i++) tmp.push_back(scenes[i]);
		//	scenes = tmp;
		//}
		//scenes = { "ScanNet-2016-07-24_01-55-26" };//debugging
		
		if (maxNumScenes == 0 || scenes.size() < maxNumScenes) maxNumScenes = (unsigned int)scenes.size();
		std::cout << maxNumScenes << " scenes" << std::endl;
		std::vector<ReprojError> errors(maxNumScenes);
		std::vector< std::vector<float> > sampledCameraDists(maxNumScenes);

		ReprojError total;
		std::ofstream s(logFile); const std::string splitter = ",";
		s << "scene" << splitter << "#corrs" << splitter << "depth l1" << splitter << "depth l2" << splitter << "intensity l1" << splitter << "intensity grad l1" << std::endl;
		for (unsigned int i = 0; i < maxNumScenes; i++) {
			if (!util::directoryExists(dataPath + scenes[i])) {
				std::cout << "WARNING: " << (dataPath + scenes[i]) << " does not exist, skipping" << std::endl;
				continue;
			}
			Timer t;

			//load sens
			ScannedScene scene(dataPath + scenes[i], scenes[i]);
			ReprojError e = scene.computeReprojection(numFramePairsToSample, maxNumSampleTries, maxDepth, sampledCameraDists[i]);
			errors[i] = e;
			if (e.numCorrs > 0) {
				s << scenes[i] << splitter << e.numCorrs << splitter << e.depthL1 << splitter << e.depthL2 << splitter << e.intensityL1 << splitter << e.intensityGradL1 << std::endl;
				total += e;
			}

			t.stop();
			std::cout << "processed [ " << (i + 1) << " | " << maxNumScenes << " ] -> time " << t.getElapsedTime() << " s" << std::endl;
		}
		s << "TOTAL" << splitter << total.numCorrs << splitter << total.depthL1 << splitter << total.depthL2 << splitter << total.intensityL1 << splitter << total.intensityGradL1 << std::endl;
		total.normalize();
		s << "TOTAL NORM" << splitter << splitter << total.depthL1 << splitter << total.depthL2 << splitter << total.intensityL1 << splitter << total.intensityGradL1 << std::endl;
		s.close();
		{ //cameras
			std::ofstream ofs("cameras.csv");
			ofs << "scene,camera distance" << std::endl;
			for (unsigned int i = 0; i < maxNumScenes; i++)
				for (unsigned int j = 0; j < sampledCameraDists[i].size(); j++)
					ofs << scenes[i] << "," << sampledCameraDists[i][j] << std::endl;
			ofs.close();
		}
		
		std::cout << "TOTAL NORM:" << std::endl;
		std::cout << "orig #corrs = " << total.numCorrs << std::endl;
		std::cout << "depth l1 = " << total.depthL1 << std::endl;
		std::cout << "depth l2 = " << total.depthL2 << std::endl;
		std::cout << "intensity l1 = " << total.intensityL1 << std::endl;
		std::cout << "intensity grad l1 = " << total.intensityGradL1 << std::endl;

		std::cout << "DONE" << std::endl;
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


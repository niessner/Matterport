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

		const size_t numFramePairsToSample = 100;
		const size_t maxNumSampleTries = 10000;
		const std::string dataPath = "W:/data/matterport/v1_converted/";
		//const std::string dataPath = "W:/data/scan-net/scans/checked/";
		const std::string logFile = "log.csv";
		unsigned int maxNumScenes = 1;

		if (!util::directoryExists(dataPath)) throw MLIB_EXCEPTION("data path (" + dataPath + ") does not exist!");

		Directory srcDir(dataPath);
		const auto scenes = srcDir.getDirectories();
		if (maxNumScenes != 0 && scenes.size() < maxNumScenes) maxNumScenes = (unsigned int)scenes.size();
		std::vector<ReprojError> errors(maxNumScenes);

//#pragma omp parallel for
		for (int i = 0; i < (int)maxNumScenes; i++) {
			if (scenes[i] == "archive") continue;
			Timer t;

			//load sens
			ScannedScene scene(srcDir.getPath() + scenes[i], scenes[i]);
			errors[i] = scene.computeReprojection(numFramePairsToSample, maxNumSampleTries);

			t.stop();
			std::cout << "processed [ " << (i + 1) << " | " << maxNumScenes << " ] -> time " << t.getElapsedTime() << " s" << std::endl;
		}
		ReprojError total;
		std::ofstream s(logFile); const std::string splitter = ",";
		s << "scene" << splitter << "#corrs" << splitter << "depth l1" << splitter << "depth l2" << splitter << "intensity l1" << splitter << "intensity grad l1" << std::endl;
		for (unsigned int i = 0; i < maxNumScenes; i++) {
			const auto& e = errors[i];
			if (e.numCorrs > 0) {
				s << scenes[i] << splitter << e.numCorrs << splitter << e.depthL1 << splitter << e.depthL2 << splitter << e.intensityL1 << splitter << e.intensityGradL1 << std::endl;
				total += e;
			}
		}
		s << "TOTAL" << splitter << total.numCorrs << splitter << total.depthL1 << splitter << total.depthL2 << splitter << total.intensityL1 << splitter << total.intensityGradL1 << std::endl;
		total.normalize();
		s << "TOTAL NORM" << splitter << splitter << total.depthL1 << splitter << total.depthL2 << splitter << total.intensityL1 << splitter << total.intensityGradL1 << std::endl;
		s.close();
		
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


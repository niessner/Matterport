// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>

#include "keyPoint.h"
#include "imageHelper.h"
#include "scannedScene.h"

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(1042);
	try {
		const std::string srcPath = "W:/data/matterport/v1_converted";

		Directory rootDir(srcPath);

		for (const std::string& s : rootDir.getDirectories()) {
			if (s == "archive") continue;

			std::cout << s << std::endl;
			const std::string path = srcPath + "/" + s;
			 
			ScannedScene ss(path, s);
			ss.findKeyPoints();
			ss.matchKeyPoints();
			
			break;
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


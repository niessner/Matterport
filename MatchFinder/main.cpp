// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>


class KeyPoint {
public:
	unsigned int m_sensorIdx;
	unsigned int m_frameIdx;
	vec2f m_pixelPos;
	float m_depth;
	vec3f m_worldPos;
};

class KeyPointMatch {
public:
	KeyPoint m_kp0;
	KeyPoint m_kp1;

	//std::vector<KeyPoint> m_keyPoints;
};


class ScannedScene {
public:
	ScannedScene(const std::string& path, const std::string& name) {
		load(path, name);
	}
	~ScannedScene() {
		for (auto* sd : m_sds) {
			SAFE_DELETE(sd);
		}
	}

	void load(const std::string& path, const std::string& name) {
		
		m_name = name;

		Directory dir(path);
		const auto& files = dir.getFilesWithSuffix(".sens");
		for (auto& f : files) {
			m_sds.push_back(new SensorData);
			SensorData* sd = m_sds.back();
			sd->loadFromFile(path + "/" + f);
			std::cout << *sd << std::endl;
		}
	}

	void findKeyPoints() {
		for (auto* sd : m_sds) {
			for (size_t i = 0; i < sd->m_frames.size(); i++) {
				ColorImageR8G8B8 c = sd->computeColorImage(i);
				DepthImage16 d = sd->computeDepthImage(i);

				//TODO find key points
			}
		}
	}
private:
	std::list<SensorData*> m_sds;
	std::string m_name;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_keyPointMatches;

};

int main(int argc, char* argv[])
{
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
	//_CrtSetBreakAlloc(7545);
	try {

		const std::string srcPath = "W:/data/matterport/v1_converted";

		Directory rootDir(srcPath);

		for (const std::string& s : rootDir.getDirectories()) {
			if (s == "archive") continue;

			std::cout << s << std::endl;
			const std::string path = srcPath + "/" + s;

			ScannedScene ss(path, s);
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


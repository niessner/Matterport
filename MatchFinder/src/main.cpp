// sensorFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iomanip>

#include "keyPoint.h"
#include "imageHelper.h"

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
		auto& files = dir.getFilesWithSuffix(".sens");
		std::sort(files.begin(), files.end());

		for (auto& f : files) {
			m_sds.push_back(new SensorData);
			SensorData* sd = m_sds.back();
			sd->loadFromFile(path + "/" + f);
			std::cout << *sd << std::endl;
		}
	}

	void findKeyPoints() {
		for (size_t sensorIdx = 0; sensorIdx < m_sds.size(); sensorIdx++) {
			SensorData* sd = m_sds[sensorIdx];
			const mat4f intrinsicInv = sd->m_calibrationDepth.m_intrinsic.getInverse();

			for (size_t imageIdx = 0; imageIdx < sd->m_frames.size(); imageIdx++) {
				ColorImageR8G8B8 c = sd->computeColorImage(imageIdx);
				DepthImage16 d = sd->computeDepthImage(imageIdx);

				const mat4f& camToWorld = sd->m_frames[imageIdx].getCameraToWorld();

				const unsigned int maxNumKeyPoints = 512;
				std::vector<vec3f> rawKeyPoints = KeyPointFinder::findKeyPoints(c, maxNumKeyPoints);

				for (vec3f& rawkp : rawKeyPoints) {
					vec2ui loc = math::round(vec2f(rawkp.x, rawkp.y));
					if (d.isValid(loc)) {
						KeyPoint kp;
						kp.m_depth = d(loc) / sd->m_depthShift;
						kp.m_frameIdx = (unsigned int)imageIdx;
						kp.m_sensorIdx = (unsigned int)sensorIdx;
						//kp.m_pixelPos = vec2f(rawkp.x, rawkp.y);
						kp.m_pixelPos = vec2f(loc);

						vec3f cameraPos = (intrinsicInv*vec4f(kp.m_pixelPos.x*kp.m_depth, kp.m_pixelPos.y*kp.m_depth, kp.m_depth, 0.0f)).getVec3();
						kp.m_worldPos = camToWorld * cameraPos;
					}
				}
			}
		}
	}
private:
	std::vector<SensorData*> m_sds;
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
			ss.findKeyPoints();
			
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


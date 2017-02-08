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
	vec2f m_offset;	//re-projection offset when m_kp0 is projected into m_kp1;

	//std::vector<KeyPoint> m_keyPoints;
};


inline float gaussR(float sigma, float dist)
{
	return exp(-(dist*dist) / (2.0f*sigma*sigma));
}
inline float gaussR(float sigma, const vec3f& d)
{
	float dist = d.length();
	return exp(-(dist*dist) / (2.0f*sigma*sigma));
}
inline float gaussR(float sigma, const vec3uc& d)
{
	vec3f _d(d);	//_d /= 255.0f;
	float dist = _d.length();
	return exp(-(dist*dist) / (2.0f*sigma*sigma));
}
 

inline float linearR(float sigma, float dist)
{
	return std::max(1.0f, std::min(0.0f, 1.0f - (dist*dist) / (2.0f*sigma*sigma)));
}

inline float gaussD(float sigma, int x, int y)
{
	return exp(-((x*x + y*y) / (2.0f*sigma*sigma)));
}

inline float gaussD(float sigma, int x)
{
	return exp(-((x*x) / (2.0f*sigma*sigma)));
}

void bilateralFilter(BaseImage<vec3uc>& img, float sigmaD, float sigmaR) {
		
	BaseImage<vec3uc> res(img.getDimensions());
	res.setInvalidValue(img.getInvalidValue());
 
	const int kernelRadius = (int)ceil(2.0*sigmaD);
	for (unsigned int y = 0; y < img.getHeight(); y++) {
		for (unsigned int x = 0; x < img.getWidth(); x++) {
			
			res.setInvalid(x, y);

			vec3f sum = vec3f(0.0f);
			float sumWeight = 0.0f;
			
			if (img.isValid(x, y)) {
				const vec3uc& center = img(x, y);

				for (int m = x - kernelRadius; m <= (int)x + kernelRadius; m++) {
					for (int n = y - kernelRadius; n <= (int)y + kernelRadius; n++) {
						if (m >= 0 && n >= 0 && m < (int)img.getWidth() && n < (int)img.getHeight()) {
							if (img.isValid(m, n)) {
								const vec3uc& current = img(m, n);
								const float weight = gaussD(sigmaD, m - x, n - y)*gaussR(sigmaR, current - center);
								sumWeight += weight;
								sum += weight*vec3f(current);								
							} 
						}
					}
				}
				if (sumWeight > 0.0f) res(x, y) = math::round(sum / sumWeight);
			}
		}
	}
	img = res;
}

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

			break;
		}
	}

	void findKeyPoints() {
		for (size_t sensorIdx = 0; sensorIdx < m_sds.size(); sensorIdx++) {
			SensorData* sd = m_sds[sensorIdx];
			const mat4f intrinsicInv = sd->m_calibrationDepth.m_intrinsic.getInverse();

			for (size_t imageIdx = 0; imageIdx < sd->m_frames.size(); imageIdx++) {
				ColorImageR8G8B8 c = sd->computeColorImage(imageIdx);
				DepthImage16 d = sd->computeDepthImage(imageIdx);

				//float sigmaD = 2.0f;	
				//float sigmaR = 0.1f;
				//sigmaD = 10.0f;
				//sigmaR = 10.0f;
				//FreeImageWrapper::saveImage("_before.png", c);
				//bilateralFilter(c, sigmaD, sigmaR);
				//FreeImageWrapper::saveImage("_after.png", c);
				//std::cout << "here" << std::endl;
				//getchar();

				const mat4f& camToWorld = sd->m_frames[imageIdx].getCameraToWorld();

				const unsigned int maxNumKeyPoints = 512;
				std::vector<vec3f> rawKeyPoints = KeyPointFinder::findKeyPoints(c, maxNumKeyPoints);

				size_t validKeyPoints = 0;
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

						validKeyPoints++;

						m_keyPoints.push_back(kp);
					} 
				}
				std::cout << "\tfound " << validKeyPoints << " keypoints for image " << sensorIdx << "|" << imageIdx << std::endl;

				if (imageIdx == 5) break;
			} 
		}
	}

	
	void matchKeyPoints() {
		const float radius = 0.1f;	//10 cm
		const unsigned int maxK = 10;
		NearestNeighborSearchFLANNf nn(4*maxK, 1);
		 
		std::vector<float> points(3 * m_keyPoints.size());
		for (size_t i = 0; i < m_keyPoints.size(); i++) {
			points[3 * i + 0] = m_keyPoints[i].m_worldPos.x;
			points[3 * i + 1] = m_keyPoints[i].m_worldPos.y;
			points[3 * i + 2] = m_keyPoints[i].m_worldPos.z;
		}
		nn.init(points.data(), m_keyPoints.size(), 3, maxK);

		size_t numMatches = 0;
		for (size_t i = 0; i < m_keyPoints.size(); i++) {
			const float* query = &points[3*i];
			std::vector<unsigned int> res = nn.fixedRadius(query, maxK, radius);
			auto resPair = nn.fixedRadiusDist(query, maxK, radius);
			if (res.size() == 0 || res[0] != i) throw MLIB_EXCEPTION("should find itself...");
			for (size_t j = 1; j < res.size(); j++) {
				
				KeyPointMatch m;
				m.m_kp0 = m_keyPoints[i];
				m.m_kp1 = m_keyPoints[res[j]];

				mat4f worldToCam_kp1 = m_sds[m.m_kp1.m_sensorIdx]->m_frames[m.m_kp1.m_frameIdx].getCameraToWorld().getInverse();
				vec3f p = worldToCam_kp1 * m.m_kp0.m_worldPos;
				p = m_sds[m.m_kp1.m_sensorIdx]->m_calibrationDepth.cameraToProj(p);
				m.m_offset = m.m_kp1.m_pixelPos - vec2f(p.x, p.y);
				m_keyPointMatches.push_back(m);
				numMatches++;
			}
			std::cout << "numMatches: " << numMatches << std::endl;
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
	//_CrtSetBreakAlloc(1489);
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


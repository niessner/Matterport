
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "matchVisualization.h"
#include "globalAppState.h"

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

		//for (auto& f : files) {
		for (size_t i = 0; i < files.size(); i++) {
			auto& f = files[i];
			m_sds.push_back(new SensorData);
			SensorData* sd = m_sds.back();
			sd->loadFromFile(path + "/" + f);
			std::cout << *sd << std::endl;

			if (GAS::get().s_maxNumSensFiles > 0 && i + 1 >= GAS::get().s_maxNumSensFiles) break;
		}
	}

	//for a loaded scene, finds all key points in the images
	void findKeyPoints();

	void negativeKeyPoints();

	//matches all previously found key points between all images and loaded sens files
	void matchKeyPoints();

	void saveMatches(const std::string& filename, const std::vector<KeyPointMatch>& matches, bool torch = true) const {
		std::ofstream outFile(filename);

		if (!torch) {	//human readable one
			outFile << "SceneName " << m_name << " ( " << matches.size() << " matches )\n";
			outFile << "\n";
			for (size_t i = 0; i < matches.size(); i++) {
				outFile << "matchIdx " << i << "\n";
				outFile << matches[i] << "\n";
			}
		}
		else {
			outFile << "SceneName " << m_name << " ( " << matches.size() << " matches )\n";
			outFile << "\n";

			const std::string sep = "\t";
			outFile << "matchIdx" << sep << "m_sensorIdx" << sep << "m_imageIdx" << sep << "m_pixelPos" << sep << "m_depth" << sep << "m_worldPos" << sep << "m_offset" << sep
				<< "m_size" << sep << "m_angle" << sep << "m_octave" << sep << "m_scale" << sep << "m_opencvPackOctave" << "\n";
			for (size_t i = 0; i < matches.size(); i++) {
				const KeyPoint& k0 = matches[i].m_kp0;
				const KeyPoint& k1 = matches[i].m_kp1;
				outFile << i << sep << k0.m_sensorIdx << sep << k0.m_imageIdx << sep << k0.m_pixelPos << sep << k0.m_depth << sep << k0.m_worldPos << sep << vec2f(0.0f) << sep << k0.m_size << sep << k0.m_angle << sep << k0.m_octave << sep << k0.m_scale << sep << k0.m_opencvPackOctave << "\n";
				outFile << i << sep << k1.m_sensorIdx << sep << k1.m_imageIdx << sep << k1.m_pixelPos << sep << k1.m_depth << sep << k1.m_worldPos << sep << matches[i].m_offset << sep << k1.m_size << sep << k1.m_angle << sep << k1.m_octave << sep << k1.m_scale << sep << k1.m_opencvPackOctave << "\n";
			}
		}


	}

	void visulizeMatches(size_t numPairs = 10, size_t minMatches = 1) {
		MatchVisualization mv;
		mv.visulizeMatches(m_sds, m_keyPointMatches, numPairs, minMatches);
	}

	const std::vector<KeyPointMatch>& getMatches() const {
		return m_keyPointMatches;
	}
	const std::vector<KeyPointMatch>& getNegatives() const {
		return m_keyPointNegatives;
	}

	void saveImages(const std::string& outPath) const;

private:
	//static vec3f computeCameraSpacePosition(const mat4f& intrinsicsInv, const DepthImage32& depthImage, unsigned int x, unsigned int y)
	//{
	//	MLIB_ASSERT(x < depthImage.getWidth() && y < depthImage.getHeight());

	//	float depth = depthImage(x, y);
	//	if (depth != -std::numeric_limits<float>::infinity()) {
	//		const vec3f res = intrinsicsInv*vec3f(x*depth, y*depth, depth);
	//		return res;
	//	}
	//	return vec3f(-std::numeric_limits<float>::infinity());
	//}
	//static vec3f computeNormal(const mat4f& intrinsicsInv, const DepthImage32& depthImage, unsigned int x, unsigned int y)
	//{
	//	MLIB_ASSERT(x < depthImage.getWidth() && y < depthImage.getHeight());

	//	vec3f normal = vec3f(-std::numeric_limits<float>::infinity());

	//	if (x > 0 && x + 1 < depthImage.getWidth() && y > 0 && y + 1 < depthImage.getHeight()) {
	//		const vec3f CC = computeCameraSpacePosition(intrinsicsInv, depthImage, x + 0, y + 0); //d_input[(y + 0)*width + (x + 0)];
	//		const vec3f PC = computeCameraSpacePosition(intrinsicsInv, depthImage, x + 0, y + 1); //d_input[(y + 1)*width + (x + 0)];
	//		const vec3f CP = computeCameraSpacePosition(intrinsicsInv, depthImage, x + 1, y + 0); //d_input[(y + 0)*width + (x + 1)];
	//		const vec3f MC = computeCameraSpacePosition(intrinsicsInv, depthImage, x + 0, y - 1); //d_input[(y - 1)*width + (x + 0)];
	//		const vec3f CM = computeCameraSpacePosition(intrinsicsInv, depthImage, x - 1, y + 0); //d_input[(y + 0)*width + (x - 1)];

	//		if (CC.x != -std::numeric_limits<float>::infinity() && PC.x != -std::numeric_limits<float>::infinity() &&
	//			CP.x != -std::numeric_limits<float>::infinity() && MC.x != -std::numeric_limits<float>::infinity() &&
	//			CM.x != -std::numeric_limits<float>::infinity()) {
	//			const vec3f n = (PC - MC) ^ (CP - CM);
	//			const float l = n.length();
	//			if (l > 0.0f) normal = n / -l; //d_output[y*width + x] = make_float4(n / -l, 1.0f);
	//		}
	//	}
	//	return normal;
	//}

	std::vector<SensorData*> m_sds;
	std::string m_name;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_keyPointMatches;

	std::vector<KeyPointMatch>	m_keyPointNegatives;
};

#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "matchVisualization.h"
#include "globalAppState.h"

class ScannedScene {
public:
	enum NORMAL_TYPE {
		DEPTH_NORMALS,
		MESH_NORMALS
	};
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
		m_path = path;

		Directory dir(path);
		auto& files = dir.getFilesWithSuffix(".sens");
		std::sort(files.begin(), files.end());

		//scannet debugging
		if (files.size() > GAS::get().s_maxNumSensFiles)
			throw MLIB_EXCEPTION("error too many sens files(" + std::to_string(files.size()) + "): " + name);

		//for (auto& f : files) {
		for (size_t i = 0; i < files.size(); i++) {
			auto& f = files[i];
			m_sds.push_back(new SensorData);
			SensorData* sd = m_sds.back();
			sd->loadFromFile(path + "/" + f);
			std::cout << *sd << std::endl;

			if (GAS::get().s_maxNumSensFiles > 0 && i + 1 >= GAS::get().s_maxNumSensFiles) break;
		}

		//frames to consider
		m_frameIndices.clear();
		for (unsigned int sensIdx = 0; sensIdx < m_sds.size(); sensIdx++) {
			m_frameIndices.push_back(std::vector<unsigned int>(m_sds[sensIdx]->m_frames.size()));
			for (unsigned int frameIdx = 0; frameIdx < m_sds[sensIdx]->m_frames.size(); frameIdx++)
				m_frameIndices[sensIdx][frameIdx] = frameIdx;
		}
		const unsigned int maxNumFrames = GAS::get().s_maxNumFramesPerScene;
		if (maxNumFrames > 0) {
			const unsigned int maxNumFramesPerSens = maxNumFrames / (unsigned int)m_sds.size();
			for (auto& frameIndices : m_frameIndices) {
				if (frameIndices.size() > maxNumFramesPerSens) {
					std::random_shuffle(frameIndices.begin(), frameIndices.end()); frameIndices.resize(maxNumFramesPerSens);
					std::sort(frameIndices.begin(), frameIndices.end());
				}
			}
		}
	}

	//for a loaded scene, finds all key points in the images
	void findKeyPoints();

	void negativeKeyPoints();

	//matches all previously found key points between all images and loaded sens files
	void matchKeyPoints();

	void computeNormals(NORMAL_TYPE type);

	void saveMatches(const std::string& filename, const std::vector<KeyPointMatch>& matches, bool torch = true) const {
		if (util::fileExists(filename)) {
			std::cout << "warning: match file (" << filename << ") already exists, press key to continue" << std::endl;
			getchar();
		}
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
	void saveNormalImages(const std::string& outPath) const;

	static void debug();
	void debugMatch() const;

private:
	void computeDepthNormals();
	void computeMeshNormals();

	std::vector<SensorData*> m_sds;
	std::string m_name;
	std::string m_path;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_keyPointMatches;

	std::vector<KeyPointMatch>	m_keyPointNegatives;

	//normals
	std::vector<std::vector<PointImage>>	m_normals;

	//limit set
	std::vector<std::vector<unsigned int>> m_frameIndices;
};
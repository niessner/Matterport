
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "matchVisualization.h"


class ScannedScene {
public:
	ScannedScene(const std::string& path, const std::string& name, size_t maxNumSensFiles = 0, size_t maxNumImages = 0) {
		load(path, name, maxNumSensFiles, maxNumImages);
	}
	~ScannedScene() {
		for (auto* sd : m_sds) {
			SAFE_DELETE(sd);
		}
	} 

	void load(const std::string& path, const std::string& name, size_t maxNumSensFiles = 0, size_t maxNumImages = 0) {

		m_maxNumSensFiles = maxNumSensFiles;
		m_maxNumImages = maxNumImages;

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

			if (maxNumSensFiles > 0 && i + 1 >= maxNumSensFiles) break;
		}
	}

	//for a loaded scene, finds all key points in the images
	void findKeyPoints();
	
	//matches all previously found key points between all images and loaded sens files
	void matchKeyPoints();

	void saveMatches(const std::string& filename) {
		std::ofstream outFile(filename);

		outFile << "SceneName " << m_name << "\n";
		for (size_t i = 0; i < m_keyPointMatches.size(); i++) {
			outFile << m_keyPointMatches[i] << "\n";
		}		
	}

	void visulizeMatches(size_t numPairs = 10) {
		MatchVisualization mv;
		mv.visulizeMatches(m_sds, m_keyPointMatches, numPairs);
	}
private:
	std::vector<SensorData*> m_sds;
	std::string m_name;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_keyPointMatches;


	size_t m_maxNumSensFiles;
	size_t m_maxNumImages;
};
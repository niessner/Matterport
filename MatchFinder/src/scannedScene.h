
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

private:
	std::vector<SensorData*> m_sds;
	std::string m_name;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_keyPointMatches;

	std::vector<KeyPointMatch>	m_keyPointNegatives;
};
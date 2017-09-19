
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "matchVisualization.h"
#include "globalAppState.h"

class Images {
public:
	Images(const std::string& path, const std::string& name) {
		m_keyDetectImageWidth = GAS::get().s_keyDetectWidth;
		m_keyDetectImageHeight = GAS::get().s_keyDetectHeight;
		load(path, name);
	}
	~Images() {
	}

	void load(std::string path, const std::string& name) {
		if (!util::directoryExists(path)) throw MLIB_EXCEPTION("directory " + path + " does not exist");
		if (!(path.back() == '/' || path.back() == '\\')) path.push_back('/');

		m_name = name;

		std::cout << "loading images..." << std::endl;
		Directory dir(path + "images");
		auto& files = dir.getFilesWithSuffix(".jpg");
		//std::sort(files.begin(), files.end());

		const unsigned int maxNumImages = GAS::get().s_maxNumImages;
		m_images.resize(GAS::get().s_maxNumSensors);
		for (size_t i = 0; i < files.size(); i++) {
			const auto& f = files[i];
			const auto parts = util::split(util::removeExtensions(f), '-');
			if (parts.size() != 3)
				throw MLIB_EXCEPTION("invalid image file: " + f);
			const unsigned int sensIdx = util::convertTo<unsigned int>(parts[1]);
			const unsigned int imgIdx = util::convertTo<unsigned int>(parts[2]);
			if (maxNumImages == 0 || imgIdx < maxNumImages) {
				if (m_images[sensIdx].size() < imgIdx + 1) m_images[sensIdx].resize(imgIdx + 1);
				FreeImageWrapper::loadImage(dir.getPath() + f, m_images[sensIdx][imgIdx]);
			}
			std::cout << "\r[ " << (i + 1) << " | " << files.size() << " ]";
		}
		std::cout << std::endl << "done!" << std::endl;

		//KeyPointMatcher::debug(m_images); //debug the sift matching
	}

	//load key points from ground truth matches
	void loadGTMatches(const std::string& filenamePos, const std::string& filenameNeg);

	//matches all previously found key points between all images and loaded sens files
	void matchKeyPoints(const std::string& logFile, bool dumpTemporaryBinary, const std::string& dataPath = "") const;

	void computeKeyMatchStatistics(const std::string& logFile) const;

	//from tmp directory data
	static void computePrecisionRecall(const std::string& logFile, const std::string& dataPath = "");

private:

	void loadMatchFile(const std::string& filename, std::vector<KeyPointMatch>& matches,
		std::vector<KeyPoint>& keypoints, std::vector<std::unordered_map<KeyPoint, size_t>>& keymaps) const;
	static vec4ui evaluatePosMatch(const std::vector<float>& dists, const float thresh);
	static vec4ui evaluateNegMatch(const std::vector<float>& dists, const float thresh);
	static void computePrecisionRecall(const std::vector<float> &matchDistsPos, const std::vector<float> &matchDistsNeg, const std::string& logFile);

	std::vector<std::vector<ColorImageR8G8B8>> m_images;
	std::string m_name;

	std::vector<KeyPoint>		m_keyPoints;
	std::vector<KeyPointMatch>	m_gtKeyPointMatchesPos;
	std::vector<KeyPointMatch>	m_gtKeyPointMatchesNeg;
	
	unsigned int m_keyDetectImageWidth;
	unsigned int m_keyDetectImageHeight;
};
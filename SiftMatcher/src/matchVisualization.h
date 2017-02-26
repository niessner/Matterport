
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "imageHelper.h"

class MatchVisualization {
public:
	static void visulizeMatches(const std::vector<SensorData*>& sds, const std::vector<KeyPointMatch>& matches, 
		const std::vector<KeyPoint>& keypoints, size_t numPairs, size_t minMatches = 1) {
		
		size_t currMatch = 0;
		for (size_t i = 0; i < numPairs;) {

			std::vector<KeyPointMatch> curr;
			for (;; currMatch++) {
				if (currMatch >= matches.size()) break;

				if (curr.size() == 0)
					curr.push_back(matches[currMatch]);
				else if (curr.front().isSameImagePair(matches[currMatch], keypoints))
					curr.push_back(matches[currMatch]);
				else break;
			}
			if (curr.size() >= minMatches) {
				const KeyPointMatch& r = curr.front();
				ColorImageR8G8B8 img0 = sds[keypoints[r.m_kp0].m_sensorIdx]->computeColorImage(keypoints[r.m_kp0].m_imageIdx);
				ColorImageR8G8B8 img1 = sds[keypoints[r.m_kp1].m_sensorIdx]->computeColorImage(keypoints[r.m_kp1].m_imageIdx);
				ColorImageR8G8B8 m = match(img0, img1, curr, keypoints);
				const std::string filename = "matches_" + std::to_string(keypoints[r.m_kp0].m_imageIdx) + "_" + std::to_string(keypoints[r.m_kp1].m_imageIdx) + ".png";
				std::cout << "creating debug file: " << filename << " ( " << curr.size() << " matches ) " <<  std::endl;
				FreeImageWrapper::saveImage(filename, m);
				i++;
			}

			if (currMatch >= matches.size()) break;
		}
	}
	static void visulizeMatches(const std::vector<std::vector<ColorImageR8G8B8>>& images, const std::vector<KeyPointMatch>& matches,
		const std::vector<KeyPoint>& keypoints, size_t numPairs, size_t minMatches = 1) {

		size_t currMatch = 0;
		for (size_t i = 0; i < numPairs;) {

			std::vector<KeyPointMatch> curr;
			for (;; currMatch++) {
				if (currMatch >= matches.size()) break;

				if (curr.size() == 0)
					curr.push_back(matches[currMatch]);
				else if (curr.front().isSameImagePair(matches[currMatch], keypoints))
					curr.push_back(matches[currMatch]);
				else break;
			}
			if (curr.size() >= minMatches) {
				const KeyPointMatch& r = curr.front();
				const ColorImageR8G8B8& img0 = images[keypoints[r.m_kp0].m_sensorIdx][keypoints[r.m_kp0].m_imageIdx];
				const ColorImageR8G8B8& img1 = images[keypoints[r.m_kp1].m_sensorIdx][keypoints[r.m_kp1].m_imageIdx];
				ColorImageR8G8B8 m = match(img0, img1, curr, keypoints);
				const std::string filename = "matches_" + std::to_string(keypoints[r.m_kp0].m_imageIdx) + "_" + std::to_string(keypoints[r.m_kp1].m_imageIdx) + ".png";
				std::cout << "creating debug file: " << filename << " ( " << curr.size() << " matches ) " << std::endl;
				FreeImageWrapper::saveImage(filename, m);
				i++;
			}

			if (currMatch >= matches.size()) break;
		}
	}
private:

	static ColorImageR8G8B8 match(const ColorImageR8G8B8& img0, const ColorImageR8G8B8& img1, const std::vector<KeyPointMatch>& matches, const std::vector<KeyPoint>& keypoints) {

		ColorImageR8G8B8 matchImage(img0.getWidth() * 2, img0.getHeight());
		matchImage.copyIntoImage(img0, 0, 0);
		matchImage.copyIntoImage(img1, img1.getWidth(), 0);

		RGBColor lowColor = ml::RGBColor::Blue;
		RGBColor highColor = ml::RGBColor::Red;
		for (size_t i = 0; i < matches.size(); i++) {
			const KeyPointMatch& kp = matches[i];

			//RGBColor c = RGBColor::interpolate(lowColor, highColor, 0.5f);
			RGBColor c = RGBColor::randomColor();
			vec2i p0 = ml::math::round(ml::vec2f(keypoints[kp.m_kp0].m_pixelPos.x, keypoints[kp.m_kp0].m_pixelPos.y));
			vec2i p1 = ml::math::round(ml::vec2f(keypoints[kp.m_kp1].m_pixelPos.x, keypoints[kp.m_kp1].m_pixelPos.y));
			p1.x += img0.getWidth();

			float size0 = keypoints[kp.m_kp0].m_size;
			float size1 = keypoints[kp.m_kp1].m_size;
			if (size0 == -std::numeric_limits<float>::infinity()) {
				size0 = 3.0f;
				std::cout << "warning: no key size available, using default size " << size0 << std::endl;
			}
			if (size1 == -std::numeric_limits<float>::infinity()) {
				size1 = 3.0f;
				std::cout << "warning: no key size available, using default size " << size1 << std::endl;
			}
			ImageHelper::drawCircle(matchImage, p0, ml::math::round(size0), c.getVec3());
			ImageHelper::drawCircle(matchImage, p1, ml::math::round(size1), c.getVec3());
			ImageHelper::drawLine(matchImage, p0, p1, c.getVec3());
		}
		return matchImage;
	}

};
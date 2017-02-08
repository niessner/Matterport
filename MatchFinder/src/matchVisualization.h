
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "imageHelper.h"

class MatchVisualization {
public:
	void visulizeMatches(const std::vector<SensorData*>& sds, const std::vector<KeyPointMatch>& matches, size_t numPairs) {
		
		size_t currMatch = 0;
		for (size_t i = 0; i < numPairs;) {

			std::vector<KeyPointMatch> curr;
			for (;; currMatch++) {
				if (currMatch >= matches.size()) return;	//nothing else to do here...

				if (curr.size() == 0)
					curr.push_back(matches[currMatch]);
				else if (curr.front().isSameImagePair(matches[currMatch]))
					curr.push_back(matches[currMatch]);
				else break;
			}
			if (curr.size()) {
				const KeyPointMatch& r = curr.front();
				ColorImageR8G8B8 img0 = sds[r.m_kp0.m_sensorIdx]->computeColorImage(r.m_kp0.m_imageIdx);
				ColorImageR8G8B8 img1 = sds[r.m_kp1.m_sensorIdx]->computeColorImage(r.m_kp1.m_imageIdx);
				ColorImageR8G8B8 m = match(img0, img1, curr);
				const std::string filename = "matches_" + std::to_string(r.m_kp0.m_imageIdx) + "_" + std::to_string(r.m_kp1.m_imageIdx) + ".png";
				std::cout << "creating debug file: " << filename << " ( " << curr.size() << " matches ) " <<  std::endl;
				FreeImageWrapper::saveImage(filename, m);
				i++;
			}
		}
	}
private:

	ColorImageR8G8B8 match(ColorImageR8G8B8& img0, ColorImageR8G8B8& img1, const std::vector<KeyPointMatch>& matches) {

		ColorImageR8G8B8 matchImage(img0.getWidth() * 2, img0.getHeight());
		matchImage.copyIntoImage(img0, 0, 0);
		matchImage.copyIntoImage(img1, img1.getWidth(), 0);

		RGBColor lowColor = ml::RGBColor::Blue;
		RGBColor highColor = ml::RGBColor::Red;
		for (size_t i = 0; i < matches.size(); i++) {
			const KeyPointMatch& kp = matches[i];

			//RGBColor c = RGBColor::interpolate(lowColor, highColor, 0.5f);
			RGBColor c = RGBColor::randomColor();
			vec2i p0 = ml::math::round(ml::vec2f(kp.m_kp0.m_pixelPos.x, kp.m_kp0.m_pixelPos.y));
			vec2i p1 = ml::math::round(ml::vec2f(kp.m_kp1.m_pixelPos.x, kp.m_kp1.m_pixelPos.y));
			p1.x += img0.getWidth();

			ImageHelper::drawCircle(matchImage, p0, ml::math::round(kp.m_kp0.m_size), c.getVec3());
			ImageHelper::drawCircle(matchImage, p1, ml::math::round(kp.m_kp1.m_size), c.getVec3());
			ImageHelper::drawLine(matchImage, p0, p1, c.getVec3());
		}
		return matchImage;
	}

};
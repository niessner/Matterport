
#pragma once

#include "mLibInclude.h"

#include "keyPoint.h"
#include "imageHelper.h"

class MatchVisualization {
public:
	ColorImageR8G8B8 match(ColorImageR8G8B8& img0, ColorImageR8G8B8& img1, const std::vector<KeyPointMatch>& matches) {

		ColorImageR8G8B8 matchImage(img0.getWidth() * 2, img0.getHeight());
		matchImage.copyIntoImage(img0, 0, 0);
		matchImage.copyIntoImage(img1, img1.getWidth(), 0);

		RGBColor lowColor = ml::RGBColor::Blue;
		RGBColor highColor = ml::RGBColor::Red;
		for (size_t i = 0; i < matches.size(); i++) {
			const KeyPointMatch& kp = matches[i];

			RGBColor c = RGBColor::interpolate(lowColor, highColor, 0.5f);
			vec3f color(c.r / 255.0f, c.g / 255.0f, c.b / 255.0f);
			vec2i p0 = ml::math::round(ml::vec2f(kp.m_kp0.m_pixelPos.x, kp.m_kp0.m_pixelPos.y));
			vec2i p1 = ml::math::round(ml::vec2f(kp.m_kp1.m_pixelPos.x, kp.m_kp1.m_pixelPos.y));
			ImageHelper::drawCircle(matchImage, p0, ml::math::round(kp.m_kp0.m_size), color);
			ImageHelper::drawCircle(matchImage, p1, ml::math::round(kp.m_kp1.m_size), color);
			ImageHelper::drawLine(matchImage, p0, p1, color);
		}
		//std::cout << "(" << imageIndices << "): max match distance = " << maxMatchDistance << std::endl;
		//FreeImageWrapper::saveImage(filename, matchImage);

		return matchImage;
	}
private:
};
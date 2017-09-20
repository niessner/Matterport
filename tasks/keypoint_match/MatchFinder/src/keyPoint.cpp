
#include "stdafx.h"

#include "keyPoint.h"


// TODO: reference additional headers your program requires here
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/imgproc/imgproc.hpp"



template<>
struct std::hash<vec2ui> : public std::unary_function < vec2ui, size_t > {
	size_t operator()(const vec2ui& v) const {
		//TODO larger prime number (64 bit) to match size_t
		const size_t p0 = 73856093;
		const size_t p1 = 19349669;
		//const size_t p2 = 83492791;
		const size_t res = ((size_t)v.x * p0) ^ ((size_t)v.y * p1);// ^ ((size_t)v.z * p2);
		return res;
	}
};

template<>
struct std::hash<vec3ui> : public std::unary_function < vec3ui, size_t > {
	size_t operator()(const vec3ui& v) const {
		//TODO larger prime number (64 bit) to match size_t
		const size_t p0 = 73856093;
		const size_t p1 = 19349669;
		const size_t p2 = 83492791;
		const size_t res = ((size_t)v.x * p0) ^ ((size_t)v.y * p1) ^ ((size_t)v.z * p2);
		return res;
	}
};


float computeLuminance(const vec3uc& c) {
	return (0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z) / 255.f;
}


void toItensityImageAsMat(const ColorImageR8G8B8& image, cv::Mat& im)
{
	im.create(image.getHeight(), image.getWidth(), CV_8UC1);

	//unsigned char* data = new unsigned char[image.size()];
	for (unsigned int i = 0; i < image.size(); i++) {
		float lum = computeLuminance(image.getData()[i]);
		//data[i] = (unsigned char)(255.f * lum);
		im.data[i] = math::round((255.f * lum));
	}
	//im = cv::Mat(image.getHeight(), image.getWidth(), CV_8UC1, data);
	//MLIB_ASSERT(im.data != nullptr);

	//SAFE_DELETE_ARRAY(data);
}

std::vector<KeyPoint> KeyPointFinder::findKeyPoints(const vec2ui& sensImIdxs, const ColorImageR8G8B8& image, unsigned int maxNumKeyPoints, float minResponse) {
	std::vector<KeyPoint> res;

	cv::Mat im;
	toItensityImageAsMat(image, im);

	//cv::imwrite("Gray_Image.jpg", im);
	//cv::namedWindow("Gray image", CV_WINDOW_AUTOSIZE);
	//cv::imshow("Gray image", im);

	// sift detector
	cv::SiftFeatureDetector detector;
	detector.set("nFeatures", (int)maxNumKeyPoints); //TODO set first octave to 1 here?
	std::vector<cv::KeyPoint> keypoints;
	
	detector.detect(im, keypoints);
	res.reserve(keypoints.size());
	
	// remove keypoints with no corresponding depth, same keypoints
	std::unordered_set<vec2ui, std::hash<vec2ui>> keyPointLocations;
	unsigned int sameCount = 0;
	unsigned int noDepthCount = 0;
	for (unsigned int i = 0; i < keypoints.size(); i++) {
		const cv::KeyPoint& k = keypoints[i];

		//std::cout << "response: " << k.response << std::endl;
		if (k.response < minResponse) continue;

		vec2ui loc((unsigned int)std::round(k.pt.x), (unsigned int)std::round(k.pt.y));
		if (keyPointLocations.find(loc) == keyPointLocations.end()) {
			keyPointLocations.insert(loc);
			//need to unpack octave from keypoint
			int octave = k.octave & 255;
			int layer = (k.octave >> 8) & 255;
			octave = octave < 128 ? octave : (-128 | octave);
			float scale = octave >= 0 ? 1.f / (1 << octave) : (float)(1 << -octave);
			res.push_back(KeyPoint(sensImIdxs.x, sensImIdxs.y, vec2f(k.pt.x, k.pt.y), k.size, k.angle, octave, scale, k.response, k.octave)); 
		}
	}

	return res;
}

#include "stdafx.h"

#include "keyPoint.h"
#include "matchVisualization.h"
#include "omp.h"

// TODO: reference additional headers your program requires here
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/imgproc/imgproc.hpp"


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


void KeyPointMatcher::matchKeyPoints(const std::vector<std::vector<ColorImageR8G8B8>>& images, const std::vector<KeyPoint>& keyPoints,
	const std::vector<KeyPointMatch>& keysToMatch, std::vector<float>& matchDists, const std::string& featureType)
{
	cv::initModule_nonfree(); //otherwise errors happen

	//compute descriptors for each image's keypoints
	std::unordered_map<vec2ui, std::vector<unsigned int>> keyIndicesPerSensImage;
	for (unsigned int i = 0; i < keyPoints.size(); i++) {
		const auto& k = keyPoints[i];
		const vec2ui sensImId(k.m_sensorIdx, k.m_imageIdx);
		auto it = keyIndicesPerSensImage.find(sensImId);
		if (it == keyIndicesPerSensImage.end()) {
			keyIndicesPerSensImage[sensImId] = std::vector<unsigned int>(1, i);
		}
		else {
			it->second.push_back(i);
		}
	}

	cv::Ptr<cv::DescriptorExtractor> extractor = cv::DescriptorExtractor::create(featureType);
	MLIB_ASSERT(extractor != NULL);
	//if (featureType == "ORB") { extractor->set("edgeThreshold", 10); }
	//(int nfeatures = 500, float scaleFactor = 1.2f, int nlevels = 8, int edgeThreshold = 31, int firstLevel = 0, int WTA_K = 2, int scoreType = ORB::HARRIS_SCORE, int patchSize = 31 )
	//( int nfeatures=0, int nOctaveLayers=3,double contrastThreshold = 0.04, double edgeThreshold = 10, double sigma = 1.6);

	//cv::SiftDescriptorExtractor extractor;
	std::vector<cv::Mat> descriptors(keyPoints.size());
	////debugging (only with 1 sensor)
	//std::vector<cv::Mat> descPerImage(50);
	//std::vector<std::vector<unsigned int>> keyIndicesPerImage(descPerImage.size());
	//for (const auto& im : keyIndicesPerSensImage) {
	//	descPerImage[im.first.y].create(im.second.size(), 128, CV_32F);
	//	keyIndicesPerImage[im.first.y] = im.second;
	//}
	////debugging
	for (const auto& im : keyIndicesPerSensImage) {
		std::vector<cv::KeyPoint> imKeypoints(im.second.size());
		for (unsigned int i = 0; i < im.second.size(); i++) {
			const unsigned int idx = im.second[i];
			//int octave;
			if (featureType == "SIFT")
				imKeypoints[i] = cv::KeyPoint(keyPoints[idx].m_pixelPos.x, keyPoints[idx].m_pixelPos.y, keyPoints[idx].m_size,
				keyPoints[idx].m_angle, 1.0f, keyPoints[idx].m_opencvPackOctave);
			else if (featureType == "ORB")
				imKeypoints[i] = cv::KeyPoint(keyPoints[idx].m_pixelPos.x, keyPoints[idx].m_pixelPos.y, keyPoints[idx].m_size,
				keyPoints[idx].m_angle);
			else if (featureType == "SURF")
				imKeypoints[i] = cv::KeyPoint(keyPoints[idx].m_pixelPos.x, keyPoints[idx].m_pixelPos.y, keyPoints[idx].m_size,
				keyPoints[idx].m_angle, 1.0f, keyPoints[idx].m_octave);
			////else throw MLIB_EXCEPTION("error: invalid feature type for octave");
			//imKeypoints[i] = cv::KeyPoint(keyPoints[idx].m_pixelPos.x, keyPoints[idx].m_pixelPos.y, keyPoints[idx].m_size,
			//	keyPoints[idx].m_angle);//, 1.0f, octave);
		}
		cv::Mat mat;
		toItensityImageAsMat(images[im.first.x][im.first.y], mat);

		const auto tmp = imKeypoints;

		//cv::imwrite("test.jpg", mat);
		cv::Mat imDescriptors; //rows->#keys, cols->desc size (128)
		extractor->compute(mat, imKeypoints, imDescriptors);
		if (imKeypoints.size() != im.second.size()) {
			std::cout << "warning: only computed " << imKeypoints.size() << " descriptors for " << im.second.size() << " keypoints of image " << im.first.y << " (sens idx " << im.first.x << ")" << std::endl;
			if (imKeypoints.empty())
				continue;
		}
		unsigned int ii = 0;
		for (unsigned int i = 0; i < imKeypoints.size(); i++) {
			unsigned int idx = im.second[ii];
			while (!(imKeypoints[i].pt.x == keyPoints[idx].m_pixelPos.x && imKeypoints[i].pt.y == keyPoints[idx].m_pixelPos.y)) {
				ii++;
				if (ii == im.second.size()) {
					throw MLIB_EXCEPTION("ERROR SFDLKJ");
					break;
				}
				idx = im.second[ii];
			}
			if (ii < im.second.size()) {
				imDescriptors.row(i).copyTo(descriptors[idx]); //assign desc rows to global descriptors array
				ii++;
			}

			//imDescriptors.row(i).copyTo(descPerImage[im.first.y].row(i)); //debugging
		}
	}

	////debug match 
	//{
	//	cv::FlannBasedMatcher matcher;
	//	std::vector<cv::DMatch> matches;
	//	matcher.match(descPerImage[7], descPerImage[8], matches);
	//	std::vector<KeyPointMatch> visMatches(std::min((int)matches.size(), 5));
	//	for (unsigned int i = 0; i < visMatches.size(); i++) {
	//		visMatches[i].m_kp0 = keyIndicesPerImage[7][matches[i].queryIdx];
	//		visMatches[i].m_kp1 = keyIndicesPerImage[8][matches[i].trainIdx];
	//	}
	//	MatchVisualization::visulizeMatches(images, visMatches, keyPoints, 10);
	//	int a = 5;
	//}
	////debug match 

	matchDists.resize(keysToMatch.size());
	//std::vector<cv::Ptr<cv::DescriptorMatcher>> matchers(omp_get_max_threads());
	//for (auto& m : matchers) {
	//	if (extractor->descriptorType() == CV_8U) m = cv::DescriptorMatcher::create("BruteForce-Hamming");
	//	else if (extractor->descriptorType() == CV_32F) m = cv::DescriptorMatcher::create("FlannBased");
	//	else throw MLIB_EXCEPTION("invalid descriptor type for descriptormatcher");
	//}
	//#pragma omp parallel for
	cv::Ptr<cv::DescriptorMatcher> matcher;
	if (extractor->descriptorType() == CV_8U) matcher = cv::DescriptorMatcher::create("BruteForce-Hamming");
	else if (extractor->descriptorType() == CV_32F) matcher = cv::DescriptorMatcher::create("FlannBased");
	else throw MLIB_EXCEPTION("invalid descriptor type for descriptormatcher");
	for (int i = 0; i < (int)keysToMatch.size(); i++) {
		//cv::FlannBasedMatcher matcher;
		//int thread = omp_get_thread_num();
		std::vector<cv::DMatch> matches;
		if (descriptors[keysToMatch[i].m_kp0].rows > 0 && descriptors[keysToMatch[i].m_kp1].rows > 0) {
			matcher->match(descriptors[keysToMatch[i].m_kp0], descriptors[keysToMatch[i].m_kp1], matches);
			matchDists[i] = matches.front().distance;

			//MatchVisualization::visulizeMatches(images, std::vector<KeyPointMatch>(1, keysToMatch[i]), keyPoints, 10);
			//int a = 5;
		}
	}
}

void KeyPointMatcher::debug(const std::vector<std::vector<ColorImageR8G8B8>>& images)
{
	//just uses the first sens file
	std::vector<cv::Mat> cvImages(images.front().size());
	for (unsigned int i = 0; i < images.front().size(); i++) {
		toItensityImageAsMat(images[0][i], cvImages[i]);
	}

	//detect and extract
	cv::SiftFeatureDetector detector;
	std::vector<std::vector<cv::KeyPoint>> keypoints(cvImages.size());
	cv::SiftDescriptorExtractor extractor;
	std::vector<cv::Mat> descriptors(cvImages.size());
	for (unsigned int i = 0; i < cvImages.size(); i++) {
		detector.detect(cvImages[i], keypoints[i]);
		extractor.compute(cvImages[i], keypoints[i], descriptors[i]);
	}

	//match
	cv::FlannBasedMatcher matcher;
	std::vector<std::pair<vec2ui, std::vector<cv::DMatch>>> matchesPerImagePair;
	std::vector<float> allMatchDists;
	for (unsigned int i = 0; i < cvImages.size(); i++) {
		for (unsigned int j = i + 1; j < cvImages.size(); j++) {
			if (descriptors[i].rows > 0 && descriptors[j].rows > 0) {
				std::vector<cv::DMatch> matches;
				matcher.match(descriptors[i], descriptors[j], matches);
				matchesPerImagePair.push_back(std::make_pair(vec2ui(i, j), matches));
				for (const auto& m : matches) allMatchDists.push_back(m.distance);

				//std::vector<cv::DMatch> goodMatches;
				//for (const auto& m : matches) {
				//	if (m.distance < 200)
				//		goodMatches.push_back(m);
				//}
				//if (!goodMatches.empty()) {
				//	std::vector<KeyPoint> keys(keypoints[i].size() + keypoints[j].size());
				//	for (unsigned int k = 0; k < keypoints[i].size(); k++) {
				//		keys[k].m_sensorIdx = 0;
				//		keys[k].m_imageIdx = i;
				//		keys[k].m_pixelPos = vec2f(keypoints[i][k].pt.x, keypoints[i][k].pt.y);
				//		keys[k].m_size = keypoints[i][k].size;
				//	}
				//	for (unsigned int k = 0; k < keypoints[j].size(); k++) {
				//		keys[keypoints[i].size() + k].m_sensorIdx = 0;
				//		keys[keypoints[i].size() + k].m_imageIdx = j;
				//		keys[keypoints[i].size() + k].m_pixelPos = vec2f(keypoints[j][k].pt.x, keypoints[j][k].pt.y);
				//		keys[keypoints[i].size() + k].m_size = keypoints[j][k].size;
				//	}
				//	std::vector<KeyPointMatch> visMatches(goodMatches.size());
				//	for (unsigned int k = 0; k < visMatches.size(); k++) {
				//		const auto& m = goodMatches[k];
				//		visMatches[k].m_kp0 = m.queryIdx;
				//		visMatches[k].m_kp1 = (unsigned int)keypoints[i].size() + m.trainIdx;
				//	}
				//	MatchVisualization::visulizeMatches(images, visMatches, keys, 10);
				//	int a = 5;
				//}
			}
		}
	}

	std::sort(allMatchDists.begin(), allMatchDists.end());
	std::cout << "match dist range: " << allMatchDists.front() << " : " << allMatchDists.back() << std::endl;
	std::cout << "avg match dist: " << std::accumulate(allMatchDists.begin(), allMatchDists.end(), 0.0f) / (float)allMatchDists.size() << std::endl;
	size_t numMatchesBelow150 = 0;
	for (const auto& m : allMatchDists) {
		if (m < 150) numMatchesBelow150++;
		else break;
	}
	std::cout << "#match dists below 150 = " << numMatchesBelow150 << " (" << allMatchDists.size() << " total)" << std::endl;

	int a = 5;
	std::cout << "waiting..." << std::endl; getchar();
}

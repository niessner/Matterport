
#pragma once

#include "mLibInclude.h"

#ifndef VAR_NAME
#define VAR_NAME(x) #x
#endif

class KeyPoint {
public:
	unsigned int m_sensorIdx;
	unsigned int m_imageIdx;
	vec2f m_pixelPos;
	float m_depth;
	vec3f m_worldPos;
	//vec3f m_worldNormal;
	float m_size;		//from the sift keypoint extractor
	float m_angle;		//from the sift keypoint extractor
	int m_octave;		//from the sift keypoint extractor
	float m_scale;		//from the sift keypoint extractor
	float m_response;	//from the sift keypoint extractor
	int m_opencvPackOctave; //from the sift keypoint extractor (makes it easier to re-run opencv)

	bool isSameImage(const KeyPoint& other) const {
		if (m_sensorIdx == other.m_sensorIdx &&
			m_imageIdx == other.m_imageIdx) return true;
		else return false;
	}

	inline bool operator==(const KeyPoint& other) const {
		if (m_sensorIdx == other.m_sensorIdx && m_imageIdx == other.m_imageIdx && 
			m_pixelPos == other.m_pixelPos && m_depth == other.m_depth &&
			m_size == other.m_size && m_angle == other.m_angle)
			return true;

		return false;
	}
};

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

//warning: doesn't take into account sensor idx
template<>
struct std::hash<KeyPoint> : public std::unary_function < KeyPoint, size_t > {
	size_t operator()(const KeyPoint& v) const {
		//TODO larger prime number (64 bit) to match size_t
		const size_t p0 = 73856093;
		const size_t p1 = 19349669;
		const size_t p2 = 83492791;
		const size_t res = ((size_t)v.m_imageIdx * p0) ^ ((size_t)math::round(v.m_pixelPos.x) * p1) ^ ((size_t)math::round(v.m_pixelPos.y) * p2);
		return res;
	}
};


inline std::ostream& operator<<(std::ostream& os, const KeyPoint& kp) {
	os 
		<< "\t" << VAR_NAME(kp.m_sensorIdx) << " = " << kp.m_sensorIdx << "\n"
		<< "\t" << VAR_NAME(kp.m_imageIdx) << " = " << kp.m_imageIdx << "\n"
		<< "\t" << VAR_NAME(kp.m_pixelPos) << " = " << kp.m_pixelPos << "\n"
		<< "\t" << VAR_NAME(kp.m_depth) << " = " << kp.m_depth << "\n"
		<< "\t" << VAR_NAME(kp.m_worldPos) << " = " << kp.m_worldPos << std::endl;
	return os;
}


class KeyPointMatch {
public:
	size_t m_kp0; //indexes into the keypoints array
	size_t m_kp1; //indexes into the keypoints array
	vec2f m_offset;	//re-projection offset when m_kp0 is projected into m_kp1;

	KeyPointMatch() {
		m_kp0 = (size_t)-1;
		m_kp1 = (size_t)-1;
	}

	bool isSameImagePair(const KeyPointMatch& other, const std::vector<KeyPoint>& keys) const {
		if (keys[m_kp0].m_sensorIdx == keys[other.m_kp0].m_sensorIdx &&
			keys[m_kp0].m_imageIdx == keys[other.m_kp0].m_imageIdx &&
			keys[m_kp1].m_sensorIdx == keys[other.m_kp1].m_sensorIdx &&
			keys[m_kp1].m_imageIdx == keys[other.m_kp1].m_imageIdx)
			return true;
		else return false;
	}
};

class KeyPointMatcher {
public:
	static void matchKeyPoints(const std::vector<std::vector<ColorImageR8G8B8>>& images, const std::vector<KeyPoint>& keyPoints,
		const std::vector<KeyPointMatch>& keysToMatch, std::vector<float>& matchDists, const std::string& featureType);

	//run detect-extract-match on the images
	static void debug(const std::vector<std::vector<ColorImageR8G8B8>>& images);
private:
};

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
	vec3f m_worldNormal;
	float m_size;		//from the sift keypoint extractor
	float m_angle;		//from the sift keypoint extractor
	int m_octave;		//from the sift keypoint extractor
	float m_scale;		//from the sift keypoint extractor
	float m_response;	//from the sift keypoint extractor
	int m_opencvPackOctave; //from the sift keypoint extractor (makes it easier to re-run opencv)

	KeyPoint()
	{
		m_sensorIdx = (unsigned int)-1;
		m_imageIdx = (unsigned int)-1;
		m_pixelPos = vec2f(-std::numeric_limits<float>::infinity());
		m_depth = -std::numeric_limits<float>::infinity();
		m_worldPos = vec3f(-std::numeric_limits<float>::infinity());
		m_worldNormal = vec3f(-std::numeric_limits<float>::infinity());
		m_size = -std::numeric_limits<float>::infinity();
		m_angle = -std::numeric_limits<float>::infinity();
		m_octave = -1;
		m_opencvPackOctave = -1;
		m_scale = -std::numeric_limits<float>::infinity();
		m_response = -std::numeric_limits<float>::infinity();
	}

	KeyPoint(unsigned int sensorIdx, unsigned int imageIdx, const vec2f& pixelPos,
		float size, float angle, int octave, float scale, float response, int packOctave)
		: m_sensorIdx(sensorIdx), m_imageIdx(imageIdx), m_pixelPos(pixelPos), m_size(size),
		m_angle(angle), m_scale(scale), m_octave(octave), m_response(response), m_opencvPackOctave(packOctave)
	{
		m_depth = -std::numeric_limits<float>::infinity();
		m_worldPos = vec3f(-std::numeric_limits<float>::infinity());
		m_worldNormal = vec3f(-std::numeric_limits<float>::infinity());
	}

	//KeyPoint(unsigned int sensorIdx, unsigned int imageIdx, const vec2f& pixelPos,
	//	float depth, const vec3f& worldPos, const vec3f& worldNormal,
	//	float size, float angle, int octave, float scale, float response, int packOctave)
	//	: m_sensorIdx(sensorIdx), m_imageIdx(imageIdx), m_pixelPos(pixelPos),
	//	m_depth(depth), m_worldPos(worldPos), m_worldNormal(worldNormal),
	//	m_size(size), m_angle(angle), m_scale(scale), m_octave(octave), m_response(response), m_opencvPackOctave(packOctave)
	//{}

	bool isSameImage(const KeyPoint& other) const {
		if (m_sensorIdx == other.m_sensorIdx &&
			m_imageIdx == other.m_imageIdx) return true;
		else return false;
	}
};

class KeyPointFinder {
public:
	//! returns the key points (x, y, size, response)
	static std::vector<KeyPoint> findKeyPoints(const vec2ui& sensImIdxs, const ColorImageR8G8B8& image, unsigned int maxNumKeyPoints, float minResponse);
private:
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
	KeyPoint m_kp0;
	KeyPoint m_kp1;
	vec2f m_offset;	//re-projection offset when m_kp0 is projected into m_kp1;

	bool isSameImagePair(const KeyPointMatch& other) const {
		if (m_kp0.m_sensorIdx == other.m_kp0.m_sensorIdx &&
			m_kp0.m_imageIdx == other.m_kp0.m_imageIdx &&
			m_kp1.m_sensorIdx == other.m_kp1.m_sensorIdx &&
			m_kp1.m_imageIdx == other.m_kp1.m_imageIdx)
			return true;
		else return false;
	}
};

inline std::ostream& operator<<(std::ostream& os, const KeyPointMatch& kpm) {
	os << VAR_NAME(kpm.m_kp0) << "\n" << kpm.m_kp0;
	os << VAR_NAME(kpm.m_kp1) << "\n" << kpm.m_kp1;
	os << VAR_NAME(kpm.m_offset) << " = " << kpm.m_offset << std::endl;
	return os;
}


#pragma once

class ImageHelper
{
public:
	ImageHelper() {}
	~ImageHelper() {}


	static inline float gaussR(float sigma, float dist)
	{
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}
	static inline float gaussR(float sigma, const vec3f& d)
	{
		float dist = d.length();
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}
	static inline float gaussR(float sigma, const vec3uc& d)
	{
		vec3f _d(d);	//_d /= 255.0f;
		float dist = _d.length();
		return exp(-(dist*dist) / (2.0f*sigma*sigma));
	}


	static inline float linearR(float sigma, float dist)
	{
		return std::max(1.0f, std::min(0.0f, 1.0f - (dist*dist) / (2.0f*sigma*sigma)));
	}

	static inline float gaussD(float sigma, int x, int y)
	{
		return exp(-((x*x + y*y) / (2.0f*sigma*sigma)));
	}

	static inline float gaussD(float sigma, int x)
	{
		return exp(-((x*x) / (2.0f*sigma*sigma)));
	}

	static void bilateralFilter(BaseImage<vec3uc>& img, float sigmaD, float sigmaR) {

		BaseImage<vec3uc> res(img.getDimensions());
		res.setInvalidValue(img.getInvalidValue());

		const int kernelRadius = (int)ceil(2.0*sigmaD);
		for (unsigned int y = 0; y < img.getHeight(); y++) {
			for (unsigned int x = 0; x < img.getWidth(); x++) {

				res.setInvalid(x, y);

				vec3f sum = vec3f(0.0f);
				float sumWeight = 0.0f;

				if (img.isValid(x, y)) {
					const vec3uc& center = img(x, y);

					for (int m = x - kernelRadius; m <= (int)x + kernelRadius; m++) {
						for (int n = y - kernelRadius; n <= (int)y + kernelRadius; n++) {
							if (m >= 0 && n >= 0 && m < (int)img.getWidth() && n < (int)img.getHeight()) {
								if (img.isValid(m, n)) {
									const vec3uc& current = img(m, n);
									const float weight = gaussD(sigmaD, m - x, n - y)*gaussR(sigmaR, current - center);
									sumWeight += weight;
									sum += weight*vec3f(current);
								}
							}
						}
					}
					if (sumWeight > 0.0f) res(x, y) = math::round(sum / sumWeight);
				}
			}
		}
		img = res;
	}
	static void bilateralFilter(DepthImage32& d, float sigmaD, float sigmaR) {

		DepthImage32 res(d.getWidth(), d.getHeight());
		res.setInvalidValue(d.getInvalidValue());

		const int kernelRadius = (int)ceil(2.0*sigmaD);
		for (unsigned int y = 0; y < d.getHeight(); y++) {
			for (unsigned int x = 0; x < d.getWidth(); x++) {
				res.setInvalid(x, y);

				float sum = 0.0f;
				float sumWeight = 0.0f;

				if (d.isValid(x, y)) {
					const float center = d(x, y);

					for (int m = x - kernelRadius; m <= (int)x + kernelRadius; m++) {
						for (int n = y - kernelRadius; n <= (int)y + kernelRadius; n++) {
							if (m >= 0 && n >= 0 && m < (int)d.getWidth() && n < (int)d.getHeight()) {
								if (d.isValid(m, n)) {
									const float current = d(m, n);
									const float weight = gaussD(sigmaD, m - x, n - y)*gaussR(sigmaR, current - center);
									sumWeight += weight;
									sum += weight*current;
								}
							}
						}
					}
					if (sumWeight > 0.0f) res(x, y) = sum / sumWeight;
				}
			}
		}
		d = res;
	}

	static void erode(DepthImage32& depth, unsigned int numIter = 2) {
		numIter = 2 * ((numIter + 1) / 2);
		DepthImage32 tmp; tmp.setInvalidValue(depth.getInvalidValue()); 
		for (unsigned int i = 0; i < numIter; i++) {
			if (i % 2 == 0) {
				erode(tmp, depth, 3, 0.05f, 0.3f);
			}
			else {
				erode(depth, tmp, 3, 0.05f, 0.3f);
			}
		}
	}

	// draw functions
	template<typename T>
	static void drawCircle(ml::BaseImage<T>& image, const ml::vec2f& center, float radius, const T& color) {
		drawCircle(image, ml::math::round(center), ml::math::round(radius), color);
	}
	template<typename T>
	static void drawCircle(ml::BaseImage<T>& image, const ml::vec2i& center, int radius, const T& color) {
		int x = radius;
		int y = 0;
		int radiusError = 1 - x;

		while (x >= y) {
			if (image.isValidCoordinate(center.x + x, center.y + y)) image(center.x + x, center.y + y) = color;
			if (image.isValidCoordinate(center.x + y, center.y + x)) image(center.x + y, center.y + x) = color;
			if (image.isValidCoordinate(center.x - x, center.y + y)) image(center.x - x, center.y + y) = color;
			if (image.isValidCoordinate(center.x - y, center.y + x)) image(center.x - y, center.y + x) = color;
			if (image.isValidCoordinate(center.x - x, center.y - y)) image(center.x - x, center.y - y) = color;
			if (image.isValidCoordinate(center.x - y, center.y - x)) image(center.x - y, center.y - x) = color;
			if (image.isValidCoordinate(center.x + x, center.y - y)) image(center.x + x, center.y - y) = color;
			if (image.isValidCoordinate(center.x + y, center.y - x)) image(center.x + y, center.y - x) = color;
			y++;
			if (radiusError < 0) {
				radiusError += 2 * y + 1;
			}
			else {
				x--;
				radiusError += 2 * (y - x) + 1;
			}
		}
	}
	template<typename T>
	static void drawLine(ml::BaseImage<T>& image, const ml::vec2i& start, const ml::vec2i& end, const T& color) {
		ml::vec2i s, e;
		s.x = ml::math::clamp(start.x, 0, (int)image.getWidth() - 1);
		s.y = ml::math::clamp(start.y, 0, (int)image.getHeight() - 1);
		e.x = ml::math::clamp(end.x, 0, (int)image.getWidth() - 1);
		e.y = ml::math::clamp(end.y, 0, (int)image.getHeight() - 1);
		MLIB_ASSERT(image.isValidCoordinate(s.x, s.y) && image.isValidCoordinate(e.x, e.y));

		if (end.x >= start.x)
			drawLineInternal(image, s, e, color);
		else
			drawLineInternal(image, e, s, color);
	}
	 
private:
	//! assumes end.x >= start.x
	template<typename T>
	static void drawLineInternal(ml::BaseImage<T>& image, const ml::vec2i& start, const ml::vec2i& end, const T& color) {
		//MLIB_ASSERT(image.isValidCoordinate(start.x, start.y) && image.isValidCoordinate(end.x, end.y));

		float dx = (float)(end.x - start.x);
		if (dx == 0.0f) { // vertical line
			for (int y = start.y; y < end.y; y++) image(start.x, y) = color;
			return;
		}
		float dy = (float)(end.y - start.y);
		float error = 0.0f;
		float dErr = ml::math::abs(dy / dx);
		int dir = ml::math::sign(end.y - start.y);
		int y = start.y;
		for (int x = start.x; x < end.x; x++) {
			image(x, y) = color;
			error += dErr;
			while (error >= 0.5f) {
				image(x, y) = color;
				y += dir;
				error--;
			}
		}
	}

	static void erode(DepthImage32& output, const DepthImage32& input, int structureSize, float dThresh, float fracReq)
	{
		output.allocate(input.getWidth(), input.getHeight());

		for (unsigned int y = 0; y < input.getHeight(); y++) {
			for (unsigned int x = 0; x < input.getWidth(); x++) {
				unsigned int count = 0;
				float oldDepth = input(x, y);
				for (int i = -structureSize; i <= structureSize; i++) {
					for (int j = -structureSize; j <= structureSize; j++) {
						if (x + j >= 0 && x + j < input.getWidth() && y + i >= 0 && y + i < input.getHeight()) {
							float depth = input(x + j, y + i);
							if (depth == input.getInvalidValue() || depth == 0.0f || fabs(depth - oldDepth) > dThresh) {
								count++;
							}
						}
					}
				}

				unsigned int sum = (2 * structureSize + 1)*(2 * structureSize + 1);
				if ((float)count / (float)sum >= fracReq) {
					output(x, y) = input.getInvalidValue();
				}
				else {
					output(x, y) = input(x, y);
				}
			}
		}
	}
};


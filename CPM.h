/*

Code of the Coarse-to-Fine PatchMatch, published at CVPR 2016 in
"Efficient Coarse-to-Fine PatchMatch for Large Displacement Optical Flow"
by Yinlin.Hu, Rui Song and Yunsong Li.

Email: huyinlin@gmail.com

Version 1.2

Copyright (C) 2016 Yinlin.Hu

Usages:

The program "cpm.exe" has been built and tested on Windows 7.

USAGE: cpm.exe img1Name img2Name outMatchName

Explanations:

The output of the program is a text file, which is in the format of "x1,y1,x2,y2"
corresponding to one match per line.

*/

#ifndef _CPM_H_
#define _CPM_H_

#include "ImagePyramid.h"

class CPM
{
public:
	CPM();
	~CPM();

	int Matching(FImage& img1, FImage& img2, FImage& outMatches);
	void SetStereoFlag(int needStereo);
	void SetStep(int step);

private:
	void imDaisy(FImage& img, UCImage& outFtImg);
	void CrossCheck(IntImage& seeds, FImage& seedsFlow, FImage& seedsFlow2, IntImage& kLabel2, int* valid, float th);
	float MatchCost(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, int x1, int y1, int x2, int y2);

	// a good initialization is already stored in bestU & bestV
	int Propogate(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* pyd1f, UCImage* pyd2f, int level, float* radius, int iterCnt, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow, float* bestCosts);
	void PyramidRandomSearch(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow);
	void OnePass(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage& seeds, IntImage& neighbors, FImage* pydSeedsFlow);
	void UpdateSearchRadius(IntImage& neighbors, FImage* pydSeedsFlow, int level, float* outRadius);

	// minimum circle
	struct Point{
		double x, y;
	};
	double dist(Point a, Point b);
	Point intersection(Point u1, Point u2, Point v1, Point v2);
	Point circumcenter(Point a, Point b, Point c);
	// return the radius of the minimal circle
	float MinimalCircle(float* x, float*y, int n, float* centerX = NULL, float* centerY = NULL);

	//
	int _step;
	int _maxIters;
	float _stopIterRatio;
	float _pydRatio;

	int _isStereo;
	int _maxDisplacement;
	float _checkThreshold;
	int _borderWidth;

	IntImage _kLabels, _kLabels2;

	FImagePyramid _pyd1;
	FImagePyramid _pyd2;

	UCImage* _im1f;
	UCImage* _im2f;

	FImage* _pydSeedsFlow;
	FImage* _pydSeedsFlow2;

	IntImage _seeds;
	IntImage _seeds2;
	IntImage _neighbors;
	IntImage _neighbors2;

};

#endif // _CPM_H_
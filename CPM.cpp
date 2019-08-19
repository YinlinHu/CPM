#include "CPM.h"
#include "ImageFeature.h"

#ifdef USE_DAISY
#include "opencv2/xfeatures2d.hpp" // for "DAISY" descriptor
#endif

// [4/6/2017 Yinlin.Hu]

#define UNKNOWN_FLOW 1e10

CPM::CPM()
{
	// default parameters
	_step = 3;
	_isStereo = false;

	_maxIters = 8;
	_stopIterRatio = 0.05;
	_pydRatio = 0.5;

	_maxDisplacement = 400;
	_checkThreshold = 3;
	_borderWidth = 5;

	_im1f = NULL;
	_im2f = NULL;
	_pydSeedsFlow = NULL;
	_pydSeedsFlow2 = NULL;
}

CPM::~CPM()
{
	if (_im1f)
		delete[] _im1f;
	if (_im2f)
		delete[] _im2f;
	if (_pydSeedsFlow)
		delete[] _pydSeedsFlow;
	if (_pydSeedsFlow2)
		delete[] _pydSeedsFlow2;
}

void CPM::SetStereoFlag(int needStereo)
{
	_isStereo = needStereo;
}

void CPM::SetStep(int step)
{
	_step = step;
}

int CPM::Matching(FImage& img1, FImage& img2, FImage& outMatches)
{
	CTimer t;

	int w = img1.width();
	int h = img1.height();

	_pyd1.ConstructPyramid(img1, _pydRatio, 30);
	_pyd2.ConstructPyramid(img2, _pydRatio, 30);

	int nLevels = _pyd1.nlevels();

	if (_im1f)
		delete[] _im1f;
	if (_im2f)
		delete[] _im2f;

	_im1f = new UCImage[nLevels];
	_im2f = new UCImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		#ifdef USE_DAISY
			imDaisy(_pyd1[i], _im1f[i]);
			imDaisy(_pyd2[i], _im2f[i]);
		#else
			ImageFeature::imSIFT(_pyd1[i], _im1f[i], 2, 1, true, 8);
			ImageFeature::imSIFT(_pyd2[i], _im2f[i], 2, 1, true, 8);
		#endif
	}
	t.toc("get feature: ");

	int step = _step;
	int gridw = w / step;
	int gridh = h / step;
	int xoffset = (w - (gridw - 1)*step) / 2;
	int yoffset = (h - (gridh - 1)*step) / 2;
	int numV = gridw * gridh;
	int numV2 = numV;

	if (_pydSeedsFlow)
		delete[] _pydSeedsFlow;
	if (_pydSeedsFlow2)
		delete[] _pydSeedsFlow2;
	_pydSeedsFlow = new FImage[nLevels];
	_pydSeedsFlow2 = new FImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		_pydSeedsFlow[i].allocate(2, numV);
		_pydSeedsFlow2[i].allocate(2, numV2);
	}

	_seeds.allocate(2, numV);
	_neighbors.allocate(12, numV);
	_neighbors.setValue(-1);
	int nbOffset[8][2] = { { 0, -1 }, { 0, 1 }, { 1, 0 }, { -1, 0 }, { -1, -1 }, { -1, 1 }, { 1, -1 }, { 1, 1 } };
	for (int i = 0; i < numV; i++){
		int gridX = i % gridw;
		int gridY = i / gridw;
		_seeds[2 * i] = gridX * step + xoffset;
		_seeds[2 * i + 1] = gridY * step + yoffset;
		int nbIdx = 0;
		for (int j = 0; j < 8; j++){
			int nbGridX = gridX + nbOffset[j][0];
			int nbGridY = gridY + nbOffset[j][1];
			if (nbGridX < 0 || nbGridX >= gridw || nbGridY < 0 || nbGridY >= gridh)
				continue;
			_neighbors[i*_neighbors.width() + nbIdx] = nbGridY*gridw + nbGridX;
			nbIdx++;
		}
	}
	_seeds2.copy(_seeds);
	_neighbors2.copy(_neighbors);

	FImage seedsFlow(2, numV);

	_kLabels.allocate(w, h);
	for (int i = 0; i < numV; i++){
		int x = _seeds[2 * i];
		int y = _seeds[2 * i + 1];
		int r = step / 2;
		for (int ii = -r; ii <= r; ii++){
			for (int jj = -r; jj <= r; jj++){
				int xx = ImageProcessing::EnforceRange(x + ii, w);
				int yy = ImageProcessing::EnforceRange(y + jj, h);
				_kLabels[yy*w + xx] = i;
			}
		}
	}
	_kLabels2.copy(_kLabels);
	//kLabels.imshow("kLabels", 0);

	//t.toc("generate seeds: ");

	t.tic();
	OnePass(_pyd1, _pyd2, _im1f, _im2f, _seeds, _neighbors, _pydSeedsFlow);
	t.toc("forward matching: ");
	OnePass(_pyd2, _pyd1, _im2f, _im1f, _seeds2, _neighbors2, _pydSeedsFlow2);
	t.toc("backward matching: ");

	// cross check
	int* validFlag = new int[numV];
	CrossCheck(_seeds, _pydSeedsFlow[0], _pydSeedsFlow2[0], _kLabels2, validFlag, _checkThreshold);
	seedsFlow.copyData(_pydSeedsFlow[0]);
	for (int i = 0; i < numV; i++){
		if (!validFlag[i]){
			seedsFlow[2 * i] = UNKNOWN_FLOW;
			seedsFlow[2 * i + 1] = UNKNOWN_FLOW;
		}
	}
	delete[] validFlag;

	// flow 2 match
	FImage tmpMatch(4, numV);
	tmpMatch.setValue(-1);
	int validMatCnt = 0;
	for (int i = 0; i < numV; i++){
		int x = _seeds[2 * i];
		int y = _seeds[2 * i + 1];
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		float x2 = x + u;
		float y2 = y + v;
		if (abs(u) < UNKNOWN_FLOW && abs(v) < UNKNOWN_FLOW){
			tmpMatch[4 * i + 0] = x;
			tmpMatch[4 * i + 1] = y;
			tmpMatch[4 * i + 2] = x2;
			tmpMatch[4 * i + 3] = y2;
			validMatCnt++;
		}
	}
	if (!outMatches.matchDimension(4, validMatCnt, 1)){
		outMatches.allocate(4, validMatCnt, 1);
	}
	int tmpIdx = 0;
	for (int i = 0; i < numV; i++){
		if (tmpMatch[4 * i + 0] >= 0){
			memcpy(outMatches.rowPtr(tmpIdx), tmpMatch.rowPtr(i), sizeof(int) * 4);
			tmpIdx++;
		}
	}

	return validMatCnt;
}

#ifdef USE_DAISY
void CPM::imDaisy(FImage& img, UCImage& outFtImg)
{
	FImage imgray;
	img.desaturate(imgray);

	int w = imgray.width();
	int h = imgray.height();

	// use the version in OpenCV
	cv::Ptr<cv::xfeatures2d::DAISY> daisy =
		cv::xfeatures2d::DAISY::create(5, 3, 4, 8, 
		cv::xfeatures2d::DAISY::NRM_FULL, cv::noArray(), false, false);
	cv::Mat cvImg(h, w, CV_8UC1);
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			cvImg.at<unsigned char>(i, j) = imgray[i*w + j] * 255;
		}
	}
	cv::Mat outFeatures;
	daisy->compute(cvImg, outFeatures);

	int itSize = outFeatures.cols;
	outFtImg.allocate(w, h, itSize);
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;
			for (int k = 0; k < itSize; k++){
				outFtImg.pData[idx*itSize + k] = outFeatures.at<float>(idx, k) * 255;
			}
		}
	}
}
#endif

void CPM::CrossCheck(IntImage& seeds, FImage& seedsFlow, FImage& seedsFlow2, IntImage& kLabel2, int* valid, float th)
{
	int w = kLabel2.width();
	int h = kLabel2.height();
	int numV = seeds.height();
	for (int i = 0; i < numV; i++){
		valid[i] = 1;
	}

	// cross check (1st step)
	int b = _borderWidth;
	for (int i = 0; i < numV; i++){
		float u = seedsFlow[2 * i];
		float v = seedsFlow[2 * i + 1];
		int x = seeds[2 * i];
		int y = seeds[2 * i + 1];
		int x2 = x + u;
		int y2 = y + v;
		if (x < b || x >= w - b || y < b || y >= h - b
			|| x2 < b || x2 >= w - b || y2 < b || y2 >= h - b
			|| sqrt(u*u + v*v)>_maxDisplacement){
			valid[i] = 0;
			continue;
		}

		int idx2 = kLabel2[y2*w + x2];
		float u2 = seedsFlow2[2 * idx2];
		float v2 = seedsFlow2[2 * idx2 + 1];
		float diff = sqrt((u + u2)*(u + u2) + (v + v2)*(v + v2));
		if (diff > th){
			valid[i] = 0;
		}

	}
}

float CPM::MatchCost(FImage& img1, FImage& img2, UCImage* im1f, UCImage* im2f, int x1, int y1, int x2, int y2)
{
	int w = im1f->width();
	int h = im1f->height();
	int ch = im1f->nchannels();
	float totalDiff;

	// fast
	x1 = ImageProcessing::EnforceRange(x1, w);
	x2 = ImageProcessing::EnforceRange(x2, w);
	y1 = ImageProcessing::EnforceRange(y1, h);
	y2 = ImageProcessing::EnforceRange(y2, h);

	unsigned char* p1 = im1f->pixPtr(y1, x1);
	unsigned char* p2 = im2f->pixPtr(y2, x2);

	totalDiff = 0;

#ifdef WITH_SSE
	// SSE2
	unsigned char *_p1 = p1, *_p2 = p2;
	hu_m128 r1, r2, r3;
	int iterCnt = ch / 16;
	int idx = 0;
	int sum0 = 0;
	int sum1 = 0;
	for (idx = 0; idx < iterCnt; idx++){
		memcpy(&r1, _p1, sizeof(hu_m128));
		memcpy(&r2, _p2, sizeof(hu_m128));
		_p1 += sizeof(hu_m128);
		_p2 += sizeof(hu_m128);
		r3.mi = _mm_sad_epu8(r1.mi, r2.mi);
		sum0 += r3.m128i_u16[0];
		sum1 += r3.m128i_u16[4];
	}
	totalDiff += sum0;
	totalDiff += sum1;
	// add the left
	for (idx *= 16; idx < ch; idx++){
		totalDiff += abs(p1[idx] - p2[idx]);
	}
#else
	totalDiff = 0;
	for (int idx = 0; idx < ch; idx++){
		totalDiff += abs(p1[idx] - p2[idx]);
	}
#endif

	return totalDiff;
}

int CPM::Propogate(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* pyd1f, UCImage* pyd2f, int level, float* radius, int iterCnt, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow, float* bestCosts)
{
	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	FImage im1 = pyd1[level];
	FImage im2 = pyd2[level];
	UCImage* im1f = pyd1f + level;
	UCImage* im2f = pyd2f + level;
	IntImage* seeds = pydSeeds + level;
	FImage* seedsFlow = pydSeedsFlow + level;

	int w = im1.width();
	int h = im1.height();
	int ptNum = seeds->height();

	int maxNb = neighbors.width();
	int* vFlags = new int[ptNum];

	// init cost
	for (int i = 0; i < ptNum; i++){
		int x = seeds->pData[2 * i];
		int y = seeds->pData[2 * i + 1];
		float u = seedsFlow->pData[2 * i];
		float v = seedsFlow->pData[2 * i + 1];
		bestCosts[i] = MatchCost(im1, im2, im1f, im2f, x, y, x + u, y + v);
	}

	int iter = 0;
	float lastUpdateRatio = 2;
	for (iter = 0; iter < _maxIters; iter++)
	{
		int updateCount = 0;

		memset(vFlags, 0, sizeof(int)*ptNum);

		int startPos = 0, endPos = ptNum, step = 1;
		if (iter % 2 == 1){
			startPos = ptNum - 1; endPos = -1; step = -1;
		}
		for (int pos = startPos; pos != endPos; pos += step){
			bool updateFlag = false;

			int idx = pos;

			int x = seeds->pData[2 * idx];
			int y = seeds->pData[2 * idx + 1];

			int* nbIdx = neighbors.rowPtr(idx);
			// Propagation: Improve current guess by trying instead correspondences from neighbors
			for (int i = 0; i < maxNb; i++){
				if (nbIdx[i] < 0){
					break;
				}
				if (!vFlags[nbIdx[i]]){ // unvisited yet
					continue;
				}
				float tu = seedsFlow->pData[2 * nbIdx[i]];
				float tv = seedsFlow->pData[2 * nbIdx[i] + 1];
				float cu = seedsFlow->pData[2 * idx];
				float cv = seedsFlow->pData[2 * idx + 1];
				if (abs(tu - cu) < 1e-6 && abs(tv - cv) < 1e-6){
					continue;
				}
				float tc = MatchCost(im1, im2, im1f, im2f, x, y, x + tu, y + tv);
				if (tc < bestCosts[idx]){
					bestCosts[idx] = tc;
					seedsFlow->pData[2 * idx] = tu;
					seedsFlow->pData[2 * idx + 1] = tv;
					updateFlag = true;
				}
			}

			// Random search: Improve current guess by searching in boxes
			// of exponentially decreasing size around the current best guess.
			for (int mag = radius[idx] + 0.5; mag >= 1; mag /= 2) {
				/* Sampling window */
				float tu = seedsFlow->pData[2 * idx] + rand() % (2 * mag + 1) - mag;

				float tv = 0;
				if (!_isStereo){
					tv = seedsFlow->pData[2 * idx + 1] + rand() % (2 * mag + 1) - mag;
				}

				float cu = seedsFlow->pData[2 * idx];
				float cv = seedsFlow->pData[2 * idx + 1];
				if (abs(tu - cu) < 1e-6 && abs(tv - cv) < 1e-6){
					continue;
				}

				float tc = MatchCost(im1, im2, im1f, im2f, x, y, x + tu, y + tv);
				if (tc < bestCosts[idx]){
					bestCosts[idx] = tc;
					seedsFlow->pData[2 * idx] = tu;
					seedsFlow->pData[2 * idx + 1] = tv;
					updateFlag = true;
				}
			}
			vFlags[idx] = 1;
			//ShowSuperPixelFlow(spt, img1, bestU, bestV, ptNum);

			if (updateFlag){
				updateCount++;
			}
		}
		//printf("iter %d: %f [s]\n", iter, t.toc());

		float updateRatio = float(updateCount) / ptNum;
		//printf("Update ratio: %f\n", updateRatio);
		if (updateRatio < _stopIterRatio || lastUpdateRatio - updateRatio < 0.01){
			iter++;
			break;
		}
		lastUpdateRatio = updateRatio;
	}

	delete[] vFlags;

	return iter;
}

void CPM::PyramidRandomSearch(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage* pydSeeds, IntImage& neighbors, FImage* pydSeedsFlow)
{
	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	FImage rawImg1 = pyd1[0];
	FImage rawImg2 = pyd2[0];
	srand(0);

	int w = rawImg1.width();
	int h = rawImg1.height();
	int numV = pydSeeds[0].height();

	float* bestCosts = new float[numV];
	float* searchRadius = new float[numV];

	// random Initialization on coarsest level
	int initR = _maxDisplacement * pow(ratio, nLevels - 1) + 0.5;
	for (int i = 0; i < numV; i++){
		pydSeedsFlow[nLevels - 1][2 * i] = rand() % (2 * initR + 1) - initR;
		if (_isStereo){
			pydSeedsFlow[nLevels - 1][2 * i + 1] = 0;
		}else{
			pydSeedsFlow[nLevels - 1][2 * i + 1] = rand() % (2 * initR + 1) - initR;
		}
	}

	// set the radius of coarsest level
	for (int i = 0; i < numV; i++){
		searchRadius[i] = initR;
	}

	int* iterCnts = new int[nLevels];
	for (int i = 0; i < nLevels; i++){
		iterCnts[i] = _maxIters;
	}

	for (int l = nLevels - 1; l >= 0; l--){ // coarse-to-fine
		int iCnt = Propogate(pyd1, pyd2, im1f, im2f, l, searchRadius, iterCnts[l], pydSeeds, neighbors, pydSeedsFlow, bestCosts);
		if (l > 0){
			UpdateSearchRadius(neighbors, pydSeedsFlow, l, searchRadius);

			// scale the radius accordingly
			int maxR = __min(32, _maxDisplacement * pow(ratio, l) + 0.5);
			for (int i = 0; i < numV; i++){
				searchRadius[i] = __max(__min(searchRadius[i], maxR), 1);
				searchRadius[i] *= (1. / _pydRatio);
			}

			pydSeedsFlow[l - 1].copyData(pydSeedsFlow[l]);
			pydSeedsFlow[l - 1].Multiplywith(1. / ratio);
		}
	}

	delete[] searchRadius;
	delete[] bestCosts;
	delete[] iterCnts;
}

void CPM::OnePass(FImagePyramid& pyd1, FImagePyramid& pyd2, UCImage* im1f, UCImage* im2f, IntImage& seeds, IntImage& neighbors, FImage* pydSeedsFlow)
{
	FImage rawImg1 = pyd1[0];
	FImage rawImg2 = pyd2[0];

	int nLevels = pyd1.nlevels();
	float ratio = pyd1.ratio();

	int numV = seeds.height();

	IntImage* pydSeeds = new IntImage[nLevels];
	for (int i = 0; i < nLevels; i++){
		pydSeeds[i].allocate(2, numV);
		int sw = pyd1[i].width();
		int sh = pyd1[i].height();
		for (int n = 0; n < numV; n++){
			pydSeeds[i][2 * n] = ImageProcessing::EnforceRange(seeds[2 * n] * pow(ratio, i), sw);
			pydSeeds[i][2 * n + 1] = ImageProcessing::EnforceRange(seeds[2 * n + 1] * pow(ratio, i), sh);
		}
	}

	PyramidRandomSearch(pyd1, pyd2, im1f, im2f, pydSeeds, neighbors, pydSeedsFlow);

	// scale
	int b = _borderWidth;
	for (int i = 0; i < nLevels; i++){
		pydSeedsFlow[i].Multiplywith(pow(1. / ratio, i));
	}

	delete[] pydSeeds;
}

void CPM::UpdateSearchRadius(IntImage& neighbors, FImage* pydSeedsFlow, int level, float* outRadius)
{
	FImage* seedsFlow = pydSeedsFlow + level;
	int maxNb = neighbors.width();

	float x[32], y[32]; // for minimal circle
	assert(maxNb < 32);

	int sCnt = seedsFlow->height();
	for (int i = 0; i < sCnt; i++){
		// add itself
		x[0] = seedsFlow->pData[2 * i];
		y[0] = seedsFlow->pData[2 * i + 1];
		int nbCnt = 1;

		// add neighbors
		int* nbIdx = neighbors.rowPtr(i);
		for (int n = 0; n < maxNb; n++){
			if (nbIdx[n] < 0){
				break;
			}

			x[nbCnt] = seedsFlow->pData[2 * nbIdx[n]];
			y[nbCnt] = seedsFlow->pData[2 * nbIdx[n] + 1];
			nbCnt++;
		}

		float circleR = MinimalCircle(x, y, nbCnt);
		outRadius[i] = circleR;
	}
}

double CPM::dist(Point a, Point b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// intersection between two lines
CPM::Point CPM::intersection(Point u1, Point u2, Point v1, Point v2)
{
	Point ans = u1;
	double t = ((u1.x - v1.x) * (v1.y - v2.y) - (u1.y - v1.y) * (v1.x - v2.x)) /
		((u1.x - u2.x) * (v1.y - v2.y) - (u1.y - u2.y) * (v1.x - v2.x));
	ans.x += (u2.x - u1.x) * t;
	ans.y += (u2.y - u1.y) * t;
	return ans;
}

// circle center containing a triangular
CPM::Point CPM::circumcenter(Point a, Point b, Point c)
{
	Point ua, ub, va, vb;
	ua.x = (a.x + b.x) / 2;
	ua.y = (a.y + b.y) / 2;
	ub.x = ua.x - a.y + b.y;
	ub.y = ua.y + a.x - b.x;
	va.x = (a.x + c.x) / 2;
	va.y = (a.y + c.y) / 2;
	vb.x = va.x - a.y + c.y;
	vb.y = va.y + a.x - c.x;
	return intersection(ua, ub, va, vb);
}

float CPM::MinimalCircle(float* x, float*y, int n, float* centerX, float* centerY)
{
	static double eps = 1e-6;

	// prepare data
	Point p[20];
	assert(n < 20);
	for (int i = 0; i < n; i++){
		p[i].x = x[i];
		p[i].y = y[i];
	}
	// center and radius of the circle
	Point o;
	double r;

	int i, j, k;
	o = p[0];
	r = 0;
	for (i = 1; i < n; i++)
	{
		if (dist(p[i], o) - r > eps)
		{
			o = p[i];
			r = 0;

			for (j = 0; j < i; j++)
			{
				if (dist(p[j], o) - r > eps)
				{
					o.x = (p[i].x + p[j].x) / 2.0;
					o.y = (p[i].y + p[j].y) / 2.0;

					r = dist(o, p[j]);

					for (k = 0; k < j; k++)
					{
						if (dist(o, p[k]) - r > eps)
						{
							o = circumcenter(p[i], p[j], p[k]);
							r = dist(o, p[k]);
						}
					}
				}
			}
		}
	}

	if (centerX){
		*centerX = o.x;
	}
	if (centerY){
		*centerY = o.y;
	}
	return r;
}


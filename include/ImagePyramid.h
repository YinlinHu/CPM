#ifndef _GaussianPyramid_h
#define _GaussianPyramid_h

#include "Image.h"

template <class T>
class ImagePyramid
{
private:
	Image<T>* ImPyramid;
	int nLevels;
	float fRatio;
public:
	ImagePyramid(void){ ImPyramid = NULL; };
	~ImagePyramid(void){if(ImPyramid != NULL) delete[]ImPyramid;};
	inline Image<T>& operator[](int level) { return ImPyramid[level]; };
	void ConstructPyramid(const FImage& image, float ratio = 0.8, int minWidth = 30);
	void ConstructPyramidLevels(const FImage& image, float ratio = 0.8, int _nLevels = 2);
	void displayTop(const char* filename){ ImPyramid[nLevels - 1].imwrite(filename); };
	inline int nlevels() const {return nLevels;};
	inline float ratio() const { return fRatio; };
};

typedef ImagePyramid<float> FImagePyramid;

//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the slow way
//---------------------------------------------------------------------------------------
/*void GaussianPyramid::ConstructPyramid(const DImage &image, float ratio, int minWidth)
{
// the ratio cannot be arbitrary numbers
if(ratio>0.98 || ratio<0.4)
ratio=0.75;
// first decide how many levels
nLevels=log((float)minWidth/image.width())/log(ratio);
if(ImPyramid!=NULL)
delete []ImPyramid;
ImPyramid=new DImage[nLevels];
ImPyramid[0].copyData(image);
float baseSigma=(1/ratio-1);
for(int i=1;i<nLevels;i++)
{
DImage foo;
float sigma=baseSigma*i;
image.GaussianSmoothing(foo,sigma,sigma*2.5);
foo.imresize(ImPyramid[i],pow(ratio,i));
}
}//*/

//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the fast way
//---------------------------------------------------------------------------------------
template <class T>
void ImagePyramid<T>::ConstructPyramid(const FImage& image, float ratio /*= 0.8*/, int minWidth /*= 30*/)
{
	// the ratio cannot be arbitrary numbers
	if (ratio>0.98 || ratio<0.4)
		ratio = 0.75;
	// first decide how many levels
	nLevels = log((float)minWidth / image.width()) / log(ratio);
	fRatio = ratio;
	if (ImPyramid != NULL)
		delete[]ImPyramid;
	ImPyramid = new FImage[nLevels];
	ImPyramid[0].copyData(image);
	float baseSigma = (1 / ratio - 1);
	int n = log(0.25) / log(ratio);
	float nSigma = baseSigma*n;
	for (int i = 1; i<nLevels; i++)
	{
		FImage foo;
		if (i <= n)
		{
			float sigma = baseSigma*i;
			image.GaussianSmoothing(foo, sigma, sigma * 3);
			foo.imresize(ImPyramid[i], pow(ratio, i));
		}
		else
		{
			ImPyramid[i - n].GaussianSmoothing(foo, nSigma, nSigma * 3);
			float rate = (float)pow(ratio, i)*image.width() / foo.width();
			foo.imresize(ImPyramid[i], rate);
		}
	}
}

template <class T>
void ImagePyramid<T>::ConstructPyramidLevels(const FImage& image, float ratio /*= 0.8*/, int _nLevels /*= 2*/)
{
	// the ratio cannot be arbitrary numbers
	if (ratio>0.98 || ratio<0.4)
		ratio = 0.75;
	nLevels = _nLevels;
	fRatio = ratio;
	if (ImPyramid != NULL)
		delete[]ImPyramid;
	ImPyramid = new FImage[nLevels];
	ImPyramid[0].copyData(image);
	float baseSigma = (1 / ratio - 1);
	int n = log(0.25) / log(ratio);
	float nSigma = baseSigma*n;
	for (int i = 1; i<nLevels; i++)
	{
		FImage foo;
		if (i <= n)
		{
			float sigma = baseSigma*i;
			image.GaussianSmoothing(foo, sigma, sigma * 3);
			foo.imresize(ImPyramid[i], pow(ratio, i));
		}
		else
		{
			ImPyramid[i - n].GaussianSmoothing(foo, nSigma, nSigma * 3);
			float rate = (float)pow(ratio, i)*image.width() / foo.width();
			foo.imresize(ImPyramid[i], rate);
		}
	}
}

#endif
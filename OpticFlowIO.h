// flowIO.h

#ifndef _OpticFlowIO_H
#define _OpticFlowIO_H

#include <math.h>
#include <memory.h>
#include <string.h>

#include <stdint.h>
#include "opencv2/opencv.hpp" // for KITTI

// read and write our simple .flo flow file format

// ".flo" file format used for optical flow evaluation
//
// Stores 2-band float image for horizontal (u) and vertical (v) flow components.
// Floats are stored in little-endian order.
// A flow value is considered "unknown" if either |u| or |v| is greater than 1e9.
//
//  bytes  contents
//
//  0-3     tag: "PIEH" in ASCII, which in little endian happens to be the float 202021.25
//          (just a sanity check that floats are represented correctly)
//  4-7     width as an integer
//  8-11    height as an integer
//  12-end  data (width*height*2*4 bytes total)
//          the float values for u and v, interleaved, in row order, i.e.,
//          u[row0,col0], v[row0,col0], u[row0,col1], v[row0,col1], ...
//

// value to use to represent unknown flow
#define UNKNOWN_FLOW 1e10

typedef struct 
{
	double aee; // Average Endpoint Error
	double aae; // Average Angular Error
}FlowErr;

class OpticFlowIO
{
public:
	// return whether flow vector is unknown
	template <class T>
	static bool unknown_flow(T u, T v);
	template <class T>
	static bool unknown_flow(T *f);

	// read a flow file into 2-band image
	template <class T>
	static int ReadFlowFile(T* U, T* V, int* w, int* h, const char* filename);

	// write a 2-band image into flow file 
	template <class T>
	static int WriteFlowFile(T* U, T* V, int w, int h, const char* filename);

	// read a KITTI flow file into 2-band image
	template <class T>
	static int ReadKittiFlowFile(T* U, T* V, int* w, int* h, const char* filename);

	// write a 2-band image into KITTI flow file 
	template <class T>
	static int WriteKittiFlowFile(T* U, T* V, int w, int h, const char* filename);

	// render the motion to a 4-band BGRA color image
	template <class T>
	static double MotionToColor(unsigned char* fillPix, T* U, T* V, int w, int h, float range = -1);

	template <class T>
	static float ShowFlow(const char* winname, T* U, T* V, int w, int h, float range = -1, int waittime = 1);
	template <class T>
	static void SaveFlowAsImage(const char* imgName, T* U, T* V, int w, int h, float range = -1);

	template <class T>
	static float ErrorImage(unsigned char* fillPix, T* u1, T* v1, T* u2, T* v2, int w, int h);
	template <class T>
	static float ErrorImage(unsigned char* fillPix, T* u1, T* v1, char* gtName, int w, int h);
	template <class T>
	static float ShowErrorImage(const char* winname, T* U, T* V, char* gtName, int w, int h, int waittime = 1);
	template <class T>
	static float SaveErrorImage(const char* imgName, T* U, T* V, char* gtName, int w, int h);

	template <class T1, class T2>
	static FlowErr CalcFlowError(T1* u1, T1* v1, T2* u2, T2*v2, int w, int h);

private:
	// first four bytes, should be the same in little endian
	#define TAG_FLOAT 202021.25  // check for this when READING the file
	#define TAG_STRING "PIEH"    // use this when WRITING the file

	#define M_PI       3.14159265358979323846

	// the "official" threshold - if the absolute value of either 
	// flow component is greater, it's considered unknown
	#define UNKNOWN_FLOW_THRESH 1e9

	#define NUM_BANDS 2

	// Color encoding of flow vectors
	// adapted from the color circle idea described at
	//   http://members.shaw.ca/quadibloc/other/colint.htm
	//
	// Daniel Scharstein, 4/2007
	// added tick marks and out-of-range coding 6/05/07

	#define MAXWHEELCOLS 60
	template <class T>
	static void setcols(T* colorwheel, int r, int g, int b, int k);
	template <class T>
	static int makecolorwheel(T* colorwheel);
	template <class T>
	static void computeColor(double fx, double fy, unsigned char *pix, T* colorwheel, int ncols);
};

template <class T>
int OpticFlowIO::ReadKittiFlowFile(T* U, T* V, int* w, int* h, const char* filename)
{
	if (filename == NULL){
		printf("ReadKittiFlowFile: empty filename\n");
		return -1;
	}

	const char *dot = strrchr(filename, '.');
	if (strcmp(dot, ".png") != 0){
		printf("ReadKittiFlowFile (%s): extension .png expected\n", filename);
		return -1;
	}

	IplImage* img = cvLoadImage(filename, CV_LOAD_IMAGE_COLOR | CV_LOAD_IMAGE_ANYDEPTH);
	if (img == NULL){
		printf("ReadKittiFlowFile: could not open %s\n", filename);
		return -1;
	}

	int width = img->width;
	int height = img->height;
	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			uint16_t* rowImgData = (uint16_t*)(img->imageData + i*img->widthStep);
			uint16_t validFlag = rowImgData[j*img->nChannels];
			if(validFlag > 0){
				U[i*width+j] = (rowImgData[j*img->nChannels + 2] - 32768.0f)/64.0f;
				V[i*width+j] = (rowImgData[j*img->nChannels + 1] - 32768.0f)/64.0f;
			}else{
				U[i*width+j] = UNKNOWN_FLOW;
				V[i*width+j] = UNKNOWN_FLOW;
			}
		}
	}

	*w = width;
	*h = height;
	cvReleaseImage(&img);
	return 0;
}

template <class T>
int OpticFlowIO::WriteKittiFlowFile(T* U, T* V, int w, int h, const char* filename)
{
	if (filename == NULL){
		printf("WriteKittiFlowFile: empty filename\n");
		return -1;
	}

	const char *dot = strrchr(filename, '.');
	if (dot == NULL){
		printf("WriteKittiFlowFile: extension required in filename '%s'\n", filename);
		return -1;
	}

	if (strcmp(dot, ".png") != 0){
		printf("WriteKittiFlowFile: filename '%s' should have extension '.png'\n", filename);
		return -1;
	}

	int width = w, height = h;

	IplImage* img = cvCreateImage(cvSize(w,h), IPL_DEPTH_16U, 3);
	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			double u,v;
			u = U[i*width+j];
			v = V[i*width+j];
			uint16_t* rowImgData = (uint16_t*)(img->imageData + i*img->widthStep);
			if(!unknown_flow(u,v)){
				rowImgData[j*img->nChannels + 2] = __max(__min(U[i*width+j]*64.0f+32768.0f, 65535), 0);
				rowImgData[j*img->nChannels + 1] = __max(__min(V[i*width+j]*64.0f+32768.0f, 65535), 0);
				rowImgData[j*img->nChannels] = 1;
			}else{
				rowImgData[j*img->nChannels + 2] = 0;
				rowImgData[j*img->nChannels + 1] = 0;
				rowImgData[j*img->nChannels] = 0;
			}
		}
	}

	const int params[2]={CV_IMWRITE_PNG_COMPRESSION, 1};
	cvSaveImage(filename, img, params); // slight lossy PNG
	cvReleaseImage(&img);
	return 0;
}

template <class T>
bool OpticFlowIO::unknown_flow(T u, T v)
{
	return (abs(u) > UNKNOWN_FLOW_THRESH) 
		|| (abs(v) > UNKNOWN_FLOW_THRESH)
		|| u != u || v != v;	// isnan()
}

template <class T>
bool OpticFlowIO::unknown_flow(T *f)
{
	return unknown_flow(f[0], f[1]);
}

template <class T>
int OpticFlowIO::ReadFlowFile(T* U, T* V, int* w, int* h, const char* filename)
{
	if (filename == NULL){
		printf("ReadFlowFile: empty filename\n");
		return -1;
	}

	const char *dot = strrchr(filename, '.');
	if (strcmp(dot, ".flo") != 0){
		printf("ReadFlowFile (%s): extension .flo expected\n", filename);
		return -1;
	}

	FILE *stream = fopen(filename, "rb");
	if (stream == 0){
		printf("ReadFlowFile: could not open %s\n", filename);
		return -1;
	}

	int width, height;
	float tag;

	if ((int)fread(&tag,    sizeof(float), 1, stream) != 1 
		||(int)fread(&width,  sizeof(int),   1, stream) != 1 
		||(int)fread(&height, sizeof(int),   1, stream) != 1)
	{
		printf("ReadFlowFile: problem reading file %s\n", filename);
		return -1;
	}

	if (tag != TAG_FLOAT) // simple test for correct endian-ness
	{
		printf("ReadFlowFile(%s): wrong tag (possibly due to big-endian machine?)\n", filename);
		return -1;
	}

	// another sanity check to see that integers were read correctly (99999 should do the trick...)
	if (width < 1 || width > 99999){
		printf("ReadFlowFile(%s): illegal width %d\n", filename, width);
		return -1;
	}

	if (height < 1 || height > 99999){
		printf("ReadFlowFile(%s): illegal height %d\n", filename, height);
		return -1;
	}

	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			float tmp[NUM_BANDS];
			if ((int)fread(tmp, sizeof(float), NUM_BANDS, stream) != NUM_BANDS){
				printf("ReadFlowFile(%s): file is too short\n", filename);
				return -1;
			}
			U[i*width+j] = tmp[0];
			V[i*width+j] = tmp[1];
		}
	}

	if (fgetc(stream) != EOF){
		printf("ReadFlowFile(%s): file is too long\n", filename);
		return -1;
	}

	*w = width;
	*h = height;

	fclose(stream);
	return 0;
}

template <class T>
int OpticFlowIO::WriteFlowFile(T* U, T* V, int w, int h, const char* filename)
{
	if (filename == NULL){
		printf("WriteFlowFile: empty filename\n");
		return -1;
	}

	const char *dot = strrchr(filename, '.');
	if (dot == NULL){
		printf("WriteFlowFile: extension required in filename '%s'\n", filename);
		return -1;
	}

	if (strcmp(dot, ".flo") != 0){
		printf("WriteFlowFile: filename '%s' should have extension '.flo'\n", filename);
		return -1;
	}

	int width = w, height = h;

	FILE *stream = fopen(filename, "wb");
	if (stream == 0){
		printf("WriteFlowFile: could not open %s\n", filename);
		return -1;
	}

	// write the header
	fprintf(stream, TAG_STRING);
	if ((int)fwrite(&width,  sizeof(int),   1, stream) != 1 
		||(int)fwrite(&height, sizeof(int),   1, stream) != 1)
	{
		printf("WriteFlowFile(%s): problem writing header\n", filename);
		return -1;
	}

	// write the rows
	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			float tmp[NUM_BANDS];
			tmp[0] = U[i*width+j];
			tmp[1] = V[i*width+j];
			if ((int)fwrite(tmp, sizeof(float), NUM_BANDS, stream) != NUM_BANDS){
				printf("WriteFlowFile(%s): problem writing data\n", filename);
				return -1;
			}
		}
	}

	fclose(stream);
	return 0;
}

template <class T>
double OpticFlowIO::MotionToColor(unsigned char* fillPix, T* U, T* V, int w, int h, float range /*= -1*/)
{
	// determine motion range:
	double maxrad;

	if (range > 0) {
		maxrad = range;
	}else{	// obtain the motion range according to the max flow
		double maxu = -999, maxv = -999;
		double minu = 999, minv = 999;
		maxrad = -1;
		for (int i = 0; i < h; i++){
			for (int j = 0; j < w; j++){
				double u = U[i*w + j];
				double v = V[i*w + j];
				if (unknown_flow(u, v))
					continue;
				maxu = __max(maxu, u);
				maxv = __max(maxv, v);
				minu = __min(minu, u);
				minv = __min(minv, v);
				double rad = sqrt(u * u + v * v);
				maxrad = __max(maxrad, rad);
			}
		}
		if (maxrad == 0) // if flow == 0 everywhere
			maxrad = 1;
	}

	//printf("max motion: %.2f  motion range: u = [%.2f,%.2f];  v = [%.2f,%.2f]\n",
	//	maxrad, minu, maxu, minv, maxv);

	int colorwheel[MAXWHEELCOLS*3];
	int ncols = makecolorwheel(colorwheel);

	for(int i=0; i<h; i++){
		for(int j=0; j<w; j++){
			int idx = i*w+j;
			double u = U[idx];
			double v = V[idx];
			if (unknown_flow(u, v)){
				memset(fillPix+idx*4, 0, 4);
				fillPix[idx*4 + 3] = 0xff; // alpha channel, only for alignment
			}else{
				double dx = __min(__max(u / maxrad, -1), 1);
				double dy = __min(__max(v / maxrad, -1), 1);
				computeColor(dx, dy, (unsigned char*)(fillPix + idx * 4), colorwheel, ncols);
			}
		}
	}

	return maxrad;
}

template <class T>
float OpticFlowIO::ShowFlow(const char* winname, T* U, T* V, int w, int h, 
	float range /*= -1*/, int waittime /*= 1*/)
{
	cv::Mat img(h, w, CV_8UC4);
	float maxFlow = OpticFlowIO::MotionToColor(img.data, U, V, w, h, range);

#if 0
	// get corner color
	int x = 10, y = 20;
	unsigned char color[4];
	unsigned char* pSrc = img.data + y*img.step + x * 4;
	color[0] = 255 - pSrc[0];
	color[1] = 255 - pSrc[1];
	color[2] = 255 - pSrc[2];
	char info[256];
	sprintf(info, "max: %.1f", maxFlow);
	cv::putText(img, info, cvPoint(x, y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cvScalar(color[0], color[1], color[2]));
#endif

	cv::imshow(winname, img);
	cv::waitKey(waittime);

	return maxFlow;
}

template <class T>
void OpticFlowIO::SaveFlowAsImage(const char* imgName, T* U, T* V, int w, int h, float range /*= -1*/)
{
	cv::Mat img(h, w, CV_8UC4);
	float maxFlow = OpticFlowIO::MotionToColor(img.data, U, V, w, h, range);

#if 1
	// get corner color
	int x = 10, y = 20;
	unsigned char color[3];
	unsigned char* pSrc = img.data + y*img.step + x * 4;
	color[0] = 255 - pSrc[0];
	color[1] = 255 - pSrc[1];
	color[2] = 255 - pSrc[2];
	char info[256];
	sprintf(info, "max: %.1f", maxFlow);
	cv::putText(img, info, cvPoint(x, y), CV_FONT_HERSHEY_SIMPLEX, 0.5, cvScalar(color[0], color[1], color[2]));
#endif

	cv::imwrite(imgName, img);
}

template <class T>
float OpticFlowIO::ErrorImage(unsigned char* fillPix, T* u1, T* v1, T* u2, T* v2, int w, int h)
{
	unsigned char pix[4];

//#define LOG_COLOR
#ifdef LOG_COLOR
	float LC[10][5] =
	{ { 0, 0.0625, 49, 54, 149 },
	{ 0.0625, 0.125, 69, 117, 180 },
	{ 0.125, 0.25, 116, 173, 209 },
	{ 0.25, 0.5, 171, 217, 233 },
	{ 0.5, 1, 224, 243, 248 },
	{ 1, 2, 254, 224, 144 },
	{ 2, 4, 253, 174, 97 },
	{ 4, 8, 244, 109, 67 },
	{ 8, 16, 215, 48, 39 },
	{ 16, 1000000000.0, 165, 0, 38 } };
#endif

	int totalCnt = 0;
	int validCnt = 0;
	for (int i = 0; i < h; i++){
		for (int j = 0; j < w; j++){
			int idx = i*w + j;

			if (unknown_flow(u1[idx], v1[idx]) || unknown_flow(u2[idx], v2[idx])){
				// red: occlusion
				pix[0] = 0; pix[1] = 0; pix[2] = 255;
				pix[2] = 0; // TODO
				pix[3] = 0xFF; // only for alignment
				memcpy(fillPix + idx * 4, pix, 4);
				continue;
			}

			float endPtErr = sqrt(pow(u1[idx] - u2[idx], 2) + pow(v1[idx] - v2[idx], 2));

#ifdef LOG_COLOR
			float f_err = endPtErr;
			float f_mag = sqrt(u2[idx] * u2[idx] + v2[idx] * v2[idx]);
			float n_err = std::min(f_err / 3.0, 20.0*f_err / f_mag);
			for (int i = 0; i < 10; i++) {
				if (n_err >= LC[i][0] && n_err < LC[i][1]) {
					pix[3] = 0xFF; // only for alignment
					pix[2] = (uint8_t)LC[i][2];
					pix[1] = (uint8_t)LC[i][3];
					pix[0] = (uint8_t)LC[i][4];
				}
			}
			if (unknown_flow(u1[idx], v1[idx]) || unknown_flow(u2[idx], v2[idx])) {
				pix[2] *= 0.5;
				pix[1] *= 0.5;
				pix[0] *= 0.5;
			}
#else
			float v = __min(endPtErr, 5.0) / 5.0;
			//float v = __min(endPtErr, 3.0) / 3.0;
			pix[0] = v * 255; pix[1] = v * 255; pix[2] = v * 255;
			pix[3] = 0xFF; // only for alignment

#if 1
			if (endPtErr > 3.0){ // red
				pix[0] = pix[1] = 0;
				pix[2] = 255;
			}
#endif
#endif
			if (endPtErr < 3.0){
				validCnt++;
			}
			memcpy(fillPix + idx * 4, pix, 4);
			totalCnt++;
		}
	}
	return 1. - (float)validCnt / totalCnt;
}

template <class T>
float OpticFlowIO::ErrorImage(unsigned char* fillPix, T* u1, T* v1, char* gtName, int w, int h)
{
	int gtw, gth;
	T* u2 = new T[w*h];
	T* v2 = new T[w*h];
	ReadFlowFile(u2, v2, &gtw, &gth, gtName);
	assert(w == gtw&&h == gth);

	float r = ErrorImage(fillPix, u1, v1, u2, v2, w, h);
	delete[] u2;
	delete[] v2;
	return r;
}

template <class T>
float OpticFlowIO::ShowErrorImage(const char* winname, T* U, T* V, char* gtName, int w, int h, int waittime /*= 1*/)
{
	cv::Mat img(h, w, CV_8UC4);
	float r = OpticFlowIO::ErrorImage(img.data, U, V, gtName, w, h);

	cv::imshow(winname, img);
	cv::waitKey(waittime);

	return r;
}

template <class T>
float OpticFlowIO::SaveErrorImage(const char* imgName, T* U, T* V, char* gtName, int w, int h)
{
	cv::Mat img(h, w, CV_8UC4);
	float r = OpticFlowIO::ErrorImage(img.data, U, V, gtName, w, h);

	cv::imwrite(imgName, img);
	return r;
}

template <class T1, class T2>
FlowErr OpticFlowIO::CalcFlowError(T1* u1, T1* v1, T2* u2, T2*v2, int w, int h)
{
	FlowErr stat;
	memset(&stat, 0, sizeof(stat));

	double endPtErr = 0;
	double angErr = 0;
	int n = 0;
	for(int i=0; i<h; i++){
		for(int j=0; j<w; j++){
			int idx = i*w+j;
			if(unknown_flow(u1[idx], v1[idx]) || unknown_flow(u2[idx], v2[idx])){
				continue;
			}
			endPtErr += sqrt(pow(u1[idx]-u2[idx],2)+pow(v1[idx]-v2[idx],2));

			double tmp = (1.0+u1[idx]*u2[idx]+v1[idx]*v2[idx])
				/(sqrt(1.0+pow(u1[idx],2)+pow(v1[idx],2))*sqrt(1.0+pow(u2[idx],2)+pow(v2[idx],2)));
			
			if(tmp < -1.0)	tmp = -1.0;
			if(tmp > 1.0)	tmp = 1.0;

			angErr += acos(tmp);
			n++;
		}
	}

	stat.aae = (angErr/n) * 180/M_PI;
	stat.aee = endPtErr/n;

	return stat;
}

template <class T>
void OpticFlowIO::setcols(T* colorwheel, int r, int g, int b, int k)
{
	colorwheel[k*3+0] = r;
	colorwheel[k*3+1] = g;
	colorwheel[k*3+2] = b;
}

template <class T>
int OpticFlowIO::makecolorwheel(T* colorwheel)
{
	// relative lengths of color transitions:
	// these are chosen based on perceptual similarity
	// (e.g. one can distinguish more shades between red and yellow 
	//  than between yellow and green)
	int RY = 15;
	int YG = 6;
	int GC = 4;
	int CB = 11;
	int BM = 13;
	int MR = 6;
	int ncols = RY + YG + GC + CB + BM + MR;
	//printf("ncols = %d\n", ncols);
	if (ncols > MAXWHEELCOLS){
		printf("Too Many Columns in ColorWheel!\n");
		//exit(1);
	}
	int i;
	int k = 0;
	for (i = 0; i < RY; i++) setcols(colorwheel, 255,	   255*i/RY,	 0,	       k++);
	for (i = 0; i < YG; i++) setcols(colorwheel, 255-255*i/YG, 255,		 0,	       k++);
	for (i = 0; i < GC; i++) setcols(colorwheel, 0,		   255,		 255*i/GC,     k++);
	for (i = 0; i < CB; i++) setcols(colorwheel, 0,		   255-255*i/CB, 255,	       k++);
	for (i = 0; i < BM; i++) setcols(colorwheel, 255*i/BM,	   0,		 255,	       k++);
	for (i = 0; i < MR; i++) setcols(colorwheel, 255,	   0,		 255-255*i/MR, k++);

	return ncols;
}

template <class T>
void OpticFlowIO::computeColor(double fx, double fy, unsigned char *pix, T* colorwheel, int ncols)
{
	double rad = sqrt(fx * fx + fy * fy);
	double a = atan2(-fy, -fx) / M_PI;
	double fk = (a + 1.0) / 2.0 * (ncols-1);
	int k0 = (int)fk;
	int k1 = (k0 + 1) % ncols;
	double f = fk - k0;
	//f = 0; // uncomment to see original color wheel
	for (int b = 0; b < 3; b++) {
		double col0 = colorwheel[k0*3+b] / 255.0;
		double col1 = colorwheel[k1*3+b] / 255.0;
		double col = (1 - f) * col0 + f * col1;
		if (rad <= 1)
			col = 1 - rad * (1 - col); // increase saturation with radius
		else
			col *= .75; // out of range
		pix[2 - b] = (int)(255.0 * col);
	}
	pix[3] = 0xff; // alpha channel, only for alignment
}

#endif //_OpticFlowIO_H

#include "CPM.h"
#include "OpticFlowIO.h"

// draw each match as a 3x3 color block
void Match2Flow(FImage& inMat, FImage& ou, FImage& ov, int w, int h)
{
	if (!ou.matchDimension(w, h, 1)){
		ou.allocate(w, h, 1);
	}
	if (!ov.matchDimension(w, h, 1)){
		ov.allocate(w, h, 1);
	}
	ou.setValue(UNKNOWN_FLOW);
	ov.setValue(UNKNOWN_FLOW);
	int cnt = inMat.height();
	for (int i = 0; i < cnt; i++){
		float* p = inMat.rowPtr(i);
		float x = p[0];
		float y = p[1];
		float u = p[2] - p[0];
		float v = p[3] - p[1];
		for (int di = -1; di <= 1; di++){
			for (int dj = -1; dj <= 1; dj++){
				int tx = ImageProcessing::EnforceRange(x + dj, w);
				int ty = ImageProcessing::EnforceRange(y + di, h);
				ou[ty*w + tx] = u;
				ov[ty*w + tx] = v;
			}
		}
	}
}

void WriteMatches(const char *filename, FImage& inMat)
{
	int len = inMat.height();
	FILE *fid = fopen(filename, "w");
	for (int i = 0; i < len; i++){
		float x1 = inMat[4 * i + 0];
		float y1 = inMat[4 * i + 1];
		float x2 = inMat[4 * i + 2];
		float y2 = inMat[4 * i + 3];
		fprintf(fid, "%.0f %.0f %.0f %.0f\n", x1, y1, x2, y2);
		//fprintf(fid, "%.3f %.3f %.3f %.3f 1 100\n", x1, y1, x2, y2);
	}
	fclose(fid);
}

int main(int argc, char** argv)
{
	if (argc < 4){
		printf("USAGE: CPM image1 image2 outMatchText <step>\n");
		return -1;
	}

	FImage img1, img2;

	img1.imread(argv[1]);
	img2.imread(argv[2]);
	char* outMatName = argv[3];
	int step = 3;
	if (argc >= 5){
		step = atoi(argv[4]);
	}

	int w = img1.width();
	int h = img1.height();
	if (img2.width() != w || img2.height() != h){
		printf("CPM can only handle images with the same dimension!\n");
		return -1;
	}

	CTimer totalT;
	FImage matches;

	CPM cpm;
	cpm.SetStep(step);
	cpm.Matching(img1, img2, matches);

	totalT.toc("total time: ");

	FImage u, v;
	char tmpName[256];
	strcpy(tmpName, outMatName);
	strcat(tmpName, ".png");
	Match2Flow(matches, u, v, w, h);
	OpticFlowIO::SaveFlowAsImage(tmpName, u.pData, v.pData, w, h);

	WriteMatches(outMatName, matches);

	return 0;
}

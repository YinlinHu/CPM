
#include "CPM.h"

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

void main(int argc, char** argv)
{
	if (argc < 4){
		printf("USAGE: cpm.exe image1 image2 outMatchText\n");
		return;
	}

	FImage img1, img2;

	img1.imread(argv[1]);
	img2.imread(argv[2]);
	char* outMatName = argv[3];

	int w = img1.width();
	int h = img1.height();
	if (img2.width() != w || img2.height() != h){
		printf("CPM can only handle images with the same dimension!\n");
		return;
	}

	CTimer totalT;
	FImage matches;

	CPM cpm;
	cpm.Matching(img1, img2, matches);

	totalT.toc("total time: ");

	WriteMatches(outMatName, matches);
}

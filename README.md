# CPM
Code for

《Efficient Coarse-to-Fine PatchMatch for Large Displacement Optical Flow》 CVPR 2016

《Coarse-to-fine PatchMatch for dense correspondence》 T-CSVT

The program has been built and tested on Windows 7 and Ubuntu 16.04.

Note that, the code here is a little different from that when the paper accepted. Rather that SIFT-FLOW descriptor, we found DAISY descriptor is more suitable for practical usage (high density with a small sacrifice in accuracy). If you want to reproduce the result on KITTI and MPI accurately, please use the SIFT-Flow descriptor.

> USAGE: cpm img1Name img2Name outMatchName

Explanations:

The output of the program is a text file, which is in the format of "x1,y1,x2,y2" corresponding to one match per line.

huyinlin@gmail.com

Yinlin

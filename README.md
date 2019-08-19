# Introduction

This repository contains the code for the paper **Efficient Coarse-to-Fine PatchMatch for Large Displacement Optical Flow**. [Yinlin Hu](http://yinlinhu.github.io), Rui Song, and Yunsong Li. CVPR. 2016. [\[Paper\]](https://zpascal.net/cvpr2016/Hu_Efficient_Coarse-To-Fine_PatchMatch_CVPR_2016_paper.pdf)

# How to Use

It is assumed that the OpenCV has been installed correctly.

```
$ cmake .
$ make
```

Then, you can play with the examples:

```
$ bash demo.sh
```

# Notes
For Windows system support, you can explore 'CMakeLists_win.txt' and prebuilt binaries in the directory 'win32'.

# Citing

```
@inproceedings{hu2016cpm,
  title={Efficient Coarse-to-Fine PatchMatch for Large Displacement Optical Flow},
  author={Yinlin Hu and Rui Song and Yunsong Li},
  booktitle={CVPR},
  year={2016}
}
```
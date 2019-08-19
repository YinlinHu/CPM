CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(PROJ_NAME CPM)
PROJECT (${PROJ_NAME})

ADD_DEFINITIONS(-DWITH_SSE)

IF(MSVC)
INCLUDE_DIRECTORIES(
./
./include
C:/OpenCV/include
C:/OpenCV/include/opencv
C:/OpenCV/include/opencv2
)

LINK_DIRECTORIES(
C:/OpenCV/lib/x86/Debug
C:/OpenCV/lib/x86/Release
)
ELSE(MSVC)
FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(./include)
ENDIF(MSVC)

ADD_EXECUTABLE(${PROJ_NAME} main.cpp CPM.cpp)

IF(MSVC)
TARGET_LINK_LIBRARIES(${PROJ_NAME}
debug opencv_core310d.lib
debug opencv_highgui310d.lib
debug opencv_imgcodecs310d.lib
debug opencv_imgproc310d.lib
debug opencv_xfeatures2d310d.lib
optimized opencv_core310.lib
optimized opencv_highgui310.lib
optimized opencv_imgcodecs310.lib
optimized opencv_imgproc310.lib
optimized opencv_xfeatures2d310.lib)
ELSE(MSVC)
TARGET_LINK_LIBRARIES(${PROJ_NAME} ${OpenCV_LIBS})
ENDIF(MSVC)

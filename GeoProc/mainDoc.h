#ifndef MAINDOC_H
#define MAINDOC_H
//#include "surfprocess.h"
#include "NormalStruct.h"
#include "gradientFieldFilterRBF.h"

class mainDoc
{
public:
	mainDoc();
	~mainDoc();
	bool setTriangleGeometry(std::vector<Point>);	
	bool setTriangleTopology(std::vector<Triangle>);
	bool setTriangleTexture(std::vector<Point>);
	//bool setPointCloudsGeometry(  );
	//bool setPointCloudsTexture(  );
	bool addPositionGaussianNoise(); //加高斯噪声
	bool RBFilter();
	float calculateMeanLength(); //// 计算邻域的平均半径
protected:
	bool calculateTriangleNormal();
	bool calculatePointNormal();
	bool calculateDataBox();
	
public:
	int* m_triangleTopology;
	float* m_triangleGeometry;
	float* m_triangleGeometryWithNoise; // 加噪声后的几何位置数据
	float* m_triangleGeometryAfterFiltering; // 滤波后的几何位置数据
	float* m_triangleTexture;
	float* m_triangleNormal;
	float* m_triangleNormalWithNoise; // // 加噪声后的法向数据
	float* m_triangleNormalAfterFiltering; // 滤波后的法向数据

	gradientFieldFilterRBF * filter;

	int m_numOfVertex;
	int m_numOfTriangle;
	float box_x_min;
	float box_x_max;
	float box_y_min;
	float box_y_max;
	float box_z_min;
	float box_z_max;

	int m_itol; // 停止标准 1  2  3  4
	double m_tol; // 运行的最大误差
	int m_itmax; // 循环的最大次数
	int m_iter; //实际循环的次数
	double m_err; // 实际误差
	//int m_numOfPoints;
	int m_kNearest;

	double* m_originalNormals; // 原始法向
	double* m_filterNormals; // 滤波后的法向
	double* m_originalMajorDirection;//原始主方向（对应于最大的特征值）
	double* m_originalMinorDirection;//原始次方向（对应于第二大特征值）
	double* m_meanCurvature;//平均曲率
	//double* m_gaussCurvature; //高斯曲率
	double m_radius;//计算 m_mapKnearest时候的搜索半径
	double* m_loadVector; //载荷向量  AX=B中的B
	double m_curveThreshold;//最大曲率
	//double m_tol; // 运行的最大误差
protected:
	float m_noiseCoffecient; // 法向噪声系数

};
#endif
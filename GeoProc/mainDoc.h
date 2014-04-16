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
	bool addPositionGaussianNoise(); //�Ӹ�˹����
	bool RBFilter();
	float calculateMeanLength(); //// ���������ƽ���뾶
protected:
	bool calculateTriangleNormal();
	bool calculatePointNormal();
	bool calculateDataBox();
	
public:
	int* m_triangleTopology;
	float* m_triangleGeometry;
	float* m_triangleGeometryWithNoise; // ��������ļ���λ������
	float* m_triangleGeometryAfterFiltering; // �˲���ļ���λ������
	float* m_triangleTexture;
	float* m_triangleNormal;
	float* m_triangleNormalWithNoise; // // ��������ķ�������
	float* m_triangleNormalAfterFiltering; // �˲���ķ�������

	gradientFieldFilterRBF * filter;

	int m_numOfVertex;
	int m_numOfTriangle;
	float box_x_min;
	float box_x_max;
	float box_y_min;
	float box_y_max;
	float box_z_min;
	float box_z_max;

	int m_itol; // ֹͣ��׼ 1  2  3  4
	double m_tol; // ���е�������
	int m_itmax; // ѭ����������
	int m_iter; //ʵ��ѭ���Ĵ���
	double m_err; // ʵ�����
	//int m_numOfPoints;
	int m_kNearest;

	double* m_originalNormals; // ԭʼ����
	double* m_filterNormals; // �˲���ķ���
	double* m_originalMajorDirection;//ԭʼ�����򣨶�Ӧ����������ֵ��
	double* m_originalMinorDirection;//ԭʼ�η��򣨶�Ӧ�ڵڶ�������ֵ��
	double* m_meanCurvature;//ƽ������
	//double* m_gaussCurvature; //��˹����
	double m_radius;//���� m_mapKnearestʱ��������뾶
	double* m_loadVector; //�غ�����  AX=B�е�B
	double m_curveThreshold;//�������
	//double m_tol; // ���е�������
protected:
	float m_noiseCoffecient; // ��������ϵ��

};
#endif
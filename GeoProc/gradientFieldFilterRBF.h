#pragma once
#include "NormalStruct.h"

#if !defined(stiff_Matrix)
#define stiff_Matrix
struct StiffMatrix{
	std::vector<double> m_elementStiffMatrix; //��ŸնȾ����Ԫ��
}; // �����նȾ���Ľṹ
#endif

class gradientFieldFilterRBF
{
public:
	gradientFieldFilterRBF(void);
	~gradientFieldFilterRBF(void);


////////////////////////////////
// ��������
public:
	static const int m_numOfIntegralPoint;
	static const double m_ordinateAndWeightCoff[120];
	static const double m_coffRBFunction[5];

protected:
	double* m_originalPointSet;
	double* m_originalColors;

	int m_numOfPoints;
	int m_kNearest;

	double* m_originalNormals; // ԭʼ����
	double* m_filterNormals; // �˲���ķ���
	double* m_originalMajorDirection;//ԭʼ�����򣨶�Ӧ����������ֵ��
	double* m_originalMinorDirection;//ԭʼ�η��򣨶�Ӧ�ڵڶ�������ֵ��
	double* m_meanCurvature;//ƽ������
	//double* m_gaussCurvature; //��˹����
	double m_radius;//���� m_mapKnearestʱ��������뾶
	double* m_loadVector; //�غ�����  AX=B�е�B

	StiffMatrix* m_stiffMatrix;     //      �նȾ��������ݶȻ�������ֵ
	StiffMatrix* m_massMatrix;      //      �նȾ������ɱ�����������ֵ

protected: ///// �������йصĲ���
	int m_numOfNeighbors; //���������������ĵ㣩
	double m_localRBFradius;  //�ֲ�����������뾶
	double* m_localOrdinate; //�������ľֲ�����
	double* m_matrixNeighbor; // ��ֵ����
	double* m_matrixNeighborInv; // ��ֵ���������

protected: // ������йصĲ���
	double* m_localIntegralPoints; // ���ֵ������
	double m_localIntegralRadius; // �ֲ����ְ뾶
	double* m_testFunction; // ���Ժ�����ÿ���������ֵ
	double* m_RBFfunction; // ÿһ��Ϊÿһ�������ľ����������ÿһ�����ֵ�ĺ���ֵ
	double* m_RBFfunctionDerOne; //����
	double* m_RBFfunctionDerTwo; //����
	double* m_testFunctionDerOne;  // ����
	double* m_testFunctionDerTwo;  //����
	double* ��_scalarProductVector; // �ڲ����������пռ��У�n��RBF��������Ժ������ֵ�����
	double* m_gradientProductVector; // �ڲ����������пռ��У�n��RBF�������ݶȷֱ�����Ժ����ݶȻ��ֵ�����

protected:  // �κ������䵼��
	double* m_shapeFunction; //���������κ����ڸ��������ֵ
	double* m_shapeFunctionDerOne; // ����
	double* m_shapeFunctionDerTwo;// ����

protected: // ������ֲ�����ת���Ĳ���
	double* m_determinentMetric;//   �ڲ����㴦���α�׼�Ĵ�������ʽ
	double* m_sqrtDetMetric; // �ڲ����㴦���α�׼�Ĵ�������ʽ�Ŀ�����
	double* m_metricMatrix;//�ڲ����㴦���α�׼�ľ���
	double* m_metricMatrixInve;////�ڲ����㴦���α�׼�ľ������

public: // ���Է������һЩϵ��
	int m_itol; // ֹͣ��׼ 1  2  3  4
	double m_tol; // ���е�������
	int m_itmax; // ѭ����������
	int m_iter; //ʵ��ѭ���Ĵ���
	double m_err; // ʵ�����

public: // 
	double* m_resultPointSet;
	long nearestTime;
	long filterTime;

	std::map<int, KnearestField> m_mapKnearest;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///// ��������
public:
	void GetGradientFieldFilter(int numOfPointSet,float* pointSet,int kNearest, double radius,
		 double ebusainu, int stopN, int maxIter,double curveThreshold);

	/*
	�������ܣ�  �ӿں���������㼯�˲�
	����˵����
	int numOfPointSet                        �������
	double* pointSet							ԭʼ�㼯
	int kNearest                             k���������
	double radius                            k���������İ뾶
	double ebusainu                          �������Է�����ʱ���ֹͣ��׼
	int stopN                                �������Է�����ʱ��ѭ������
	*/

	void DeleteGradientFieldFilter();
	/*
	 *  �������ܣ� �ͷű����ڴ�
	 */

	void ComputeMapKnearest();
	/*
	 *  �������ܣ� ����ÿһ�����k����
	 */

	void ComputeNormal();
	/*
	*	�������ܣ�����ÿһ����ķ���
	*/

    	void ComputeNearestParameter(int indexPointSet);
	/*
	 *	�������ܣ� ������һ�㴦�Ļ��ְ뾶���������������ľֲ����꣬
	 *  ����˵����
	 *  int indexPointSet        ������	
	 *  double m_localRBFradius;            ����������뾶
	 *  double* m_localOrdinate             �������ľֲ�����
	 *  double* m_localIntegralPoint        ���ֵ������
	 *  double  m_localIntegralRadius       ������İ뾶
  	 */

	void ComputeTestFunction(int indexPointSet);
	/*
	 *	������Ժ����ڻ��ֵ�ĺ���ֵ���䵼��ֵ
	 */

	void ComputeMatrixNeighbor(int indexPointSet);
	/*
  	 *	�������ܣ������ֵ����������
	 *  ����˵����
	 *  double* m_matrixNeighbor; ��ֵ����
	 *  double* m_matrixNeighborInv; ��ֵ���������  
	 */

	void ComputeMatrixRelatedSampling(int indexPointSet);
	/*
	 *	�������ܣ� ������������йصľ��� 
	 *  ����˵���� 
	 *  int indexPointSet        ������
	 *  double* m_RBFfunction;        ÿһ��Ϊ��ÿһ�������ľ����������ÿһ�����ֵ�ĺ���ֵ
	 */
	
	void ComputeManifoldMetric(int indexPointSet);
	/*
	 *	�������ܣ� �������α�׼�Ĵ�������ʽ
	 *  ����˵����
	 *  int indexPointSet        ������
	 *  double* m_determinentMetric   �ڲ����㴦���α�׼�Ĵ�������ʽ
	 */
	void ComputeScalarProduct( int indexPointSet );
	/*
	 *  �������ܣ� ����RBF��������Ժ������ڻ�
	 *  ����˵����
	 *  int indexPointSet        ������
	 *  double*����_scalarProcuctVector  �ڲ����������пռ��У�n��RBF��������Ժ������ֵ�����
	 */

	void ComputeVectorProductVector( int indexPointSet );
	/*
	 * �������ܣ� ����RBF��������Ժ����ݶȵ��ڻ�
	 * ����˵����
	 * int indexPointSet        ������
	 * double* m_gradientProductVector  �ڲ����������пռ��У�n��RBF�������ݶȷֱ�����Ժ����ݶȻ��ֵ�����
	 */

	void ComputeStiffMatrix(int indexPointSet);
	/*
	 *	�������ܣ� ����նȾ���ͺ�������
	 *  ��ÿһ��������indexPointSet��Ϊ�նȾ������һ�У�Ϊ������������һ��ֵ
	 *  ����˵���� 
	 *  int indexPointSet                    ������
	 *  double* m_loadVector                  ��������
	 *  StiffMatrix* m_stiffMatrix           �նȾ��������ݶȻ�������ֵ
	 *  StiffMatrix* m_massMatrix            �նȾ������ɱ�����������ֵ
	 */

	void EstimateNormalDirection();
	/*
	 *	�������ܣ� ���Ʒ����׼ȷ����
	 */

	void CalculateLinearSystem();
	/*
  	 *	�������ܣ� �����Է���
	 */
	
	void  ComputeGradientShapeFunction();
	/*
	 *	�������ܣ� ����ÿһ�����ֵ��κ������ݶȺͲ��Ժ������ݶ�
	 *  ����˵����
	 *  double* m_gradientShapeFunction   �κ������ݶ�
	 *  double* m_gradientTestFunction    ���Ժ������ݶ�
	 */

	void CalculateLinearSystemGauss();
	/*
	*	�������ܣ� ��˹��Ԫ���ⷽ�̣�
	*/
	void AssembleStiffMatrix();
	//////////////////////////////////////////////////////////////////////////
	/*
	*	�����նȾ���
	*/
	void AdjustStiffMatrix();
	//////////////////////////////////////////////////////////////////////////
	/*
	*	pbcg ��ϡ�����Է�����
	*/
	//////////////////////////////////////////////////////////////////////////

	void linbcg(double* b,double* x,int itol,double tol,int itmax, int & iter, double &err);
	double snrm(double* sx,const int itol);
	void atimes(double* x,double* r,const int itrnsp);
	void asolve(double* b,double* x,const int itrnsp);






















};


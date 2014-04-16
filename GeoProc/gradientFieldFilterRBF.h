#pragma once
#include "NormalStruct.h"

#if !defined(stiff_Matrix)
#define stiff_Matrix
struct StiffMatrix{
	std::vector<double> m_elementStiffMatrix; //存放刚度矩阵的元素
}; // 保存点刚度矩阵的结构
#endif

class gradientFieldFilterRBF
{
public:
	gradientFieldFilterRBF(void);
	~gradientFieldFilterRBF(void);


////////////////////////////////
// 变量声明
public:
	static const int m_numOfIntegralPoint;
	static const double m_ordinateAndWeightCoff[120];
	static const double m_coffRBFunction[5];

protected:
	double* m_originalPointSet;
	double* m_originalColors;

	int m_numOfPoints;
	int m_kNearest;

	double* m_originalNormals; // 原始法向
	double* m_filterNormals; // 滤波后的法向
	double* m_originalMajorDirection;//原始主方向（对应于最大的特征值）
	double* m_originalMinorDirection;//原始次方向（对应于第二大特征值）
	double* m_meanCurvature;//平均曲率
	//double* m_gaussCurvature; //高斯曲率
	double m_radius;//计算 m_mapKnearest时候的搜索半径
	double* m_loadVector; //载荷向量  AX=B中的B

	StiffMatrix* m_stiffMatrix;     //      刚度矩阵中由梯度积产生的值
	StiffMatrix* m_massMatrix;      //      刚度矩阵中由标量积产生的值

protected: ///// 与邻域有关的参数
	int m_numOfNeighbors; //邻域数（包括中心点）
	double m_localRBFradius;  //局部径向基函数半径
	double* m_localOrdinate; //各邻域点的局部坐标
	double* m_matrixNeighbor; // 插值矩阵
	double* m_matrixNeighborInv; // 插值矩阵的逆阵

protected: // 与采样有关的参数
	double* m_localIntegralPoints; // 积分点的坐标
	double m_localIntegralRadius; // 局部积分半径
	double* m_testFunction; // 测试函数在每个采样点的值
	double* m_RBFfunction; // 每一行为每一个邻域点的径向基函数在每一个积分点的函数值
	double* m_RBFfunctionDerOne; //导数
	double* m_RBFfunctionDerTwo; //导数
	double* m_testFunctionDerOne;  // 导数
	double* m_testFunctionDerTwo;  //导数
	double* ｍ_scalarProductVector; // 在采样点所在切空间中，n个RBF函数与测试函数积分的向量
	double* m_gradientProductVector; // 在采样点所在切空间中，n个RBF函数的梯度分别与测试函数梯度积分的向量

protected:  // 形函数及其导数
	double* m_shapeFunction; //各邻域点的形函数在各采样点的值
	double* m_shapeFunctionDerOne; // 导数
	double* m_shapeFunctionDerTwo;// 导数

protected: // 流形与局部坐标转换的参数
	double* m_determinentMetric;//   在采样点处流形标准的代数行列式
	double* m_sqrtDetMetric; // 在采样点处流形标准的代数行列式的开根号
	double* m_metricMatrix;//在采样点处流形标准的矩阵
	double* m_metricMatrixInve;////在采样点处流形标准的矩阵的逆

public: // 线性方程组的一些系数
	int m_itol; // 停止标准 1  2  3  4
	double m_tol; // 运行的最大误差
	int m_itmax; // 循环的最大次数
	int m_iter; //实际循环的次数
	double m_err; // 实际误差

public: // 
	double* m_resultPointSet;
	long nearestTime;
	long filterTime;

	std::map<int, KnearestField> m_mapKnearest;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
///// 函数声明
public:
	void GetGradientFieldFilter(int numOfPointSet,float* pointSet,int kNearest, double radius,
		 double ebusainu, int stopN, int maxIter,double curveThreshold);

	/*
	函数功能：  接口函数，计算点集滤波
	变量说明：
	int numOfPointSet                        点的数量
	double* pointSet							原始点集
	int kNearest                             k邻域的数量
	double radius                            k邻域搜索的半径
	double ebusainu                          计算线性方程组时候的停止标准
	int stopN                                计算线性方程组时的循环次数
	*/

	void DeleteGradientFieldFilter();
	/*
	 *  函数功能： 释放变量内存
	 */

	void ComputeMapKnearest();
	/*
	 *  函数功能： 计算每一个点的k邻域
	 */

	void ComputeNormal();
	/*
	*	函数功能：计算每一个点的法向
	*/

    	void ComputeNearestParameter(int indexPointSet);
	/*
	 *	函数功能： 计算任一点处的积分半径，计算其邻域各点的局部坐标，
	 *  变量说明：
	 *  int indexPointSet        索引点	
	 *  double m_localRBFradius;            径向基函数半径
	 *  double* m_localOrdinate             各邻域点的局部坐标
	 *  double* m_localIntegralPoint        积分点的坐标
	 *  double  m_localIntegralRadius       积分域的半径
  	 */

	void ComputeTestFunction(int indexPointSet);
	/*
	 *	计算测试函数在积分点的函数值及其导数值
	 */

	void ComputeMatrixNeighbor(int indexPointSet);
	/*
  	 *	函数功能：计算插值矩阵及其逆阵
	 *  变量说明：
	 *  double* m_matrixNeighbor; 插值矩阵
	 *  double* m_matrixNeighborInv; 插值矩阵的逆阵  
	 */

	void ComputeMatrixRelatedSampling(int indexPointSet);
	/*
	 *	函数功能： 计算与采样点有关的矩阵 
	 *  变量说明： 
	 *  int indexPointSet        索引点
	 *  double* m_RBFfunction;        每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	 */
	
	void ComputeManifoldMetric(int indexPointSet);
	/*
	 *	函数功能： 计算流形标准的代数行列式
	 *  变量说明：
	 *  int indexPointSet        索引点
	 *  double* m_determinentMetric   在采样点处流形标准的代数行列式
	 */
	void ComputeScalarProduct( int indexPointSet );
	/*
	 *  函数功能： 计算RBF函数与测试函数的内积
	 *  变量说明：
	 *  int indexPointSet        索引点
	 *  double*　ｍ_scalarProcuctVector  在采样点所在切空间中，n个RBF函数与测试函数积分的向量
	 */

	void ComputeVectorProductVector( int indexPointSet );
	/*
	 * 函数功能： 计算RBF函数与测试函数梯度的内积
	 * 变量说明：
	 * int indexPointSet        索引点
	 * double* m_gradientProductVector  在采样点所在切空间中，n个RBF函数的梯度分别与测试函数梯度积分的向量
	 */

	void ComputeStiffMatrix(int indexPointSet);
	/*
	 *	函数功能： 计算刚度矩阵和荷载向量
	 *  在每一个索引点indexPointSet，为刚度矩阵计算一行，为荷载向量计算一个值
	 *  变量说明： 
	 *  int indexPointSet                    索引点
	 *  double* m_loadVector                  荷载向量
	 *  StiffMatrix* m_stiffMatrix           刚度矩阵中由梯度积产生的值
	 *  StiffMatrix* m_massMatrix            刚度矩阵中由标量积产生的值
	 */

	void EstimateNormalDirection();
	/*
	 *	函数功能： 估计法向的准确方向
	 */

	void CalculateLinearSystem();
	/*
  	 *	函数功能： 解线性方程
	 */
	
	void  ComputeGradientShapeFunction();
	/*
	 *	函数功能： 计算每一个积分点形函数的梯度和测试函数的梯度
	 *  变量说明：
	 *  double* m_gradientShapeFunction   形函数的梯度
	 *  double* m_gradientTestFunction    测试函数的梯度
	 */

	void CalculateLinearSystemGauss();
	/*
	*	函数功能： 高斯消元法解方程；
	*/
	void AssembleStiffMatrix();
	//////////////////////////////////////////////////////////////////////////
	/*
	*	调整刚度矩阵
	*/
	void AdjustStiffMatrix();
	//////////////////////////////////////////////////////////////////////////
	/*
	*	pbcg 求稀疏线性方程组
	*/
	//////////////////////////////////////////////////////////////////////////

	void linbcg(double* b,double* x,int itol,double tol,int itmax, int & iter, double &err);
	double snrm(double* sx,const int itol);
	void atimes(double* x,double* r,const int itrnsp);
	void asolve(double* b,double* x,const int itrnsp);






















};


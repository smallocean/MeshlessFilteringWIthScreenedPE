#include "gradientFieldFilterRBF.h"
#include <math.h>
#include "Matrix.h"
#include ".\mathlib\mathlib.h"

using namespace MATHLIB;
const int gradientFieldFilterRBF::m_numOfIntegralPoint=40;
const double gradientFieldFilterRBF::m_ordinateAndWeightCoff[120]={ 
//////////////////////////////////////////////////////////////////////////
/*  
*   第一列为半径为1时候的x坐标
*   第二列为半径为1时候的y坐标
*   第三列为两个一维高斯积分的加权值的乘积
*/
    -0.046546,-0.005837,0.023984,
	-0.037660,-0.027969,0.052688,
	-0.003759,-0.046759,0.074326,
	0.039334,-0.025562,0.085930,
	0.039334,0.025562,0.085930,
	-0.003759,0.046759,0.074326,
	-0.037660,0.027969,0.052688,
	-0.046546,0.005837,0.023984,
	-0.228972,-0.028714,0.048451,
	-0.185262,-0.137588,0.106438,
	-0.018490,-0.230023,0.150149,
	0.193496,-0.125745,0.173591,
	0.193496,0.125745,0.173591,
	-0.018490,0.230023,0.150149,
	-0.185262,0.137588,0.106438,
	-0.228972,0.028714,0.048451,
	-0.496114,-0.062215,0.057588,
	-0.401409,-0.298113,0.126510,
	-0.040063,-0.498392,0.178464,
	0.419249,-0.272453,0.206327,
	0.419249,0.272453,0.206327,
	-0.040063,0.498392,0.178464,
	-0.401409,0.298113,0.126510,
	-0.496114,0.062215,0.057588,
	-0.763256,-0.095716,0.048451,
	-0.617555,-0.458637,0.106438,
	-0.061636,-0.766761,0.150149,
	0.645001,-0.419161,0.173591,
	0.645001,0.419161,0.173591,
	-0.061636,0.766761,0.150149,
	-0.617555,0.458637,0.106438,
	-0.763256,0.095716,0.048451,
	-0.046546,-0.005837,0.023984,
	-0.037660,-0.027969,0.052688,
	-0.003759,-0.046759,0.074326,
	0.039334,-0.025562,0.085930,
	0.039334,0.025562,0.085930,
	-0.003759,0.046759,0.074326,
	-0.037660,0.027969,0.052688,
	-0.046546,0.005837,0.023984
};
const double gradientFieldFilterRBF::m_coffRBFunction[5]={8,40,48,25,5};


gradientFieldFilterRBF::gradientFieldFilterRBF(){

	/// 初始化所有数据
	m_originalPointSet = NULL;
	m_numOfPoints = NULL;
	m_numOfPoints = 0;
	m_kNearest = 0;
	m_originalNormals = NULL;
	m_radius = 0; // 计算 m_mapKnearest的搜索半径
	m_numOfNeighbors = 0; // 邻域数不包括中心点

	m_localRBFradius=0;//         局部径向基函数半径
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	m_matrixNeighbor=NULL;//            插值矩阵
	m_matrixNeighborInv=NULL;//         插值矩阵的逆阵
	m_localIntegralPoints=NULL;  //积分点的坐标
	m_localIntegralRadius=0;//      局部积分半径
	m_testFunction=NULL;//       测试函数在每一个采样点的值
	m_RBFfunction=NULL;//        每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式
	m_sqrtDetMetric = NULL;
	m_metricMatrix=NULL;//在采样点处流形标准的矩阵
	m_metricMatrixInve=NULL;////在采样点处流形标准的矩阵的逆
	m_resultPointSet=NULL;

	nearestTime=0;
	filterTime=0;
	m_RBFfunctionDerOne=NULL; 


	m_RBFfunctionDerTwo=NULL;
	m_shapeFunctionDerOne=NULL;
	m_shapeFunctionDerTwo=NULL;
	m_testFunctionDerOne=NULL;
	m_testFunctionDerTwo=NULL;
	m_stiffMatrix=NULL;
	m_massMatrix=NULL;
	ｍ_scalarProductVector = NULL;
	m_gradientProductVector = NULL;


	m_originalNormals=NULL;
	m_originalMajorDirection=NULL;
	m_originalMinorDirection=NULL;
	m_meanCurvature=NULL;
}

gradientFieldFilterRBF::~gradientFieldFilterRBF(){

	DeleteGradientFieldFilter();

}

void gradientFieldFilterRBF::GetGradientFieldFilter(int numOfPointSet, float *pointSet, int kNearest, double radius, double ebusainu, int stopN, int maxIter, double curveThreshold){

	/*
	 * 函数功能：  接口函数，计算点集滤波
	 * 变量说明：
	 * int numOfPointSet                        点的数量
	 * double* pointSet							原始点集
	 * int kNearest                             k邻域的数量
	 * double radius                            k邻域搜索的半径
	 * double ebusainu                          计算线性方程组时候的停止标准
	 * 	int stopN                                计算线性方程组时的循环次数
	 * maxIter   双共轭梯度下降法求解稀疏线性方程组时的最大迭代次数
	 */
	m_numOfPoints = numOfPointSet;
	m_kNearest = kNearest;
	m_radius = radius;
	m_itol = stopN;
	m_itmax = maxIter;
	m_tol = ebusainu;
	if ( m_originalPointSet != NULL ){
		delete[] m_originalPointSet;
		m_originalPointSet = NULL;
	}
	m_originalPointSet = new double[ m_numOfPoints * 3 ];

	if ( m_resultPointSet != NULL ){
		delete[] m_resultPointSet;
		m_resultPointSet = NULL;
	}
	m_resultPointSet = new double[ m_numOfPoints * 3 ];
	
	for( int i = 0; i < m_numOfPoints * 3; i++ ){
		m_originalPointSet[ i ] = ( double ) pointSet[ i ];
		m_resultPointSet[ i ] = ( double ) pointSet[ i ];
	}

	ComputeMapKnearest();
	ComputeNormal();

	if( m_stiffMatrix != NULL ){
		delete[] m_stiffMatrix;
		m_stiffMatrix = NULL;
	}
	m_stiffMatrix = new StiffMatrix[ m_numOfPoints ];

	if ( m_massMatrix != NULL ){
		delete[] m_massMatrix;
		m_massMatrix = NULL;
	}
	m_massMatrix = new StiffMatrix[ m_numOfPoints ];

	if ( m_loadVector != NULL ){
		delete[] m_loadVector;
		m_loadVector = NULL;
	}
	m_loadVector = new double[ m_numOfPoints ];

	for ( int i = 0; i < m_numOfPoints; i++ ){
		ComputeNearestParameter( i );
		ComputeTestFunction( i );
		ComputeMatrixNeighbor( i );
		ComputeMatrixRelatedSampling( i );
		ComputeManifoldMetric( i );
		ComputeScalarProduct( i );
		ComputeVectorProductVector( i );
		ComputeStiffMatrix( i );

		if ( m_localOrdinate != NULL ){
			delete[] m_localOrdinate;
			m_localOrdinate = NULL;
		}

		if ( m_determinentMetric != NULL ){
			delete[] m_determinentMetric;
			m_determinentMetric = NULL;
		}

		if ( m_metricMatrix != NULL ){
			delete[] m_metricMatrix;
			m_metricMatrix = NULL;
		}

		if ( m_metricMatrixInve != NULL ){
			delete[] m_metricMatrixInve;
			m_metricMatrixInve = NULL;
		}

		if ( m_testFunction != NULL ){
			delete[] m_testFunction;
			m_testFunction = NULL;
		}

		if ( m_testFunctionDerOne != NULL ){
			delete[] m_testFunctionDerOne;
			m_testFunctionDerOne = NULL;
		}

		if ( m_testFunctionDerTwo != NULL ){
			delete[] m_testFunctionDerTwo;
			m_testFunctionDerTwo = NULL;
		}

		if( m_localIntegralPoints != NULL ){
			delete[] m_localIntegralPoints;
			m_localIntegralPoints = NULL;	
		}
		if( m_matrixNeighbor != NULL ){
			delete[] m_matrixNeighbor;
			m_matrixNeighbor = NULL;
		}
		if( m_matrixNeighborInv != NULL ){
			delete[] m_matrixNeighborInv;
			m_matrixNeighborInv = NULL;
		}
		if( m_RBFfunction != NULL ){
			delete[] m_RBFfunction;
			m_RBFfunction = NULL;
		}
		if( m_RBFfunctionDerOne != NULL ){
			delete[] m_RBFfunctionDerOne;
			m_RBFfunctionDerOne=NULL;
		}
		if( m_RBFfunctionDerTwo != NULL){
			delete[] m_RBFfunctionDerTwo;
			m_RBFfunctionDerTwo = NULL;
		}

	}

	//////////////////////////////////////////
	/*
	 * 解线性系统
	 */
	//////////////////////////////////////////
	AssembleStiffMatrix();

	double* b;
	double* x;
	b = new double[ m_numOfPoints ];
	x = new double[ m_numOfPoints ];

	for ( int i = 0; i < 3; i++ ){
		for ( int j = 0; j < m_numOfPoints; j++ ){
			b[ j ] = m_loadVector[ j * 3 + i ];
			x[ j ] = 0;
		}
		int iter; 
		double err;
		linbcg( b, x, m_itol, m_tol, m_itmax, iter, err );
		for ( int j = 0; j < m_numOfPoints; j++ ){
			if( abs( x[ j ] ) > ( radius / 8 )){
				continue;
			}

			m_resultPointSet[ j * 3 + i ] = x[ j * 3 + i ];
		}
	} 

	delete[] b;
	delete[] x;

	///////////////////////////////////
	/*
	 * 释放内存
	 */
	///////////////////////////////////

	for ( int i = 0; i < m_numOfPoints; i++ ){
		m_stiffMatrix[ i ].m_elementStiffMatrix.clear();
		m_massMatrix[ i ].m_elementStiffMatrix.clear();
	}

	if(m_originalNormals!=NULL){
		delete[] m_originalNormals;
		m_originalNormals=NULL;
	}
	if(m_originalMajorDirection!=NULL){
		delete[] m_originalMajorDirection;
		m_originalMajorDirection=NULL;
	}
	if(m_originalMinorDirection!=NULL){
		delete[] m_originalMinorDirection;
		m_originalMinorDirection=NULL;
	}
	m_mapKnearest.clear();

	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_stiffMatrix!=NULL)
		delete[] m_stiffMatrix;
	m_stiffMatrix=NULL;
	if(m_loadVector!=NULL)
		delete[] m_loadVector;
	m_loadVector=NULL;

	return;
}

void gradientFieldFilterRBF::DeleteGradientFieldFilter(){
	/*
	*	函数功能： 释放变量内存
	*/
	if(m_originalPointSet!=NULL)
		delete[] m_originalPointSet;
	m_originalPointSet=NULL;
	if(m_originalColors!=NULL)
		delete[] m_originalColors;
	m_originalColors=NULL;
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	m_originalNormals=NULL;
	if(m_loadVector!=NULL)
		delete[] m_loadVector;
	m_loadVector=NULL;//荷载向量
	if(m_localOrdinate!=NULL)
		delete[] m_loadVector;
	m_localOrdinate=NULL;//             各邻域点的局部坐标
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighbor=NULL;//            插值矩阵
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighborInv=NULL;//        插值矩阵的逆阵
	if(m_localIntegralPoints!=NULL)
		delete[] m_localIntegralPoints;
	m_localIntegralPoints=NULL;  //积分点的坐标
	if(m_testFunction!=NULL)
		delete[] m_testFunction;
	m_testFunction=NULL;//       测试函数在每一个采样点的值
	if(m_RBFfunction!=NULL)
		delete[] m_RBFfunction;
	m_RBFfunction=NULL;//       每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	if(m_shapeFunction!=NULL)
		delete[] m_shapeFunction;
	m_shapeFunction=NULL;//                各邻域点的形函数在各采样点的值
	if(m_determinentMetric!=NULL)
		delete[] m_determinentMetric;
	m_determinentMetric=NULL;//   在采样点处流形标准的代数行列式

	if ( m_sqrtDetMetric != NULL ){
		delete[] m_sqrtDetMetric;
	}
	m_sqrtDetMetric = NULL;

	if( ｍ_scalarProductVector != NULL )
		delete[] ｍ_scalarProductVector;
	ｍ_scalarProductVector = NULL;

	if ( m_gradientProductVector != NULL ){
		delete[] m_gradientProductVector;
	}
	m_gradientProductVector = NULL;

	if(m_metricMatrix!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrix=NULL;//在采样点处流形标准的矩阵
	if(m_metricMatrixInve!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrixInve=NULL;////在采样点处流形标准的矩阵的逆
	if(m_resultPointSet!=NULL)
		delete[] m_resultPointSet;
	m_resultPointSet=NULL;
	//if(m_resultColor!=NULL)
	//	delete[] m_resultColor;
	//m_resultColor=NULL;
	//if(m_meanCurvature!=NULL){
	//	delete[] m_meanCurvature;
	//	m_meanCurvature=NULL;
	//}
}

void gradientFieldFilterRBF::ComputeMapKnearest()
{
	/*
	*	函数功能：计算每一个点的 k邻域
	*/
	double radius=m_radius*m_radius;//存放近邻点到所求点的距离
	std::multimap<double ,int> mapKnearest;//对每一个点建立邻域vector时建立临时map
	int numOfKnearest;//统计一个点的临时邻域大小
	if(m_mapKnearest.begin()!=m_mapKnearest.end()){
		m_mapKnearest.clear();
	}

	//////////////////////////////////////////////////////////////////////////
	//遍历所有点，并且按照距离由小到大存放存放到map中
	for(int i=0;i<m_numOfPoints;i++){
		numOfKnearest=0;
		for(int j=0;j<m_numOfPoints;j++){
			if(j==i)
				continue;
			double distancePointToPoint;
			distancePointToPoint=(m_resultPointSet[i*3]-m_resultPointSet[j*3])*(m_resultPointSet[i*3]-m_resultPointSet[j*3])
				+(m_resultPointSet[i*3+1]-m_resultPointSet[j*3+1])*(m_resultPointSet[i*3+1]-m_resultPointSet[j*3+1])
				+(m_resultPointSet[i*3+2]-m_resultPointSet[j*3+2])*(m_resultPointSet[i*3+2]-m_resultPointSet[j*3+2]);
			if(distancePointToPoint>radius)
				continue;
			mapKnearest.insert(std::multimap<double,int>::value_type(distancePointToPoint,j));
			numOfKnearest+=1;
			if(numOfKnearest>m_kNearest){
				std::multimap<double,int>::iterator mapIterator = mapKnearest.end();
				mapIterator--;
				double w;
				w=(*mapIterator).first;
				mapKnearest.erase(mapIterator);
				numOfKnearest-=1;
			}
		}
		KnearestField fieldKnearest;//临时结构
		fieldKnearest.m_numOfNearest=numOfKnearest;
		fieldKnearest.m_IfBoundary = false;
		std::multimap<double,int>::iterator mapIterator = mapKnearest.begin();
		for(;mapIterator!=mapKnearest.end();mapIterator++){
			fieldKnearest.m_nearest.push_back((*mapIterator).second);
		}
		m_mapKnearest.insert(std::map<int, KnearestField>::value_type(i,fieldKnearest));

		mapKnearest.clear();
	}	
	//////////////////////////////////////////////////////////////////////////
	//// 测试代码
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\knearest.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		std::map<int, KnearestField>::iterator mapKnearestIterator=m_mapKnearest.find(i);
	//		std::vector<int>::iterator vectorIterator=(*mapKnearestIterator).second.m_nearest.begin();
	//		fprintf(fpout,"%d ",i);
	//		for(;vectorIterator!=(*mapKnearestIterator).second.m_nearest.end();vectorIterator++){
	//			fprintf(fpout,"%d ",(*vectorIterator));
	//		}
	//		fprintf(fpout,"\n");
	//	}
	//	fclose(fpout);
	//}
	return;
}


void gradientFieldFilterRBF::ComputeNormal(){
	/*
	*	函数功能：计算每一个点的法向
	*/
	double centroidPosition[3];//重心位置
	double* localVariationVector;//重心到点的位置的向量
	int numOfNearest;

	mat_f8 mat_covaMatrix(3, 3);
	if(m_originalNormals!=NULL)
		delete[] m_originalNormals;
	if(m_originalNormals!=NULL){
		delete[] m_originalNormals;
		m_originalNormals=NULL;
	}
	if(m_originalMajorDirection!=NULL){
		delete[] m_originalMajorDirection;
		m_originalMajorDirection=NULL;
	}
	if(m_originalMinorDirection!=NULL){
		delete[] m_originalMinorDirection;
		m_originalMinorDirection=NULL;
	}
	m_originalNormals=new double[m_numOfPoints*3];
	m_originalMajorDirection=new double[m_numOfPoints*3];
	m_originalMinorDirection=new double[m_numOfPoints*3];
	if(m_meanCurvature!=NULL){
		delete[] m_meanCurvature;
		m_meanCurvature=NULL;
	}
	m_meanCurvature=new double[m_numOfPoints];
	//m_gaussCurvature=new double[m_numOfPoints*2];
	//if(m_meanAreas!=NULL)
	//	delete[] m_meanAreas;
	//m_meanAreas=new double[m_numOfPoints];
	//if(m_meanCurvature!=NULL)
	//	delete[] m_meanCurvature;
	//m_meanCurvature=new double[m_numOfPoints];
	//if(m_resultColor!=NULL)
	//	delete[] m_resultColor;
	//m_resultColor=NULL;
	//m_resultColor=new double[m_numOfPoints*3];


	for(int i=0;i<m_numOfPoints;i++){
		std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
		std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		numOfNearest=(*mapIterator).second.m_numOfNearest;
		for(int j=0;j<3;j++){
			centroidPosition[j]=m_resultPointSet[i*3+j];
		}
		for(;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++){
			centroidPosition[0]+=m_resultPointSet[(*vectorIterator)*3];
			centroidPosition[1]+=m_resultPointSet[(*vectorIterator)*3+1];
			centroidPosition[2]+=m_resultPointSet[(*vectorIterator)*3+2];
		}
		for(int j=0;j<3;j++){
			centroidPosition[j]/=(numOfNearest+1);
		}
		localVariationVector=new double[(numOfNearest+1)*3];
		localVariationVector[0]=m_resultPointSet[i*3]-centroidPosition[0];
		localVariationVector[1]=m_resultPointSet[i*3+1]-centroidPosition[1];
		localVariationVector[2]=m_resultPointSet[i*3+2]-centroidPosition[2];
		vectorIterator=(*mapIterator).second.m_nearest.begin();
		for(int j=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,j++){
			localVariationVector[j*3]=m_resultPointSet[(*vectorIterator)*3]-centroidPosition[0];
			localVariationVector[j*3+1]=m_resultPointSet[(*vectorIterator)*3+1]-centroidPosition[1];
			localVariationVector[j*3+2]=m_resultPointSet[(*vectorIterator)*3+2]-centroidPosition[2];
		}
		//求方差矩阵
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				mat_covaMatrix(j,k)=0;

				for(int m=0;m<(numOfNearest+1);m++){
					mat_covaMatrix(j,k)=mat_covaMatrix(j,k)+localVariationVector[m*3+j]*localVariationVector[m*3+k];
				}
			}
		}
		vec_f8 eval(3);
		mat_f8 evec(3, 3);
		eigen_symm(eval, evec, mat_covaMatrix);
		m_originalNormals[i*3]=evec(0,2);
		m_originalNormals[i*3+1]=evec(1,2);
		m_originalNormals[i*3+2]=evec(2,2);
		m_originalMajorDirection[i*3]=evec(0,0);
		m_originalMajorDirection[i*3+1]=evec(1,0);
		m_originalMajorDirection[i*3+2]=evec(2,0);
		m_originalMinorDirection[i*3]=evec(0,1);
		m_originalMinorDirection[i*3+1]=evec(1,1);
		m_originalMinorDirection[i*3+2]=evec(2,1);
		m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
		//m_gaussCurvature[i*2]=eval(0);
		//m_gaussCurvature[i*2+1]=eval(1);
		//m_meanAreas[i]=abs(eval(0)*eval(2));
		//	m_meanCurvature[i]=eval(2)/(eval(0)+eval(1)+eval(2));
		if(m_meanCurvature[i]*4<=0.005){
		//m_resultColor[i*3]=0;
		//m_resultColor[i*3+1]=0;
		//m_resultColor[i*3+2]=0.5+80*m_meanCurvature[i]*4;
		}
		if((m_meanCurvature[i]*4>0.005)&(m_meanCurvature[i]*4<=0.01)){
	/*	m_resultColor[i*3]=0;
		m_resultColor[i*3+1]=0.45+50*(m_meanCurvature[i]*4-0.005);
		m_resultColor[i*3+2]=0;*/
		}
		if(m_meanCurvature[i]*4>0.01){
		//m_resultColor[i*3]=0.15+5*(m_meanCurvature[i]*4-0.01);
		//m_resultColor[i*3+1]=0;
		//m_resultColor[i*3+2]=0;

		}

		delete[] localVariationVector;
	}
	EstimateNormalDirection();

	//////////////////////////////////////////////////////////////////////////
	// 写法向文件
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\normal.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalNormals[i*3],m_originalNormals[i*3+1],m_originalNormals[i*3+2]);

	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写主方向文件
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\majorDirection.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalMajorDirection[i*3],m_originalMajorDirection[i*3+1],m_originalMajorDirection[i*3+2]);

	//	}
	//	fclose(fpout);
	//}

	////////////////////////////////////////////////////////////////////////////
	//// 写次方向文件
	//filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\minorDirection.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_originalMinorDirection[i*3],m_originalMinorDirection[i*3+1],m_originalMinorDirection[i*3+2]);

	//	}
	//	fclose(fpout);
	//}


	//filename_pw = "D:\\qin\\sourcecode\\marchingcube\\ceshi\\NonShiftmeanCurvature.txt";

	//if((fpout = fopen(filename_pw, "w")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f\n",m_meanCurvature[i]);
	//	}
	//	fclose(fpout);
	//}
}

void gradientFieldFilterRBF::EstimateNormalDirection(){

	/*
	*	函数功能： 估计法向的准确方向
	*/
	//首先确定z值最大的那一点的法向
	int maxZpoint=0;
	for(int i=0;i<m_numOfPoints;i++){
		if(m_resultPointSet[i*3+2]>m_resultPointSet[maxZpoint*3+2])
			maxZpoint=i;

	}
	if(m_originalNormals[maxZpoint*3+2]<0){
		for(int i=0;i<3;i++){
			m_originalNormals[maxZpoint*3+i]=-m_originalNormals[maxZpoint*3+i];
		}
	}

	int* parent=new int[m_numOfPoints];//存放父节点

	std::multimap<double,int> valueKeySeries;//根据点的代价索引的map
	std::map<int,double> pointKeySeries;//根据点索引的map

	//初始化map
	for(int i=0;i<maxZpoint;i++){
		pointKeySeries.insert(std::map<int,double>::value_type(i,i+10));
		valueKeySeries.insert(std::map<double,int>::value_type(i+10,i));
	}
	for(int i=maxZpoint+1;i<m_numOfPoints;i++){
		pointKeySeries.insert(std::map<int,double>::value_type(i,i+10));
		valueKeySeries.insert(std::map<double,int>::value_type(i+10,i));
	}
	pointKeySeries.insert(std::map<int,double>::value_type(maxZpoint,0));
	valueKeySeries.insert(std::map<double,int>::value_type(0,maxZpoint));

	std::map<int,double>::iterator mapIterator = pointKeySeries.begin();
	std::multimap<double,int>::iterator multimapIterator=valueKeySeries.begin();

	while(multimapIterator!=valueKeySeries.end()){
		int pointExtract;
		if((*multimapIterator).first>1)
			break;
		pointExtract=(*multimapIterator).second;
		valueKeySeries.erase(multimapIterator);
		mapIterator=pointKeySeries.find(pointExtract);
		pointKeySeries.erase(mapIterator);
		std::map<int, KnearestField>::iterator mapNearestIterator=m_mapKnearest.find(pointExtract);
		std::vector<int>::iterator vectorIterator=(*mapNearestIterator).second.m_nearest.begin();


		for(;vectorIterator!=(*mapNearestIterator).second.m_nearest.end();vectorIterator++){
			int neighborPoint=(*vectorIterator);
			mapIterator=pointKeySeries.find(neighborPoint);
			if(mapIterator!=pointKeySeries.end()){
				multimapIterator=valueKeySeries.find((*mapIterator).second);
				if(neighborPoint!=(*multimapIterator).second){
					continue;
				}
				double relativeCost,normalDotProduct;
				double pointExtractNormal[3],neighborPointNormal[3];
				for(int j=0;j<3;j++){
					pointExtractNormal[j]=m_originalNormals[pointExtract*3+j];
					neighborPointNormal[j]=m_originalNormals[neighborPoint*3+j];
				}
				normalDotProduct=Vector3Vector(pointExtractNormal,neighborPointNormal);
				if(normalDotProduct<0){
					normalDotProduct=-normalDotProduct;
					for(int j=0;j<3;j++){
						m_originalNormals[neighborPoint*3+j]=-m_originalNormals[neighborPoint*3+j];
					}
				}
				relativeCost=1-normalDotProduct;
				if(relativeCost<(*mapIterator).second){
					parent[neighborPoint]=pointExtract;
					(*mapIterator).second=relativeCost;
					valueKeySeries.erase(multimapIterator);
					valueKeySeries.insert(std::map<double,int>::value_type(relativeCost,neighborPoint));
				}
			}
		}
		mapIterator=pointKeySeries.begin();
		multimapIterator=valueKeySeries.begin();
	}

	/*if(mapIterator!=pointKeySeries.end()){
	int numOfPointsLeave;
	double* pointSetLeave;
	numOfPointsLeave=pointKeySeries.size();
	pointSetLeave=new double[numOfPointsLeave*3];
	int i=0;
	while(mapIterator!=pointKeySeries.end()){
	int tempPointIndex;
	tempPointIndex=(*mapIterator).first;
	for(int j=0;j<3;j++){
	pointSetLeave[i*3+j]=m_originalPointSet[tempPointIndex*3+j];
	}			
	mapIterator++;
	i++;
	}

	pointKeySeries.clear();
	valueKeySeries.clear();
	delete[] parent;
	parent=NULL;

	EstimateNormalDirection(numOfPointsLeave,pointSetLeave,numOfEachPointNeighbor,allPointsNeighbor,allPointsNormal);
	delete[] pointSetLeave;
	}*/
	pointKeySeries.clear();
	valueKeySeries.clear();
	if(parent!=NULL)
		delete[] parent;
	return;

}

void gradientFieldFilterRBF::ComputeNearestParameter(int indexPointSet){
	/*
	*	函数功能： 计算任一点处的积分半径，计算其邻域各点的局部坐标，
	*  变量说明：
	*  int indexPointSet        索引点	
	*  double m_localRBFradius;            径向基函数半径
	*  double* m_localOrdinate             各邻域点的局部坐标
	*  double* m_localIntegralPoint        积分点的坐标
	*  double  m_localIntegralRadius       积分域的半径
	*/
	double localPoint[3], localNormal[3],localMajorDirection[3],localMinorDirection[3];

	m_localRBFradius=0;
	for(int i=0;i<3;i++){
		localMajorDirection[i]=m_originalMajorDirection[indexPointSet*3+i];
		localMinorDirection[i]=m_originalMinorDirection[indexPointSet*3+i];
		localNormal[i]=m_originalNormals[indexPointSet*3+i];
	}
	std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(indexPointSet);
	std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
	m_numOfNeighbors=(*mapIterator).second.m_nearest.size()+1;
	if(m_localOrdinate!=NULL){
		delete[] m_localOrdinate;
		m_localOrdinate=NULL;
	}
	m_localOrdinate=new double[m_numOfNeighbors*2];
	m_localOrdinate[0]=0;
	m_localOrdinate[1]=0;

	for(int i=1;vectorIterator!=(*mapIterator).second.m_nearest.end();vectorIterator++,i++){
		double vectorNeighborToPoint[3];
		double localRadius;
		double tempRadiusSquare;
		for(int j=0;j<3;j++){
			vectorNeighborToPoint[j]=m_resultPointSet[(*vectorIterator)*3+j]-m_resultPointSet[indexPointSet*3+j];
		}
		m_localOrdinate[i*2]=Vector3Vector(vectorNeighborToPoint,localMajorDirection);
		m_localOrdinate[i*2+1]=Vector3Vector(vectorNeighborToPoint,localMinorDirection); 
		//	m_localOrdinate[i*3+2]=Vector3Vector(vectorNeighborToPoint,localNormal);
		tempRadiusSquare=m_localOrdinate[i*2]*m_localOrdinate[i*2]+m_localOrdinate[i*2+1]*m_localOrdinate[i*2+1];
		if(tempRadiusSquare>m_localRBFradius)
			m_localRBFradius=tempRadiusSquare;
	}
	m_localRBFradius=sqrt(m_localRBFradius);

	//////////////////////////////////////////////////////////////////////////
	/*
	*	写局部坐标文件
	*/
	//CString filename_pw = "D:\\qin\\MeshlessFilter\\CeshiData\\localOrdinalte.txt";
	//FILE *fpout;
	//if((fpout = fopen(filename_pw, "a")) == NULL)
	//{
	//	int dkjkd;
	//	//MessageBox("can't open the file!");
	//}
	//else
	//{
	//	fprintf(fpout,"%d\n",indexPointSet);
	//	for(int i=0;i<m_numOfNeighbors;i++){
	//		fprintf(fpout,"%f %f\n",m_localOrdinate[i*3],m_localOrdinate[i*3+1]);

	//	}
	//	fclose(fpout);
	//}


	return;

}

void gradientFieldFilterRBF::ComputeTestFunction(int indexPointSet){

	/*
	*	计算测试函数在积分点的函数值及其导数值
	*/
	if(m_localIntegralPoints!=NULL){
		delete[] m_localIntegralPoints;
		m_localIntegralPoints=NULL;
	}
	m_localIntegralPoints=new double[m_numOfIntegralPoint*2];
	if(m_testFunction!=NULL){
		delete[] m_testFunction;
		m_testFunction=NULL;
	}
	m_testFunction=new double[m_numOfIntegralPoint];
	if(m_testFunctionDerOne!=NULL){
		delete[] m_testFunctionDerOne;
		m_testFunctionDerOne=NULL;
	}
	m_testFunctionDerOne=new double[m_numOfIntegralPoint];
	if(m_testFunctionDerTwo!=NULL){
		delete[] m_testFunctionDerTwo;
		m_testFunctionDerTwo=NULL;
	}
	m_testFunctionDerTwo=new double[m_numOfIntegralPoint];
	m_localIntegralRadius=m_localRBFradius/2;

	//////////////////////////////////////////////////////////////////////////
	/*
	*	计算积分点的局部坐标
	*/
	//////////////////////////////////////////////////////////////////////////
	for (int i=0;i<m_numOfIntegralPoint;i++){
		m_localIntegralPoints[i*2]=m_ordinateAndWeightCoff[i*3]*m_localIntegralRadius;
		m_localIntegralPoints[i*2+1]=m_ordinateAndWeightCoff[i*3+1]*m_localIntegralRadius;		
	}
	//////////////////////////////////////////////////////////////////////////
	/*
	*	计算形函数及其导数  选择 c和r相等
	*/
	//////////////////////////////////////////////////////////////////////////
	double tempEXP,tempOneMinusEXP;// 计算测试函数时候的一些常数
	double temp=-1;
	tempEXP=exp(temp);
	tempOneMinusEXP=1-tempEXP;
	double integralRaiusSquare, radiusCSquare;

	integralRaiusSquare=m_localIntegralRadius*m_localIntegralRadius;
	radiusCSquare=integralRaiusSquare;

	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempDistanceSquare;//积分点与原点之间的距离的平方
		tempDistanceSquare=m_localIntegralPoints[i*2]*m_localIntegralPoints[i*2]
		+m_localIntegralPoints[i*2+1]*m_localIntegralPoints[i*2+1];
		if(tempDistanceSquare>integralRaiusSquare){ 
			m_testFunction[i]=0;
			m_testFunctionDerOne[i]=0;
			m_testFunctionDerTwo[i]=0;

		}
		else{
			double tempEXPDistance;
			tempEXPDistance=exp(-tempDistanceSquare/radiusCSquare);
			m_testFunction[i]=(tempEXPDistance-tempEXP)/tempOneMinusEXP;
			m_testFunctionDerOne[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2]/radiusCSquare/tempOneMinusEXP;
			m_testFunctionDerTwo[i]=-2*tempEXPDistance*m_localIntegralPoints[i*2+1]/radiusCSquare/tempOneMinusEXP;
		}

	}
	return;

}

void gradientFieldFilterRBF::ComputeMatrixNeighbor(int indexPointSet){

	/*
	*	函数功能：计算插值矩阵及其逆阵
	*  变量说明：
	*  double* m_matrixNeighbor; 插值矩阵
	*  double* m_matrixNeighborInv; 插值矩阵的逆阵  
	*/

	if ( m_matrixNeighbor != NULL ){
		delete[] m_matrixNeighbor;
		m_matrixNeighbor = NULL;
	}

	m_matrixNeighbor = new double [ m_numOfNeighbors * m_numOfNeighbors];
	if ( m_matrixNeighborInv != NULL ){
		delete[] m_matrixNeighborInv;
		m_matrixNeighborInv = NULL;
	}

	m_matrixNeighborInv = new double [ m_numOfNeighbors * m_numOfNeighbors ];

		for ( int i = 0; i < m_numOfNeighbors; i++ ){
			// 计算每一行 在不同点不同基函数的值

			for ( int j = 0; j < m_numOfNeighbors ; j++ ){
				double tempDistance; // d/r
				double tempVector[2];
				for( int k = 0; k < 2; k++ ){
					tempVector[ k ] = m_localOrdinate[ i * 2 + k ] - m_localOrdinate[ j * 2 + k ];
				}
				tempDistance = tempVector[ 0 ] * tempVector[ 0 ] + tempVector[ 1 ] * tempVector[ 1 ];
				tempDistance = sqrt( tempDistance );
				if( tempDistance > m_localRBFradius ){
					m_matrixNeighbor[ i * m_numOfNeighbors + j ] = 0;
				}
				else{
					double oneMinusTempDistance;
					double tempOne, tempTwo;
					tempDistance /= m_localRBFradius;
					oneMinusTempDistance = 1 - tempDistance;
					tempOne = oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance;
					tempTwo = m_coffRBFunction[ 0 ] + m_coffRBFunction[ 1 ] * tempDistance
						+ m_coffRBFunction[ 2 ] * tempDistance * tempDistance
						+ m_coffRBFunction[ 3 ] * tempDistance * tempDistance * tempDistance
						+ m_coffRBFunction[ 4 ] * tempDistance * tempDistance * tempDistance * tempDistance;
					m_matrixNeighbor[ i * m_numOfNeighbors + j ] = tempOne * tempTwo;
				}
			}
		}
		mat_f8 mat_localMatrixNeighbor( m_numOfNeighbors, m_numOfNeighbors );
		mat_f8 mat_localMatrixNeighborInv( m_numOfNeighbors, m_numOfNeighbors );

		for( int i = 0; i < m_numOfNeighbors; i++ ){
			for ( int j = 0; j < m_numOfNeighbors; j++ ){
				mat_localMatrixNeighbor( i, j ) = m_matrixNeighbor[i * m_numOfNeighbors + j ];
			}
		}

		if( !matrix_inverse( mat_localMatrixNeighbor)){
			///报求解逆矩阵出错
		}

		for( int i = 0; i < m_numOfNeighbors; i++ ){
			for( int j = 0; j < m_numOfNeighbors; j++ ){
				m_matrixNeighborInv[ i * m_numOfNeighbors + j ] = mat_localMatrixNeighbor( i, j );
			}
		}

		return;

}

void gradientFieldFilterRBF::ComputeMatrixRelatedSampling(int indexPointSet){

	/*
	 * 函数功能： 计算与采样点有关的矩阵
	 * 变量说明：
	 * int indexPointSet 索引点
	 * double* m_RBFfunction; 每一行为在每一个邻域点的径向基函数在每一个积分点的函数值
	 */
	/////////////////////////////////////////////////////////////////////////
	/*
	 * 分配内存空间
	 */
	if ( m_RBFfunction != NULL ){
		delete[] m_RBFfunction;
		m_RBFfunction = NULL;
	}
	m_RBFfunction = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];

	if ( m_RBFfunctionDerOne != NULL ){
		delete[] m_RBFfunctionDerOne;
		m_RBFfunctionDerOne = NULL;
	}
	m_RBFfunctionDerOne = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];

	if( m_RBFfunctionDerTwo != NULL ){
		delete[] m_RBFfunctionDerTwo;
		m_RBFfunctionDerTwo = NULL;
	}
	m_RBFfunctionDerTwo = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];

	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){
			double tempDistanceD; /// d/r
			double tempVector[ 2 ];
			for ( int k = 0; k < 2; k++ ){
				tempVector[ k ] = m_localIntegralPoints[ j * 2 + k] - m_localOrdinate[ i * 2 + k ];
			}
			tempDistanceD = tempVector[ 0 ] * tempVector[ 0 ] + tempVector[ 1 ] * tempVector[ 1 ];
			tempDistanceD = sqrt( tempDistanceD );

			if( tempDistanceD > m_localRBFradius ){
				m_RBFfunction[ i * m_numOfIntegralPoint + j ] = 0;
				m_RBFfunctionDerTwo[ i * m_numOfIntegralPoint + j ] = 0;
				m_RBFfunctionDerOne[ i * m_numOfIntegralPoint + j ] = 0;
			}
			else{
				double oneMinusTempDistance, tempDistance;
				double tempOne, tempTwo;
				for( int k = 0; k < 2; k++ ){
					tempVector[ k ] /= tempDistanceD;
				}
				tempDistance = tempDistanceD / m_localRBFradius;
				oneMinusTempDistance = 1 - tempDistance;
				tempOne = oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance*oneMinusTempDistance;
				tempTwo = m_coffRBFunction[ 0 ] + m_coffRBFunction[ 1 ] * tempDistance
					+ m_coffRBFunction[ 2 ] * tempDistance * tempDistance
					+ m_coffRBFunction[ 3 ] * tempDistance * tempDistance * tempDistance
					+ m_coffRBFunction[ 4 ] * tempDistance * tempDistance * tempDistance * tempDistance;\
				m_RBFfunction[ i * m_numOfIntegralPoint + j ] = tempOne * tempTwo;
				m_RBFfunctionDerOne[ i * m_numOfIntegralPoint + j ] = -5 * ( tempOne / oneMinusTempDistance ) * tempTwo * tempVector[ 0 ] / m_localRBFradius;
				m_RBFfunctionDerTwo[ i * m_numOfIntegralPoint + j ] = -5 * ( tempOne / oneMinusTempDistance ) * tempTwo * tempVector[ 1 ] / m_localRBFradius;
				tempTwo=40/m_localRBFradius
					+96*tempDistanceD/m_localRBFradius/m_localRBFradius
					+75*tempDistanceD*tempDistanceD/m_localRBFradius/m_localRBFradius/m_localRBFradius
					+20*tempDistanceD*tempDistanceD*tempDistanceD/m_localRBFradius/m_localRBFradius/m_localRBFradius/m_localRBFradius;
				m_RBFfunctionDerOne[i*m_numOfIntegralPoint+j]+=tempOne*tempTwo*tempVector[0];
				m_RBFfunctionDerTwo[i*m_numOfIntegralPoint+j]+=tempOne*tempTwo*tempVector[1];

			}
		}
	}
	return;
}
void gradientFieldFilterRBF::ComputeManifoldMetric(int indexPointSet){

	/*
	 * double* m_determinentMetric;//   在采样点处流形标准的代数行列式
	 * double* m_metricMatrix;//在采样点处流形标准的矩阵
	 * double* m_metricMatrixInve;////在采样点处流形标准的矩阵的逆
	 */
	////////////////////////////////////
	/*
	 * 分配内存空间
	 */
	if ( m_determinentMetric != NULL ) {
		delete[] m_determinentMetric;
		m_determinentMetric = NULL;
	}
	m_determinentMetric = new double[ m_numOfIntegralPoint ];

	if ( m_sqrtDetMetric != NULL ){
		delete[] m_sqrtDetMetric;
		m_sqrtDetMetric = NULL;
	}
	m_sqrtDetMetric = new double[ m_numOfIntegralPoint ];

	if ( m_metricMatrix != NULL ){
		delete[] m_metricMatrix;
		m_metricMatrix = NULL;
	}
	m_metricMatrix = new double[ m_numOfIntegralPoint * 4 ];

	if ( m_metricMatrixInve != NULL ){
		delete[] m_metricMatrixInve;
		m_metricMatrixInve = NULL;
	}
	m_metricMatrixInve = new double[ m_numOfIntegralPoint * 4 ];

	std::map< int, KnearestField >::iterator mapIterator = m_mapKnearest.find( indexPointSet );
	std::vector< int >::iterator vectorIterator = ( *mapIterator ).second.m_nearest.begin();
	m_numOfNeighbors = ( *mapIterator ).second.m_nearest.size()+1;

	double x_Der_One[ 3 ]; //// 三维空间的位置坐标对局部空间第一个坐标的导数
	double x_Der_Two[ 3 ]; //// 三维空间的位置坐标对局部空间第二个坐标的导数

	for ( int i = 0; i < m_numOfIntegralPoint; i++ ){
		/// 临时变量
		for ( int j = 0; j < 3; j++ ){
			x_Der_One[ j ] = 0;
			x_Der_Two[ j ] = 0;
		}
		
		for ( int j = 0; j < 3; j++ ){

			for ( int k = 0; vectorIterator!=(*mapIterator).second.m_nearest.end(); k++, vectorIterator++ ){
				double tempVarOne = 0;
				double tempVarTwo = 0;
				for ( int m = 0; m < m_numOfNeighbors; m++ ){
					tempVarOne += m_RBFfunctionDerOne[ m * m_numOfIntegralPoint + i ] * m_matrixNeighborInv[ m * m_numOfIntegralPoint + k ];
					tempVarTwo += m_RBFfunctionDerTwo[ m * m_numOfIntegralPoint + i ] * m_matrixNeighborInv[ m * m_numOfIntegralPoint + k ];
				}
				x_Der_One[ j ] += tempVarOne *  m_originalPointSet[ ( *vectorIterator ) * 3 + j ];
				x_Der_Two[ j ] += tempVarTwo *  m_originalPointSet[ ( *vectorIterator ) * 3 + j ];			
			}
		}

		m_metricMatrix[ i * m_numOfIntegralPoint ] = x_Der_One[ 0 ] * x_Der_One[ 0 ]
												   + x_Der_One[ 1 ] * x_Der_One[ 1 ]
												   + x_Der_One[ 2 ] * x_Der_One[ 2 ];
	   	m_metricMatrix[ i * m_numOfIntegralPoint + 1 ] = x_Der_One[ 0 ] * x_Der_Two[ 0 ]
												       + x_Der_One[ 1 ] * x_Der_Two[ 1 ]
												       + x_Der_One[ 2 ] * x_Der_Two[ 2 ];

		m_metricMatrix[ i * m_numOfIntegralPoint + 2 ] = m_metricMatrix[ i * m_numOfIntegralPoint + 1 ];
	    m_metricMatrix[ i * m_numOfIntegralPoint + 1 ] = x_Der_Two[ 0 ] * x_Der_Two[ 0 ]
												       + x_Der_Two[ 1 ] * x_Der_Two[ 1 ]
												       + x_Der_Two[ 2 
													   
													   ] * x_Der_Two[ 2 ];

	   mat_f8 mat_metricMatrix( 2, 2 );
	   mat_f8 mat_metrixMatrixInv( 2, 2 );

	   for( int j = 0; i < 2; j++ ){
			for ( int k = 0; k < 2; k++ ){
				mat_metricMatrix( j, k ) = m_metricMatrix[i * m_numOfIntegralPoint + j * 2 + k ];
			}
		}

		if( !matrix_inverse( mat_metricMatrix)){
			///报求解逆矩阵出错
		}
        
		for( int j = 0; j < 2; j++ ){
			for( int k = 0; k < 2; k++ ){
				m_metricMatrixInve[ i * m_numOfIntegralPoint + j * 2 + k ] = mat_metricMatrix( j ,k );
			}
		}

		m_determinentMetric[ i ] = m_metricMatrix[ i * m_numOfIntegralPoint ] * m_metricMatrix[ i * m_numOfIntegralPoint + 3 ] 
						         - m_metricMatrix[ i * m_numOfIntegralPoint + 1 ] * m_metricMatrix[ i * m_numOfIntegralPoint + 2 ];
	    m_sqrtDetMetric[ i ] = sqrt( m_determinentMetric[ i ] );

	}
	return;
}

void gradientFieldFilterRBF::ComputeScalarProduct(int indexPointSet){

	/*
	 *  函数功能： 计算RBF函数与测试函数的内积
	 *  变量说明：
	 *  int indexPointSet        索引点
	 *  double*　ｍ_scalarProcuctVector  在采样点所在切空间中，n个RBF函数与测试函数积分的向量
	 */
	
	if ( ｍ_scalarProductVector != NULL ){
		delete[] ｍ_scalarProductVector;
		ｍ_scalarProductVector = NULL;
	}
	std::map< int, KnearestField >::iterator mapIterator = m_mapKnearest.find( indexPointSet );
	std::vector< int >::iterator vectorIterator = ( *mapIterator ).second.m_nearest.begin();
	m_numOfNeighbors = ( *mapIterator ).second.m_nearest.size()+1;

	ｍ_scalarProductVector = new double[ m_numOfNeighbors ]; 

	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		ｍ_scalarProductVector[ i ] = 0;
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){

			ｍ_scalarProductVector[ i ] += m_RBFfunction[ i * m_numOfIntegralPoint + j ] * m_testFunction[ j ]
									* m_ordinateAndWeightCoff[ j * 3 + 2 ] * m_sqrtDetMetric[ j ];
		}
	}

}

void gradientFieldFilterRBF::ComputeVectorProductVector(int indexPointSet){
	/////////////////////////////////////////////////////
	// 首先由导数计算梯度
	/////////////////////////////////////////////////////
    /*
	 * 临时变量
	 */
	double* gradient_RBF_One;
	double* gradient_RBF_Two;
	double* gradient_test_One;
	double* gradient_test_Two;
	gradient_RBF_One = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];
	gradient_RBF_Two = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];
	gradient_test_One = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];
	gradient_test_Two = new double[ m_numOfNeighbors * m_numOfIntegralPoint ];

	for ( int i = 0; i < m_numOfNeighbors * m_numOfIntegralPoint; i++ ){
		gradient_RBF_One[ i ] = 0;
		gradient_RBF_Two[ i ] = 0;
		gradient_test_One[ i ] = 0;
		gradient_test_Two[ i ] = 0;
	}

	for ( int i = 0; i < m_numOfNeighbors; i++ ){ // 有m_numOfNeighbors个RBF函数，每个函数和test function有一个积分
		                                          // 因此 每一个RBF函数在每一个积分点处要求一个梯度出来
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){ // 有m_numOfIntegralPoint个积分点
			
			gradient_RBF_One[ i * m_numOfIntegralPoint + j ] += m_metricMatrixInve[ j * 4 ] * m_RBFfunctionDerOne[ i * m_numOfIntegralPoint + j ]
															  + m_metricMatrixInve[ j * 4 + 1 ] * m_RBFfunctionDerTwo[ i * m_numOfIntegralPoint + j ];
		    gradient_RBF_Two[ i * m_numOfIntegralPoint + j ] += m_metricMatrixInve[ j * 4 + 1 ] * m_RBFfunctionDerOne[ i * m_numOfIntegralPoint + j ]
															  + m_metricMatrixInve[ j * 4 + 3 ] * m_RBFfunctionDerTwo[ i * m_numOfIntegralPoint + j ];
		}
	}

	for ( int i = 0; i < m_numOfIntegralPoint; i++ ){
		gradient_test_One[ i ] += m_metricMatrixInve[ i * 4 ] * m_testFunctionDerOne[ i ]
							    + m_metricMatrixInve[ i * 4 + 1 ] * m_testFunctionDerTwo[ i ];
		gradient_test_Two[ i ] += m_metricMatrixInve[ i * 4 + 1 ] * m_testFunctionDerOne[ i ]
								+ m_metricMatrixInve[ i * 4 + 3 ] * m_testFunctionDerTwo[ i ];

	}

	////////////////////////////////////////////////////////////////////
	/*
	 *计算每一个RBF函数与test function在流形上的积分
	 */
	////////////////////////////////////////////////////////////////////
	
	if ( m_gradientProductVector != NULL ){
		delete[] m_gradientProductVector;
		m_gradientProductVector = NULL;
	}
	m_gradientProductVector = new double[ m_numOfNeighbors ]; 
	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		m_gradientProductVector[ i ] = 0;
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){
			double temp_Integral;// 没有乘以每一个积分点的积分系数
			temp_Integral = ( gradient_RBF_One[ i * m_numOfIntegralPoint + j ] * m_metricMatrix[ j * 4 ]
							  + gradient_RBF_Two[ i * m_numOfIntegralPoint + j ] * m_metricMatrix[ j * 4 + 2 ] )
						    * gradient_test_One[ j ]
						   + ( gradient_RBF_One[ i * m_numOfIntegralPoint + j ] * m_metricMatrix[ j * 4 + 1 ]
						      + gradient_RBF_Two[ i * m_numOfIntegralPoint + j ] * m_metricMatrix[ j * 4 + 3 ] )
						    * gradient_test_Two[ j ];
		  m_gradientProductVector[ i ] += temp_Integral * m_ordinateAndWeightCoff[ j * 3 + 2 ];
		}
	}
	
	delete[] gradient_RBF_One;
	delete[] gradient_RBF_Two;
	delete[] gradient_test_One;
	delete[] gradient_test_Two;

	return;
}

void gradientFieldFilterRBF::ComputeStiffMatrix(int indexPointSet){

	/*
	 *	函数功能： 计算刚度矩阵和荷载向量
	 *  在每一个索引点indexPointSet，为刚度矩阵计算一行，为荷载向量计算一个值
	 *  变量说明： 
	 *  int indexPointSet                    索引点
	 *  double* m_loadVector                  荷载向量
	 *  StiffMatrix* m_stiffMatrix           刚度矩阵中由梯度积产生的值
	 *  StiffMatrix* m_massMatrix            刚度矩阵中由标量积产生的值
	 */
	///////////////////////////////
	/*
	 *  计算 m_massMatrix
	 */
	///////////////////////////////
	    
	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		double temp_one = 0;

		for ( int j = 0; j < m_numOfNeighbors; j++ ){
			temp_one += ｍ_scalarProductVector[ j ] * m_matrixNeighborInv[ j * m_numOfNeighbors + i ];
		}
		m_massMatrix[indexPointSet].m_elementStiffMatrix.push_back( temp_one );
	}
	///////////////////////////////
	/*
	 *  计算 m_stiffMatrix
	 */
	///////////////////////////////
	
	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		double temp_two = 0;

		for ( int j = 0; j < m_numOfNeighbors; j++ ){
			temp_two += m_gradientProductVector[ j ] * m_matrixNeighborInv[ j * m_numOfNeighbors + i ];
		}
		m_stiffMatrix[indexPointSet].m_elementStiffMatrix.push_back( temp_two );			
	}

	return;
}

void gradientFieldFilterRBF::AssembleStiffMatrix()
{
	//////////////////////////////////////////////////////////////////////////	
	/*
	 *	函数功能： 组装刚度矩阵
	 *  计算方程组右边和刚度矩阵，先计算右边，再改刚度矩阵 组合质量矩阵和刚度矩阵
	 */
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/*
	 * 计算方程组右边 计算向量
	 */
	//////////////////////////////////////////////////////////////////////////
	double beta;
	if ( m_loadVector != NULL ){
		delete[] m_loadVector;
		m_loadVector = NULL;
	}
	m_loadVector = new double[ m_numOfPoints * 3 ];


	for ( int i = 0; i < m_numOfPoints; i++ ){
			std::vector< double >::iterator iteratorStiffMatrix = m_stiffMatrix[i].m_elementStiffMatrix.begin();
			std::vector< double >::iterator iteratorMassMatrix = m_massMatrix[i].m_elementStiffMatrix.begin();
			std::map< int, KnearestField >::iterator mapIterator=m_mapKnearest.find(i);
			std::vector< int >::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
		for ( int j = 0; j < 3; j++ ){
			m_loadVector[ i * 3 + j ] = 0;
			m_loadVector[ i * 3 + j ] = ( ( *iteratorStiffMatrix ) * beta + ( *iteratorMassMatrix ) ) * m_originalPointSet[ i * 3 + j ];
		}
		iteratorStiffMatrix++;
		iteratorMassMatrix++;
		for ( iteratorStiffMatrix; 
			iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMatrix.end()
			|| vectorIterator != (*mapIterator).second.m_nearest.end(); 
		    iteratorStiffMatrix++, iteratorMassMatrix++, vectorIterator++ ){
			for ( int j = 0; j < 3; j++ ){
				m_loadVector[ i * 3 + j ] +=
					( ( *iteratorStiffMatrix ) * beta + ( *iteratorMassMatrix ) ) * m_originalPointSet[ ( *vectorIterator) * 3 + j ];				
			}
		}	

	}
	
	//////////////////////////////////////////////////////////////////////////
	/*
	 * 计算方程组左边 计算刚度矩阵
	 */
	//////////////////////////////////////////////////////////////////////////
	double alpha;
	for ( int i = 0; i < m_numOfPoints; i++ ){
		std::vector< double >::iterator iteratorStiffMatrix = m_stiffMatrix[i].m_elementStiffMatrix.begin();
		std::vector< double >::iterator iteratorMassMatrix = m_massMatrix[i].m_elementStiffMatrix.begin();
		for ( iteratorStiffMatrix;
			iteratorStiffMatrix != m_stiffMatrix[ i ].m_elementStiffMatrix.end();
			iteratorStiffMatrix++, iteratorMassMatrix++ ){
				( *iteratorStiffMatrix ) = ( *iteratorStiffMatrix ) + ( *iteratorMassMatrix ) * alpha;
		}
	}				
}


//////////////////////////////////////////////////////////////////////////
	/*
	*	pbcg 求稀疏线性方程组
	*/
	//////////////////////////////////////////////////////////////////////////

	void gradientFieldFilterRBF::linbcg(double* b,double* x,int itol,double tol,int itmax, int & iter, double &err)
	{
	double ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zminrm,znrm;
	const double EPS=1.0E-14;
	int j;
	int n=m_numOfPoints;
	double* p=new double[n];
	double* pp=new double[n];
	double* r=new double[n];
	double* rr=new double[n];
	double* z=new double[n];
	double* zz=new double[n];
	iter=0;
	atimes(x,r,0);
	for(j=0;j<n;j++){
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//atimes(r,rr,0);
	if(itol==1){
		bnrm=snrm(b,itol);
		asolve(r,z,0);
	}
	else if (itol==2){
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
	}
	else if(itol==3 || itol==4){
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
		znrm=snrm(z,itol);
	}

	while(iter <itmax){
		++iter;
		asolve(rr,zz,1);
		for(bknum=0.0,j=0;j<n;j++) bknum+=z[j]*rr[j];
		if(iter==1){
			for(j=0;j<n;j++){
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else{
			bk=bknum/bkden;
			for(j=0;j<n;j++){
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(p,z,0);
		for(akden=0.0,j=0;j<n;j++) akden+=z[j]*pp[j];
		ak=bknum/akden;
		atimes(pp,zz,1);
		for(j=0;j<n;j++){
			x[j]+=ak*p[j];
			r[j]-=ak*z[j];
			rr[j]-=ak*zz[j];
		}
		asolve(r,z,0);
		if(itol==1)
			err=snrm(r,itol)/bnrm;
		else if(itol==2)
			err=snrm(z,itol)/bnrm;
		else if (itol==3||itol||4){
			zminrm=znrm;
			znrm=snrm(z,itol);
			if(fabs(zminrm-znrm)>EPS*znrm){
				dxnrm=fabs(ak)*snrm(p,itol);
				err=znrm/fabs(zminrm-znrm)*dxnrm;
			}else{
				err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(x,itol);
			if(err<=0.5*xnrm) err/=xnrm;
			else{
				err=znrm/bnrm;
				continue;
			}

		}
		if(err<=tol) break;

	}

	}
	double gradientFieldFilterRBF::snrm(double* sx,const int itol)
	{
	int i,isamax;
	double ans;
	int n=m_numOfPoints;
	if(itol<=3){
		ans=0.0;
		for(i=0;i<n;i++) ans+=sx[i]*sx[i];
		return sqrt(ans);
	}else{
		isamax=0;
		for(i=0;i<n;i++){
			if(fabs(sx[i])>fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
	}
	void gradientFieldFilterRBF::atimes(double* x,double* r,const int itrnsp)
	{
	if(itrnsp==0){
		for (int i=0;i<m_numOfPoints;i++){
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMatrix.begin();
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			r[i]=(*iteratorStiffMatrix)*x[i];

			for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMatrix.end();iteratorStiffMatrix++,vectorIterator++){
				r[i]+=(*iteratorStiffMatrix)*x[(*vectorIterator)];				
			}		
		}
	}
	else{
		for(int i=0;i<m_numOfPoints;i++){
			r[i]=0;
		}
		for (int i=0;i<m_numOfPoints;i++){
			std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMatrix.begin();
			std::map<int, KnearestField>::iterator mapIterator=m_mapKnearest.find(i);
			std::vector<int>::iterator vectorIterator=(*mapIterator).second.m_nearest.begin();
			r[i]+=(*iteratorStiffMatrix)*x[i];

			for(iteratorStiffMatrix++;iteratorStiffMatrix!=m_stiffMatrix[i].m_elementStiffMatrix.end();iteratorStiffMatrix++,vectorIterator++){
				r[(*vectorIterator)]+=(*iteratorStiffMatrix)*x[i];				
			}		
		}

	}
	return;

	}
	void gradientFieldFilterRBF::asolve(double* b,double* x,const int itrnsp)
	{
	int i;
	int n=m_numOfPoints;
	for(i=0;i<n;i++){
		std::vector<double>::iterator iteratorStiffMatrix=m_stiffMatrix[i].m_elementStiffMatrix.begin();
		x[i]=((*iteratorStiffMatrix)!=0.0 ? b[i]/(*iteratorStiffMatrix):b[i]);
	}
	return;
	}



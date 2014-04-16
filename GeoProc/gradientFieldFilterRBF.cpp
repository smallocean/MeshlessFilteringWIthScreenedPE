#include "gradientFieldFilterRBF.h"
#include <math.h>
#include "Matrix.h"
#include ".\mathlib\mathlib.h"

using namespace MATHLIB;
const int gradientFieldFilterRBF::m_numOfIntegralPoint=40;
const double gradientFieldFilterRBF::m_ordinateAndWeightCoff[120]={ 
//////////////////////////////////////////////////////////////////////////
/*  
*   ��һ��Ϊ�뾶Ϊ1ʱ���x����
*   �ڶ���Ϊ�뾶Ϊ1ʱ���y����
*   ������Ϊ����һά��˹���ֵļ�Ȩֵ�ĳ˻�
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

	/// ��ʼ����������
	m_originalPointSet = NULL;
	m_numOfPoints = NULL;
	m_numOfPoints = 0;
	m_kNearest = 0;
	m_originalNormals = NULL;
	m_radius = 0; // ���� m_mapKnearest�������뾶
	m_numOfNeighbors = 0; // ���������������ĵ�

	m_localRBFradius=0;//         �ֲ�����������뾶
	m_localOrdinate=NULL;//             �������ľֲ�����
	m_matrixNeighbor=NULL;//            ��ֵ����
	m_matrixNeighborInv=NULL;//         ��ֵ���������
	m_localIntegralPoints=NULL;  //���ֵ������
	m_localIntegralRadius=0;//      �ֲ����ְ뾶
	m_testFunction=NULL;//       ���Ժ�����ÿһ���������ֵ
	m_RBFfunction=NULL;//        ÿһ��Ϊ��ÿһ�������ľ����������ÿһ�����ֵ�ĺ���ֵ
	m_shapeFunction=NULL;//                ���������κ����ڸ��������ֵ
	m_determinentMetric=NULL;//   �ڲ����㴦���α�׼�Ĵ�������ʽ
	m_sqrtDetMetric = NULL;
	m_metricMatrix=NULL;//�ڲ����㴦���α�׼�ľ���
	m_metricMatrixInve=NULL;////�ڲ����㴦���α�׼�ľ������
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
	��_scalarProductVector = NULL;
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
	 * �������ܣ�  �ӿں���������㼯�˲�
	 * ����˵����
	 * int numOfPointSet                        �������
	 * double* pointSet							ԭʼ�㼯
	 * int kNearest                             k���������
	 * double radius                            k���������İ뾶
	 * double ebusainu                          �������Է�����ʱ���ֹͣ��׼
	 * 	int stopN                                �������Է�����ʱ��ѭ������
	 * maxIter   ˫�����ݶ��½������ϡ�����Է�����ʱ������������
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
	 * ������ϵͳ
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
	 * �ͷ��ڴ�
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
	*	�������ܣ� �ͷű����ڴ�
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
	m_loadVector=NULL;//��������
	if(m_localOrdinate!=NULL)
		delete[] m_loadVector;
	m_localOrdinate=NULL;//             �������ľֲ�����
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighbor=NULL;//            ��ֵ����
	if(m_matrixNeighbor!=NULL)
		delete[] m_matrixNeighbor;
	m_matrixNeighborInv=NULL;//        ��ֵ���������
	if(m_localIntegralPoints!=NULL)
		delete[] m_localIntegralPoints;
	m_localIntegralPoints=NULL;  //���ֵ������
	if(m_testFunction!=NULL)
		delete[] m_testFunction;
	m_testFunction=NULL;//       ���Ժ�����ÿһ���������ֵ
	if(m_RBFfunction!=NULL)
		delete[] m_RBFfunction;
	m_RBFfunction=NULL;//       ÿһ��Ϊ��ÿһ�������ľ����������ÿһ�����ֵ�ĺ���ֵ
	if(m_shapeFunction!=NULL)
		delete[] m_shapeFunction;
	m_shapeFunction=NULL;//                ���������κ����ڸ��������ֵ
	if(m_determinentMetric!=NULL)
		delete[] m_determinentMetric;
	m_determinentMetric=NULL;//   �ڲ����㴦���α�׼�Ĵ�������ʽ

	if ( m_sqrtDetMetric != NULL ){
		delete[] m_sqrtDetMetric;
	}
	m_sqrtDetMetric = NULL;

	if( ��_scalarProductVector != NULL )
		delete[] ��_scalarProductVector;
	��_scalarProductVector = NULL;

	if ( m_gradientProductVector != NULL ){
		delete[] m_gradientProductVector;
	}
	m_gradientProductVector = NULL;

	if(m_metricMatrix!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrix=NULL;//�ڲ����㴦���α�׼�ľ���
	if(m_metricMatrixInve!=NULL)
		delete[] m_metricMatrix;
	m_metricMatrixInve=NULL;////�ڲ����㴦���α�׼�ľ������
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
	*	�������ܣ�����ÿһ����� k����
	*/
	double radius=m_radius*m_radius;//��Ž��ڵ㵽�����ľ���
	std::multimap<double ,int> mapKnearest;//��ÿһ���㽨������vectorʱ������ʱmap
	int numOfKnearest;//ͳ��һ�������ʱ�����С
	if(m_mapKnearest.begin()!=m_mapKnearest.end()){
		m_mapKnearest.clear();
	}

	//////////////////////////////////////////////////////////////////////////
	//�������е㣬���Ұ��վ�����С�����Ŵ�ŵ�map��
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
		KnearestField fieldKnearest;//��ʱ�ṹ
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
	//// ���Դ���
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
	*	�������ܣ�����ÿһ����ķ���
	*/
	double centroidPosition[3];//����λ��
	double* localVariationVector;//���ĵ����λ�õ�����
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
		//�󷽲����
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
	// д�����ļ�
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
	//// д�������ļ�
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
	//// д�η����ļ�
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
	*	�������ܣ� ���Ʒ����׼ȷ����
	*/
	//����ȷ��zֵ������һ��ķ���
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

	int* parent=new int[m_numOfPoints];//��Ÿ��ڵ�

	std::multimap<double,int> valueKeySeries;//���ݵ�Ĵ���������map
	std::map<int,double> pointKeySeries;//���ݵ�������map

	//��ʼ��map
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
	*	�������ܣ� ������һ�㴦�Ļ��ְ뾶���������������ľֲ����꣬
	*  ����˵����
	*  int indexPointSet        ������	
	*  double m_localRBFradius;            ����������뾶
	*  double* m_localOrdinate             �������ľֲ�����
	*  double* m_localIntegralPoint        ���ֵ������
	*  double  m_localIntegralRadius       ������İ뾶
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
	*	д�ֲ������ļ�
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
	*	������Ժ����ڻ��ֵ�ĺ���ֵ���䵼��ֵ
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
	*	������ֵ�ľֲ�����
	*/
	//////////////////////////////////////////////////////////////////////////
	for (int i=0;i<m_numOfIntegralPoint;i++){
		m_localIntegralPoints[i*2]=m_ordinateAndWeightCoff[i*3]*m_localIntegralRadius;
		m_localIntegralPoints[i*2+1]=m_ordinateAndWeightCoff[i*3+1]*m_localIntegralRadius;		
	}
	//////////////////////////////////////////////////////////////////////////
	/*
	*	�����κ������䵼��  ѡ�� c��r���
	*/
	//////////////////////////////////////////////////////////////////////////
	double tempEXP,tempOneMinusEXP;// ������Ժ���ʱ���һЩ����
	double temp=-1;
	tempEXP=exp(temp);
	tempOneMinusEXP=1-tempEXP;
	double integralRaiusSquare, radiusCSquare;

	integralRaiusSquare=m_localIntegralRadius*m_localIntegralRadius;
	radiusCSquare=integralRaiusSquare;

	for(int i=0;i<m_numOfIntegralPoint;i++){
		double tempDistanceSquare;//���ֵ���ԭ��֮��ľ����ƽ��
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
	*	�������ܣ������ֵ����������
	*  ����˵����
	*  double* m_matrixNeighbor; ��ֵ����
	*  double* m_matrixNeighborInv; ��ֵ���������  
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
			// ����ÿһ�� �ڲ�ͬ�㲻ͬ��������ֵ

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
			///�������������
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
	 * �������ܣ� ������������йصľ���
	 * ����˵����
	 * int indexPointSet ������
	 * double* m_RBFfunction; ÿһ��Ϊ��ÿһ�������ľ����������ÿһ�����ֵ�ĺ���ֵ
	 */
	/////////////////////////////////////////////////////////////////////////
	/*
	 * �����ڴ�ռ�
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
	 * double* m_determinentMetric;//   �ڲ����㴦���α�׼�Ĵ�������ʽ
	 * double* m_metricMatrix;//�ڲ����㴦���α�׼�ľ���
	 * double* m_metricMatrixInve;////�ڲ����㴦���α�׼�ľ������
	 */
	////////////////////////////////////
	/*
	 * �����ڴ�ռ�
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

	double x_Der_One[ 3 ]; //// ��ά�ռ��λ������Ծֲ��ռ��һ������ĵ���
	double x_Der_Two[ 3 ]; //// ��ά�ռ��λ������Ծֲ��ռ�ڶ�������ĵ���

	for ( int i = 0; i < m_numOfIntegralPoint; i++ ){
		/// ��ʱ����
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
			///�������������
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
	 *  �������ܣ� ����RBF��������Ժ������ڻ�
	 *  ����˵����
	 *  int indexPointSet        ������
	 *  double*����_scalarProcuctVector  �ڲ����������пռ��У�n��RBF��������Ժ������ֵ�����
	 */
	
	if ( ��_scalarProductVector != NULL ){
		delete[] ��_scalarProductVector;
		��_scalarProductVector = NULL;
	}
	std::map< int, KnearestField >::iterator mapIterator = m_mapKnearest.find( indexPointSet );
	std::vector< int >::iterator vectorIterator = ( *mapIterator ).second.m_nearest.begin();
	m_numOfNeighbors = ( *mapIterator ).second.m_nearest.size()+1;

	��_scalarProductVector = new double[ m_numOfNeighbors ]; 

	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		��_scalarProductVector[ i ] = 0;
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){

			��_scalarProductVector[ i ] += m_RBFfunction[ i * m_numOfIntegralPoint + j ] * m_testFunction[ j ]
									* m_ordinateAndWeightCoff[ j * 3 + 2 ] * m_sqrtDetMetric[ j ];
		}
	}

}

void gradientFieldFilterRBF::ComputeVectorProductVector(int indexPointSet){
	/////////////////////////////////////////////////////
	// �����ɵ��������ݶ�
	/////////////////////////////////////////////////////
    /*
	 * ��ʱ����
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

	for ( int i = 0; i < m_numOfNeighbors; i++ ){ // ��m_numOfNeighbors��RBF������ÿ��������test function��һ������
		                                          // ��� ÿһ��RBF������ÿһ�����ֵ㴦Ҫ��һ���ݶȳ���
		for ( int j = 0; j < m_numOfIntegralPoint; j++ ){ // ��m_numOfIntegralPoint�����ֵ�
			
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
	 *����ÿһ��RBF������test function�������ϵĻ���
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
			double temp_Integral;// û�г���ÿһ�����ֵ�Ļ���ϵ��
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
	 *	�������ܣ� ����նȾ���ͺ�������
	 *  ��ÿһ��������indexPointSet��Ϊ�նȾ������һ�У�Ϊ������������һ��ֵ
	 *  ����˵���� 
	 *  int indexPointSet                    ������
	 *  double* m_loadVector                  ��������
	 *  StiffMatrix* m_stiffMatrix           �նȾ��������ݶȻ�������ֵ
	 *  StiffMatrix* m_massMatrix            �նȾ������ɱ�����������ֵ
	 */
	///////////////////////////////
	/*
	 *  ���� m_massMatrix
	 */
	///////////////////////////////
	    
	for ( int i = 0; i < m_numOfNeighbors; i++ ){
		double temp_one = 0;

		for ( int j = 0; j < m_numOfNeighbors; j++ ){
			temp_one += ��_scalarProductVector[ j ] * m_matrixNeighborInv[ j * m_numOfNeighbors + i ];
		}
		m_massMatrix[indexPointSet].m_elementStiffMatrix.push_back( temp_one );
	}
	///////////////////////////////
	/*
	 *  ���� m_stiffMatrix
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
	 *	�������ܣ� ��װ�նȾ���
	 *  ���㷽�����ұߺ͸նȾ����ȼ����ұߣ��ٸĸնȾ��� �����������͸նȾ���
	 */
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	/*
	 * ���㷽�����ұ� ��������
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
	 * ���㷽������� ����նȾ���
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
	*	pbcg ��ϡ�����Է�����
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



#include "mainDoc.h"
#include <math.h>
#include "GenerateNoise.h"
#include "Vectors.h"
#include "qmessagebox.h"
#include <QtGui/QApplication>
mainDoc::mainDoc()
{
	m_triangleTopology = NULL;
	m_triangleGeometry = NULL;
	m_triangleTexture = NULL;
	m_triangleNormal = NULL;
	m_triangleGeometryWithNoise = NULL;
	m_triangleGeometryAfterFiltering = NULL;
	m_triangleNormalWithNoise=NULL;
	filter=NULL;
//	m_triangleNormalAfterFiltering=null
	//m_numOfVertex = tVertex.size();	




	m_noiseCoffecient = 0.1;
	
}

mainDoc::~mainDoc()
{
	if ( m_triangleTopology != NULL){
		delete[] m_triangleTopology;
		m_triangleTopology = NULL;
	}
	if( m_triangleGeometry != NULL){
		delete[] m_triangleGeometry;
		m_triangleGeometry = NULL;
	}
	if ( m_triangleTexture != NULL ){
		delete[] m_triangleTexture;
		m_triangleTexture = NULL;
	}
	if( m_triangleNormal != NULL ){
		delete[] m_triangleNormal;
		m_triangleNormal = NULL;
	}

	if ( m_triangleGeometryWithNoise != NULL ){
		delete[] m_triangleGeometryWithNoise;
		m_triangleGeometryWithNoise = NULL;
	}

	if ( m_triangleGeometryAfterFiltering != NULL ){
		delete[] m_triangleGeometryAfterFiltering;
		m_triangleGeometryAfterFiltering = NULL;
	}
	if(filter)
	{
		delete filter;
	}
	
	
}

bool mainDoc::setTriangleGeometry(std::vector<Point> tVertex){

	if( m_triangleGeometry != NULL){
		delete[] m_triangleGeometry;
		m_triangleGeometry = NULL;
	}
	
	m_numOfVertex = tVertex.size();	
	//m_triangleGeometryWithNoise=new float[m_numOfPoints*3];
	//m_originalOriginalPointSets=new float[m_numOfPoints*3];
	m_triangleGeometry = new float[ m_numOfVertex *3 ];

	std::vector<Point>::iterator tempPointIterator=tVertex.begin();

	for(int i = 0; tempPointIterator!=tVertex.end();tempPointIterator++ , i++){

		m_triangleGeometry[i*3] = (*tempPointIterator).vertexPosition[0] * 1000;
		m_triangleGeometry[i*3+1] = (*tempPointIterator).vertexPosition[1] * 1000 ;
		m_triangleGeometry[i*3+2] = (*tempPointIterator).vertexPosition[2] * 1000;
	}
	calculateDataBox();

	return true;
}

bool mainDoc::setTriangleTopology(std::vector<Triangle> tTopology){
	
	m_numOfTriangle = tTopology.size();
	if( m_triangleTopology != NULL ){
		delete[] m_triangleTopology;
		m_triangleTopology = NULL;
	}
	m_triangleTopology = new int[ m_numOfTriangle * 3 ];
    
	std::vector<Triangle>::iterator tempTriangleIterator=tTopology.begin();
	for(int i=0;tempTriangleIterator!=tTopology.end();tempTriangleIterator++,i++){
	
		m_triangleTopology[i*3]=(*tempTriangleIterator).vertexIndex[0];
		m_triangleTopology[i*3+1]=(*tempTriangleIterator).vertexIndex[1];
		m_triangleTopology[i*3+2]=(*tempTriangleIterator).vertexIndex[2];

	}

	calculateTriangleNormal();
	
	return true;
}
bool mainDoc::setTriangleTexture(std::vector<Point> tTexture){
	return true;

}

bool mainDoc::calculateTriangleNormal(){
	
	if( m_triangleNormal != NULL ){
		delete[] m_triangleNormal;
		m_triangleNormal = NULL;
	}

	m_triangleNormal= new float[m_numOfVertex*3];
	
	for (int i = 0; i < m_numOfVertex*3; i++) {
		m_triangleNormal[i]=0;
	}

	// Calculate normals.
	for (int i = 0; i < m_numOfTriangle; i++) {
		float vec1[3], vec2[3], normal[3];
		int id0, id1, id2;
		id0 = m_triangleTopology[i*3];
		id1 = m_triangleTopology[i*3+1];
		id2 = m_triangleTopology[i*3+2];
		vec1[0] = m_triangleGeometry[id1*3]- m_triangleGeometry[id0*3];
		vec1[1] = m_triangleGeometry[id1*3+1] - m_triangleGeometry[id0*3+1];
		vec1[2] = m_triangleGeometry[id1*3+2] - m_triangleGeometry[id0*3+2];
		vec2[0] = m_triangleGeometry[id2*3] - m_triangleGeometry[id0*3];
		vec2[1] = m_triangleGeometry[id2*3+1] - m_triangleGeometry[id0*3+1];
		vec2[2] = m_triangleGeometry[id2*3+2] - m_triangleGeometry[id0*3+2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
	
		m_triangleNormal[id0*3+0] += normal[0];
		m_triangleNormal[id0*3+1] += normal[1];
		m_triangleNormal[id0*3+2] += normal[2];
		m_triangleNormal[id1*3+0] += normal[0];
		m_triangleNormal[id1*3+1] += normal[1];
		m_triangleNormal[id1*3+2] += normal[2];
		m_triangleNormal[id2*3+0] += normal[0];
		m_triangleNormal[id2*3+1] += normal[1];
		m_triangleNormal[id2*3+2] += normal[2];
	}

	// Normalize normals.
	for (int i = 0; i < m_numOfVertex; i++) {
		
		float length = sqrt(m_triangleNormal[i*3+0]*m_triangleNormal[i*3+0] + m_triangleNormal[i*3+1]*m_triangleNormal[i*3+1] + m_triangleNormal[i*3+2]*m_triangleNormal[i*3+2]);
		if(length<0.000001){
				m_triangleNormal[i*3] = 0;		m_triangleNormal[i*3+1] = 0;		m_triangleNormal[i*3+2] = 1;
		}
		else{
		m_triangleNormal[i*3] /= length;
		m_triangleNormal[i*3+1] /= length;
		m_triangleNormal[i*3+2] /= length;
		}
	}
	return true;
}
bool mainDoc::calculatePointNormal(){
	return true;
}
bool mainDoc::calculateDataBox(){


	box_x_max=box_x_min=m_triangleGeometry[0];
	box_y_max=box_y_min=m_triangleGeometry[1];
	box_z_max=box_z_min=m_triangleGeometry[2];

	for(int i=0;i<m_numOfVertex;i++){
		if(m_triangleGeometry[i*3]>box_x_max)
			box_x_max=m_triangleGeometry[i*3];
		else{
			if(m_triangleGeometry[i*3]<box_x_min)
				box_x_min=m_triangleGeometry[i*3];
		}
		if(m_triangleGeometry[i*3+1]>box_y_max)
			box_y_max=m_triangleGeometry[i*3+1];
		else{
			if(m_triangleGeometry[i*3+1]<box_y_min)
				box_y_min=m_triangleGeometry[i*3+1];
		}
		if(m_triangleGeometry[i*3+2]>box_z_max)
			box_z_max=m_triangleGeometry[i*3+2];
		else{
			if(m_triangleGeometry[i*3+2]<box_z_min)
				box_z_min=m_triangleGeometry[i*3+2];
		}
	}
	return true;
}

bool mainDoc::addPositionGaussianNoise(){
	float m_meanLength = 1.0;
	int m_nNormals = m_numOfVertex;
	m_meanLength = calculateMeanLength();
	
	//////////////////////////////////////////////////////////////////////////
	/*
	 *	加法向噪声
	 */
	///////////////////////////////////////////////////////

	m_triangleGeometryWithNoise=new float[m_numOfVertex*3];
	//m_triangleGeometry=new float[m_numOfVertex*3];

	
	float noise; 
	long idum = 1; 
	for ( int i = 0; i < m_numOfVertex; i++ )
	{
		noise = gasdev( &idum );	
		for ( int j = 0; j < 3; j++ )
		{
			//
			//m_triangleGeometryWithNoise[i*3+j] = m_triangleGeometry[i*3+j]+noise*m_triangleNormal[i*3+j]*m_noiseCoffecient*m_meanLength;
			m_triangleGeometryWithNoise[i*3+j] = m_triangleGeometry[i*3+j]+noise*m_triangleNormal[i*3+j]*m_noiseCoffecient*m_meanLength;
		    m_triangleGeometry[i*3+j]=m_triangleGeometryWithNoise[i*3+j];
		}
	}


	//if((fpout = fopen("bunny_noise.", "r")) == NULL)
	//{		
	//	exit(1);
	//}
	//else
	//{
	//	while(!feof(fpout)){			
	//		
 //           Triangle tempTriangle;
	//		fscanf( fpout, "%d %d %d", &tempTriangle.vertexIndex[0], &tempTriangle.vertexIndex[1], &tempTriangle.vertexIndex[2]);
	//		tempTriangleSets.push_back( tempTriangle );	
	//	}
	//	fclose(fpout);
	//}




	//////////////////////////////////////////////////////////////////////////
	/*
	 *	计算加噪声后的法向
	 */
	///////////////////////////////////////////////////////
		
	//if ( m_triangleNormalWithNoise != NULL ){
	//delete[] m_triangleNormalWithNoise;
	//m_triangleNormalWithNoise = NULL;
	//}
	//m_triangleNormalWithNoise = new float[ m_numOfVertex * 3 ];

	//for ( int i = 0; i < m_numOfVertex * 3; i++ ){
	//	m_triangleNormalWithNoise[ i ] = 0;
	//}

	//for (  int i = 0; i < m_numOfTriangle; i++ ){
	//	VECTOR3D vec1, vec2, normal;
	//	int id0, id1, id2;
	//	id0 = m_triangleTopology[ i * 3 ];
	//	id1 = m_triangleTopology[ i * 3 + 1 ];
	//	id2 = m_triangleTopology[ i * 3 + 2 ];

	//	vec1[ 0 ] = m_triangleGeometry[ id1 * 3 ] - m_triangleGeometry[ id0 * 3 ];
	//	vec1[ 1 ] = m_triangleGeometry[ id1 * 3 + 1 ] - m_triangleGeometry[ id0 * 3 + 1 ];
	//	vec1[ 2 ] = m_triangleGeometry[ id1 * 3 + 2 ] - m_triangleGeometry[ id0 * 3 + 2 ];

	//	vec2[ 0 ] = m_triangleGeometry[ id2 * 3 ] - m_triangleGeometry[ id0 * 3 ];
	//	vec2[ 1 ] = m_triangleGeometry[ id2 * 3 + 1 ] - m_triangleGeometry[ id0 * 3 + 1 ];
	//	vec2[ 2 ] = m_triangleGeometry[ id2 * 3 + 2 ] - m_triangleGeometry[ id0 * 3 + 2 ];

	//	normal[ 0 ] = vec1[ 2 ] * vec2[ 1 ] - vec1[ 1 ] * vec2[ 2 ];
	//	normal[ 1 ] = vec1[ 0 ] * vec2[ 2 ] - vec1[ 2 ] * vec2[ 0 ];
	//	normal[ 2 ] = vec1[ 1 ] * vec2[ 0 ] - vec1[ 0 ] * vec2[ 1 ];

	//	m_triangleNormalWithNoise[ id0 * 3 + 0 ] += normal[ 0 ];
	//	m_triangleNormalWithNoise[ id0 * 3 + 1 ] += normal[ 1 ];
	//	m_triangleNormalWithNoise[ id0 * 3 + 2 ] += normal[ 2 ];

	//	m_triangleNormalWithNoise[ id1 * 3 + 0 ] += normal[ 0 ];
	//	m_triangleNormalWithNoise[ id1 * 3 + 1 ] += normal[ 1 ];
	//	m_triangleNormalWithNoise[ id1 * 3 + 2 ] += normal[ 2 ];

	//	m_triangleNormalWithNoise[ id2 * 3 + 0 ] += normal[ 0 ];
	//	m_triangleNormalWithNoise[ id2 * 3 + 1 ] += normal[ 1 ];
	//	m_triangleNormalWithNoise[ id2 * 3 + 2 ] += normal[ 2 ];
	//}

	////// Normalize normals
	//for ( int i = 0; i < m_numOfVertex; i++ ){
	//	float length = sqrt( m_triangleNormalWithNoise[ i * 3 + 0 ] * m_triangleNormalWithNoise[ i * 3 + 0 ] 
	//						+ m_triangleNormalWithNoise[ i * 3 + 1 ] * m_triangleNormalWithNoise[ i * 3 + 1 ]
	//						+ m_triangleNormalWithNoise[ i * 3 + 2 ] * m_triangleNormalWithNoise[ i * 3 + 2 ] );
	//	m_triangleNormalWithNoise[ i * 3 + 0 ] /= length;
	//	m_triangleNormalWithNoise[ i * 3 + 1 ] /= length;
	//	m_triangleNormalWithNoise[ i * 3 + 2 ] /= length;

	//}











	for ( int i = 0; i < m_numOfVertex * 3; i++ ){
		m_triangleNormal[ i ] = 0;
	}

	for (  int i = 0; i < m_numOfTriangle; i++ )
	{
		VECTOR3D vec1, vec2, normal;
		int id0, id1, id2;
		id0 = m_triangleTopology[ i * 3 ];
		id1 = m_triangleTopology[ i * 3 + 1 ];
		id2 = m_triangleTopology[ i * 3 + 2 ];

		vec1[ 0 ] = m_triangleGeometry[ id1 * 3 ] - m_triangleGeometry[ id0 * 3 ];
		vec1[ 1 ] = m_triangleGeometry[ id1 * 3 + 1 ] - m_triangleGeometry[ id0 * 3 + 1 ];
		vec1[ 2 ] = m_triangleGeometry[ id1 * 3 + 2 ] - m_triangleGeometry[ id0 * 3 + 2 ];

		vec2[ 0 ] = m_triangleGeometry[ id2 * 3 ] - m_triangleGeometry[ id0 * 3 ];
		vec2[ 1 ] = m_triangleGeometry[ id2 * 3 + 1 ] - m_triangleGeometry[ id0 * 3 + 1 ];
		vec2[ 2 ] = m_triangleGeometry[ id2 * 3 + 2 ] - m_triangleGeometry[ id0 * 3 + 2 ];

		normal[ 0 ] = vec1[ 2 ] * vec2[ 1 ] - vec1[ 1 ] * vec2[ 2 ];
		normal[ 1 ] = vec1[ 0 ] * vec2[ 2 ] - vec1[ 2 ] * vec2[ 0 ];
		normal[ 2 ] = vec1[ 1 ] * vec2[ 0 ] - vec1[ 0 ] * vec2[ 1 ];

		m_triangleNormal[ id0 * 3 + 0 ] += normal[ 0 ];
		m_triangleNormal[ id0 * 3 + 1 ] += normal[ 1 ];
		m_triangleNormal[ id0 * 3 + 2 ] += normal[ 2 ];

		m_triangleNormal[ id1 * 3 + 0 ] += normal[ 0 ];
		m_triangleNormal[ id1 * 3 + 1 ] += normal[ 1 ];
		m_triangleNormal[ id1 * 3 + 2 ] += normal[ 2 ];

		m_triangleNormal[ id2 * 3 + 0 ] += normal[ 0 ];
		m_triangleNormal[ id2 * 3 + 1 ] += normal[ 1 ];
		m_triangleNormal[ id2 * 3 + 2 ] += normal[ 2 ];
	}

	//// Normalize normals
	for ( int i = 0; i < m_numOfVertex; i++ ){
		float length = sqrt( m_triangleNormal[ i * 3 + 0 ] * m_triangleNormal[ i * 3 + 0 ]+ m_triangleNormal[ i * 3 + 1 ] * m_triangleNormal[ i * 3 + 1 ]+ m_triangleNormal[ i * 3 + 2 ] * m_triangleNormal[ i * 3 + 2 ] );
		m_triangleNormal[ i * 3 + 0 ] /= length;
		m_triangleNormal[ i * 3 + 1 ] /= length;
		m_triangleNormal[ i * 3 + 2 ] /= length;

	}
	

	
	return true;
}
float mainDoc::calculateMeanLength(){

	float _mean_length = 0;
	for ( int i = 0; i < m_numOfTriangle; i++ ){
		_mean_length += sqrt( ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 1 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 1 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 2 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 2 ] ) );
		_mean_length += sqrt( ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2] * 3 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2] * 3 + 1 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2] * 3 + 1 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2] * 3 + 2 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 2] * 3 + 2 ] ) );
		_mean_length += sqrt( ( m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 1 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 + 1 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 1 ] )
							+ ( m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 2 ] ) * (m_triangleGeometry[ m_triangleTopology[ i * 3 + 2 ] * 3 + 2 ] - m_triangleGeometry[ m_triangleTopology[ i * 3 + 1] * 3 + 2 ] ) );

	}

	_mean_length /=( 3 * m_numOfTriangle );

	return _mean_length;
}
bool mainDoc::RBFilter()
{

	float m_meanlength;	
	m_meanlength=calculateMeanLength();
	
	
//	m_dimensionColor=3;
	//filter->GetGradientFieldFilter(m_numOfVertex,m_triangleGeometry,m_kNearest,m_radius,m_tol,m_itol,m_itmax,m_curveThreshold); 
//filter->GetGradientFieldFilter(m_numOfVertex,m_triangleGeometry,6,0.5,0.5,4,6,0.9); 
	QMessageBox::information(NULL, "aaa ", "hxy ");
	//m_meshLessFilter.GetMeshlessFilter(m_numOfPoints,m_originalPointSets,m_K,m_radius*m_meanlength,m_iterativeTimes,m_timeStep,m_loadConstant,m_positionError*m_meanlength,m_stopN,m_maxIter);
	
	
	
	//
	//for(int i=0;i<m_numOfPoints*3;i++){
	//	m_triangleGeometry[i]=(float)m_resultPointSet[i];	
	//	
	//}
	//
	//int m_nNormals = m_numOfPoints;
	//m_normals=new float[m_nNormals*3];

	//// Set all normals to 0.
	//for (int i = 0; i < m_nNormals*3; i++) {
	//	m_normals[i]=0;
	//}

	//// Calculate normals.
	//for (int i = 0; i < m_numOfTriangles; i++) {
	//	VECTOR3D vec1, vec2, normal;
	//	int id0, id1, id2;
	//	id0 = m_triangles[i*3];
	//	id1 = m_triangles[i*3+1];
	//	id2 = m_triangles[i*3+2];
	//	vec1[0] = m_pointSets[id1*3]- m_pointSets[id0*3];
	//	vec1[1] = m_pointSets[id1*3+1] - m_pointSets[id0*3+1];
	//	vec1[2] = m_pointSets[id1*3+2] - m_pointSets[id0*3+2];
	//	vec2[0] = m_pointSets[id2*3] - m_pointSets[id0*3];
	//	vec2[1] = m_pointSets[id2*3+1] - m_pointSets[id0*3+1];
	//	vec2[2] = m_pointSets[id2*3+2] - m_pointSets[id0*3+2];
	//	normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
	//	normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
	//	normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
	//	m_normals[id0*3+0] += normal[0];
	//	m_normals[id0*3+1] += normal[1];
	//	m_normals[id0*3+2] += normal[2];
	//	m_normals[id1*3+0] += normal[0];
	//	m_normals[id1*3+1] += normal[1];
	//	m_normals[id1*3+2] += normal[2];
	//	m_normals[id2*3+0] += normal[0];
	//	m_normals[id2*3+1] += normal[1];
	//	m_normals[id2*3+2] += normal[2];
	//}

	//// Normalize normals.
	//for (int i = 0; i < m_nNormals; i++) {
	//	float length = sqrt(m_normals[i*3+0]*m_normals[i*3+0] + m_normals[i*3+1]*m_normals[i*3+1] + m_normals[i*3+2]*m_normals[i*3+2]);
	//	m_normals[i*3] /= length;
	//	m_normals[i*3+1] /= length;
	//	m_normals[i*3+2] /= length;
	//}



	//QString filename_pw = "E:\\resultPoints.txt";

	//FILE *fpout;
	//if((fpout=fopen(filename_pw,"w"))==NULL){
	//	int dkjkd;

	//	return;
	//}
	//else{
	//	for(int i=0;i<m_numOfPoints;i++){
	//		fprintf(fpout,"%f %f %f\n",m_triangleGeometry[i*3],m_triangleGeometry[i*3+1],m_triangleGeometry[i*3+2]);		

	//	}
	//	fclose(fpout);
 //       
	//}


	//
	//
	
	filter->DeleteGradientFieldFilter();
	
	return true;
}


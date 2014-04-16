#include "mainGLView.h"
#include <string>
#include "NormalStruct.h"
#include <iostream>
#include <fstream>
#include "mainDoc.h"
#include <qmessagebox.h>



mainGLView::mainGLView(void)
{
	//  Qt initialization]
	//makeCurrent();
	setFocusPolicy( Qt::WheelFocus);
	setFocus(Qt::MouseFocusReason);
	setMinimumSize(512,512 );
	setAttribute(Qt::WA_NoSystemBackground);
	setFormat( QGLFormat(  QGL::DepthBuffer | QGL::DoubleBuffer) );
	///////////////////////////
	char* a = getenv("PWD");
	currOpenFileDir = new QString(a ? a : "/");
	//m_pointCloudsGeometry =NULL;
	m_mainDoc = new mainDoc;

	///////////////
	// initiate private variable

	m_cameraPosition.pointVector[0] = 0.0;
	m_cameraPosition.pointVector[1] = 0.0;
	m_cameraPosition.pointVector[2] = 200.0;
	m_focusPosition.pointVector[0] = 0.0;
	m_focusPosition.pointVector[1] = 0.0;
	m_focusPosition.pointVector[2] = 0.0;
	m_upDirection.pointVector[0] = 0.0;
	m_upDirection.pointVector[1] = 1.0;
	m_upDirection.pointVector[2] = 0.0;
	m_fovy = 90;
	m_aspect = 1;
	m_near = 1;
	m_fast = 400;	
	rotationX = 0.0;
	rotationY = 0.0;
	scaleFactor = 1.0;
	if_triangle_geometry = false;
	if_triangle_topology = false;
	if_triangle_rendering = false;
}

mainGLView::~mainGLView(void)
{
	delete m_mainDoc;


}


bool mainGLView::openPointCloudsGeometry(QString fileName){

	//std::string t_filename = fileName.toAscii().data();
	//FILE *fpout;

	//std::vector<Point> tempPointSets;

	//if((fpout = fopen(t_filename.c_str(), "r")) == NULL)
	//{		
	//	exit(1);
	//}
	//else
	//{
	//	while(!feof(fpout)){
	//		Point tempPoint;
	//		float temp1,temp2;
	//		char tempChar[10];

	//		fscanf(fpout,"%f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2]);
	//	    tempPointSets.push_back(tempPoint);
	//	}
	//	fclose(fpout);
	//}

	//m_numOfPointClouds = tempPointSets.size();

	//if (m_pointCloudsGeometry !=NULL)
	//	delete[] m_pointCloudsGeometry;

	//m_pointCloudsGeometry = new float[m_numOfPointClouds*3];

	//std::vector<Point>::iterator tempPointIterator=tempPointSets.begin();

	//for(int i = 0; tempPointIterator!=tempPointSets.end();tempPointIterator++ , i++){

	//	m_pointCloudsGeometry[i*3] = (*tempPointIterator).vertexPosition[0];
	//	m_pointCloudsGeometry[i*3+1] = (*tempPointIterator).vertexPosition[1];
	//	m_pointCloudsGeometry[i*3+2] = (*tempPointIterator).vertexPosition[2];

	//}

	//tempPointSets.clear();

	return true;
}
bool mainGLView::openPointCloudsTexture( QString fileName ){
	return true;
}
bool mainGLView::openTriangleGeometry(QString fileName ){
	
	std::string t_filename = fileName.toAscii().data();
	FILE *fpout;

	std::vector<Point> tempPointSets;

	if((fpout = fopen(t_filename.c_str(), "r")) == NULL)
	{		
		exit(1);
	}
	else
	{
		while(!feof(fpout)){
			
			
			Point tempPoint;
			float temp1,temp2;
			char tempChar[10];

			fscanf(fpout,"%f %f %f %f %f",&tempPoint.vertexPosition[0],&tempPoint.vertexPosition[1],&tempPoint.vertexPosition[2], &temp1, &temp2);
		    tempPointSets.push_back(tempPoint);
		}
		fclose(fpout);
	}
  //  m_mainDoc = new mainDoc();
	m_mainDoc->setTriangleGeometry( tempPointSets );
	tempPointSets.clear();
	if_triangle_geometry = true;
	if( if_triangle_geometry && if_triangle_topology )
		if_triangle_rendering = true;

	return true;

}
bool mainGLView::openTriangleTopology( QString fileName ){
	std::string t_filename = fileName.toAscii().data();
	FILE *fpout;

	std::vector<Triangle> tempTriangleSets;

	if((fpout = fopen(t_filename.c_str(), "r")) == NULL)
	{		
		exit(1);
	}
	else
	{
		while(!feof(fpout)){			
			
            Triangle tempTriangle;
			fscanf( fpout, "%d %d %d", &tempTriangle.vertexIndex[0], &tempTriangle.vertexIndex[1], &tempTriangle.vertexIndex[2]);
			tempTriangleSets.push_back( tempTriangle );	
		}
		fclose(fpout);
	}

	m_mainDoc->setTriangleTopology(tempTriangleSets);

	tempTriangleSets.clear();
	if_triangle_topology = true;
	if( if_triangle_geometry && if_triangle_topology )
		if_triangle_rendering = true;

	this->updateGL();


	//QMessageBox message(QMessageBox::NoIcon,"Show Qt","Do you want to show Qt dialog?", QMessageBox::Yes | QMessageBox::No, NULL);  
	//if(message.exec() == QMessageBox::Yes)  
	//{  
	//	QMessageBox::aboutQt(NULL,"About Qt");  
	//}  
	return true;
}

bool mainGLView::openTriangleTexture( QString fileName ){
	return true;
}

void mainGLView::initializeGL(){

	GLfloat light0_position[] = { 0, 0, 0, 1.0 };	


	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glShadeModel( GL_SMOOTH );
	glEnable( GL_DEPTH_TEST );
	glEnable( GL_CULL_FACE );

	GLfloat mat_ambient[]={0.1,0.5,0.8,0.15};
	GLfloat mat_diffuse[] = { 0.1, 0.5, 0.8, 1.0 };
	GLfloat mat_specular[]={ 1.0,1.0,1.0,1.0 };
	GLfloat mat_shiness[]={ 15.0 };

	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shiness);

	GLfloat global_ambient[] = { AMBIENT, AMBIENT, AMBIENT, 1 };
	GLfloat light0_ambient[] = { 0, 0, 0, 1 };
	GLfloat light0_diffuse[] = { DIFFUSE, DIFFUSE, DIFFUSE, 1.0 };
	GLfloat light0_specular[] = { SPECULAR, SPECULAR, SPECULAR, 1.0};
    
	glLightfv( GL_LIGHT0, GL_AMBIENT, light0_ambient );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, light0_diffuse );
	glLightfv( GL_LIGHT0, GL_SPECULAR, light0_ambient );
//	glLightfv( GL_LIGHT0, GL_POSITION, light0_position );

	glLightModelfv( GL_LIGHT_MODEL_AMBIENT, global_ambient );
	glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );

	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );

	glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
	glEnable( GL_COLOR_MATERIAL );	


}

void mainGLView::resizeGL( int width, int height){

	glViewport( 0, 0, width, height);	
	glMatrixMode(GL_PROJECTION);	
	glLoadIdentity();
	gluPerspective( m_fovy, m_aspect, m_near, m_fast );
	glMatrixMode( GL_MODELVIEW );

}
void mainGLView::paintGL(){

	// Make sure we're rendering in the correct GL instance	
    
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	draw();		
}

void mainGLView::draw(){

	if (!if_triangle_rendering )
		return;

	glMatrixMode( GL_MODELVIEW );

	glLoadIdentity();

	GLfloat position[]={-0,-0,-20000,1};


 	glLightfv(GL_LIGHT0,GL_POSITION,position);


	gluLookAt(m_cameraPosition.pointVector[0], m_cameraPosition.pointVector[1], m_cameraPosition.pointVector[2], 
		m_focusPosition.pointVector[0], m_focusPosition.pointVector[1], m_focusPosition.pointVector[2], 
		m_upDirection.pointVector[0], m_upDirection.pointVector[1], m_upDirection.pointVector[2] );
	glRotatef( rotationX, 1.0, 0.0, 0.0 );
	glRotatef( rotationY, 0.0, 1.0, 0.0 );
	
	glScalef( scaleFactor, scaleFactor, scaleFactor) ;

	glTranslated( -(m_mainDoc->box_x_max + m_mainDoc->box_x_min) / 2, 
		- ( m_mainDoc->box_y_max + m_mainDoc->box_y_min ) / 2,
		- ( m_mainDoc->box_z_max + m_mainDoc->box_z_min ) / 2 );



	int t_numOfVertex = m_mainDoc->m_numOfVertex;

	int t_numOfTriangle = m_mainDoc->m_numOfTriangle;

	float *t_vertex = m_mainDoc->m_triangleGeometry;
	int *t_triangle = m_mainDoc->m_triangleTopology;
	float *t_normal = m_mainDoc->m_triangleNormal;
	
	for (int i = 0; i < t_numOfTriangle; i++){

		glBegin( GL_TRIANGLES );
		{
			int j = t_triangle[ i * 3 ];
			glNormal3f( t_normal[j *3], t_normal[ j * 3 + 1], t_normal[j * 3  +2] );
			//glColor3f( 1.0, 0.8, 0.50);
			glVertex3f( t_vertex[ j *3], t_vertex[j *3 + 1], t_vertex[j *3 + 2] );
			//glVertex3f(-0.10, -0.10, -5);

			j = t_triangle[ i * 3 + 1];
			glNormal3f( t_normal[ j *3], t_normal[ j * 3 + 1], t_normal[ j * 3  +2] );
			//glColor3f( 0.8, 1.0, 0.50);
			glVertex3f( t_vertex[ j *3], t_vertex[ j *3 + 1], t_vertex[ j *3 + 2] );
		//	glVertex3f(0.10, -0.10, -5);

			j = t_triangle[ i * 3 + 2 ];
			glNormal3f( t_normal[j *3], t_normal[ j * 3 + 1], t_normal[ j * 3  +2] );
			//glColor3f( 0.3, 0.4, 1.0 );
			glVertex3f( t_vertex[ j *3], t_vertex[ j *3 + 1], t_vertex[ j *3 + 2] );
		//	glVertex3f( -0.10, 0.10, -5 );
		}

		glEnd();
	}	
}
bool mainGLView::PositionGaussianNoise(){

	m_mainDoc->addPositionGaussianNoise();


	updateGL();

//	QMessageBox::information(NULL, "aaa ", "hxy ");
	return true;
}
bool mainGLView::gradientFieldFilterRBFilter()
{

	m_mainDoc->RBFilter();

	
	updateGL();

		//QMessageBox::information(NULL, "aaa ", "hxy ");
	return true;
}
void mainGLView::mousePressEvent(QMouseEvent *event){

	lastPos = event->pos();
}

void mainGLView::mouseMoveEvent( QMouseEvent *event ){
	GLfloat dx = GLfloat( event->x() - lastPos.x())/width();
	GLfloat dy = GLfloat( event->y() - lastPos.y())/height();

	if( event->buttons() & Qt::LeftButton ){
		rotationX += 180 * dy;
		rotationY += 180 * dx;
		updateGL();
	}
	//else( event->buttons() & Qt::RightButton ){
	//	updateGL();
	//}
	lastPos = event->pos();
}

void mainGLView::wheelEvent( QWheelEvent *event){

	scaleFactor = scaleFactor * ( 1 + 0.05 * event->delta() / 120.0 );
	updateGL();

}
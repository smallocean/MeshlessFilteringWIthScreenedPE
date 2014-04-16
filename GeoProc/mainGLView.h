#ifndef MAINGLVIEW_H
#define MAINGLVIEW_H

#include <QGLWidget>
#include <QFileDialog>
#include <QMouseEvent>
#include <QWHeelEvent>
#include <QString>
#include <GL/glut.h>
#include "NormalStruct.h"
#include "mainDoc.h"

static const float AMBIENT = 0.25;
static const float DIFFUSE = 0.85;
static const float SPECULAR = 0.85;

class mainGLView :
	public QGLWidget
{
	Q_OBJECT
public:
	mainGLView(void);
	~mainGLView(void);

public slots:
	
	
	bool openPointCloudsGeometry(QString);
	bool openPointCloudsTexture( QString );
	bool openTriangleGeometry( QString );
	bool openTriangleTopology( QString );
	bool openTriangleTexture( QString );
	bool PositionGaussianNoise(); //¼Ó¸ßË¹ÔëÉù

	bool gradientFieldFilterRBFilter();  //RBFÂË²¨
protected:
	void initializeGL();
	void resizeGL( int, int );
	void paintGL();
	void draw();
	void mousePressEvent( QMouseEvent *event);
	void mouseMoveEvent( QMouseEvent *event );
	void wheelEvent(QWheelEvent *event);

protected:
	QString* currOpenFileDir;
	float* m_pointCloudsGeometry;
	float* m_pointCloudsTexture;
	float* m_triangleGeometry;
	int* m_triangleTopology;
	float* m_triangelTexture;

	mainDoc *m_mainDoc;

	int m_numOfPointClouds;
	int m_numOfVertex;
	int m_numOfTriangle;

	POINTVECTOR3D m_cameraPosition;
	POINTVECTOR3D m_focusPosition;
	POINTVECTOR3D m_upDirection;
	float m_fovy;
	float m_aspect;
	float m_near;
	float m_fast;
	GLfloat rotationX;
	GLfloat rotationY;
	GLfloat rotationZ;
	GLfloat scaleFactor;
	QPoint lastPos;
	bool if_triangle_rendering;
	bool if_triangle_geometry;
	bool if_triangle_topology;

private:
	QAction *openActionVertex;
	QAction *openActionTopology;

};
#endif 

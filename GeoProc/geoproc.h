#ifndef GEOPROC_H
#define GEOPROC_H

#include <QtGui/QMainWindow>
#include "mainGLView.h"
#include "ui_geoproc.h"
#include<QGridLayout>

class GeoProc : public QMainWindow
{
	Q_OBJECT

public:
	GeoProc( int argc = 0, char **argv = NULL );
	~GeoProc();


public slots:
	
	void openTriangleGeometryName();
	void openTriangleTopologyName();
	void addGaussionNoise();
//	void openTriangleTextureName();

signals:
	void openTriangleGeometrySig(QString);
	void openTriangleTopologySig(QString);
	void addGaussionNoiseSig(QString);
//	void openTriangleTextureSig(QString);


protected:
	void createActions();
	void createMenus();
	mainDoc* pdc;
private:
	//Ui::GeoProcClass ui;
	QString* currOpenFileDir;
	mainGLView *_mainGLView;  
	QWidget* centralWidget;
	QAction *_open_Vertex_Action;
	QAction *_open_Topology_Action;
	QAction *_add_noise;
	QAction * _pro_noise;
	QMenu *fileMenu;
	 

};

#endif // GEOPROC_H

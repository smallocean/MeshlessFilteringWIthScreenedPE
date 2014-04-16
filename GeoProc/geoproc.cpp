#include "geoproc.h"
#include <QGLWidget>
#include <GenerateNoise.h>
#include <mainDoc.h>
#include <qmessagebox.h>
GeoProc::GeoProc( int argc, char** argv )
	
{
	//ui.setupUi(this);
	char* a = getenv("PWD");
	currOpenFileDir = new QString(a ? a : "/");

	centralWidget = new QWidget;
	this->setCentralWidget( centralWidget );

	_mainGLView = new mainGLView;
	_mainGLView->setMinimumSize( 512,512 );

	QGridLayout *centralLayout = new QGridLayout;

	centralLayout->addWidget( _mainGLView, 0, 0 );
	centralLayout->setSizeConstraint(QLayout::SetMinimumSize);
    
	centralWidget->setLayout( centralLayout );
    pdc=NULL;
	//setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
    // setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);
    createActions();
	createMenus();
	setWindowTitle("Geometry Processing");
	
	move(100,100);

}

GeoProc::~GeoProc()
{

}

void GeoProc::createActions(){
	_open_Vertex_Action = new QAction( tr("&OpenV"), this );
	_open_Vertex_Action->setStatusTip( tr(" Open a vertex file") );
	connect( _open_Vertex_Action, SIGNAL(triggered()), this, SLOT(openTriangleGeometryName()));
	connect( this, SIGNAL(openTriangleGeometrySig(QString)), _mainGLView, SLOT(openTriangleGeometry(QString)));

	_open_Topology_Action = new QAction( tr("&OpenT"), this );
	_open_Topology_Action->setStatusTip( tr("Open a topology file"));
	connect( _open_Topology_Action, SIGNAL(triggered()), this, SLOT(openTriangleTopologyName()));
	connect( this, SIGNAL(openTriangleTopologySig(QString)), _mainGLView, SLOT(openTriangleTopology(QString)));

	_add_noise = new QAction( tr( "&Noise" ), this );
	_add_noise->setStatusTip( tr( "Add Gaussion noise on points" ) );
	connect( _add_noise, SIGNAL( triggered() ), _mainGLView, SLOT(PositionGaussianNoise()));
	//connect(_add_noise,SIGNAL(triggered()),this,SLOT(addGaussionNoise()));
  //   connect(_add_noise,SIGNAL(triggered()),this,SLOT(addGaussionNoise()));
	_pro_noise = new QAction( tr( "&Pro_Noise" ), this );
	_pro_noise->setStatusTip( tr( "Process Noise" ) );
	connect( _pro_noise, SIGNAL( triggered() ), _mainGLView, SLOT(gradientFieldFilterRBFilter()));

	//connect( _add_noise, SIGNAL( triggered() ), _mainGLView, SLOT(addGaussionNoise() ) );
//	connect(_add_noise,SIGNAL(triggered()),this,SLOT(addGaussionNoise()));
	




	


}

void GeoProc::createMenus(){
	fileMenu = menuBar()->addMenu( tr("&File"));
	fileMenu->addAction( _open_Vertex_Action );
	fileMenu->addAction( _open_Topology_Action );

	fileMenu = menuBar()->addMenu( tr("&GeoP") );
	fileMenu->addAction( _add_noise ); 
	fileMenu->addAction( _pro_noise ); 
}
void GeoProc::openTriangleGeometryName() {

   QString s = QFileDialog::getOpenFileName(this, "Choose a file name",
                                           *currOpenFileDir);
  if ( !(s.isEmpty() || s.isNull()) ){
    *currOpenFileDir = s.section('/',0,-2);
   emit(openTriangleGeometrySig(s));
  }


}



void GeoProc::openTriangleTopologyName()
{
  QString s = QFileDialog::getOpenFileName(this, "Choose a file name",
                                           *currOpenFileDir);
  if ( !(s.isEmpty() || s.isNull()) ){
    *currOpenFileDir = s.section('/',0,-2);
   emit(openTriangleTopologySig(s));
  }
}


void GeoProc::addGaussionNoise()
{
	
  //QGLWidget::update();
	//Position


	//QMessageBox message(QMessageBox::NoIcon,"Show Qt","Do you want to show Qt dialog?", QMessageBox::Yes | QMessageBox::No, NULL);  
	//if(message.exec() == QMessageBox::Yes)  
	//{  
	//	QMessageBox::aboutQt(NULL,"About Qt");  
	//}  
}
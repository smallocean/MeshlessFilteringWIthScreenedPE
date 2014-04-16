#include "geoproc.h"
#include <QtGui/QApplication>
#include <qmessagebox.h>
int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	GeoProc w;



	w.show();
	return a.exec();
}

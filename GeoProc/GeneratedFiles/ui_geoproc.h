/********************************************************************************
** Form generated from reading UI file 'geoproc.ui'
**
** Created: Fri Jun 1 14:39:07 2012
**      by: Qt User Interface Compiler version 4.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GEOPROC_H
#define UI_GEOPROC_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_GeoProcClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *GeoProcClass)
    {
        if (GeoProcClass->objectName().isEmpty())
            GeoProcClass->setObjectName(QString::fromUtf8("GeoProcClass"));
        GeoProcClass->resize(600, 400);
        menuBar = new QMenuBar(GeoProcClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        GeoProcClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(GeoProcClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        GeoProcClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(GeoProcClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        GeoProcClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(GeoProcClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        GeoProcClass->setStatusBar(statusBar);

        retranslateUi(GeoProcClass);

        QMetaObject::connectSlotsByName(GeoProcClass);
    } // setupUi

    void retranslateUi(QMainWindow *GeoProcClass)
    {
        GeoProcClass->setWindowTitle(QApplication::translate("GeoProcClass", "GeoProc", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class GeoProcClass: public Ui_GeoProcClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GEOPROC_H

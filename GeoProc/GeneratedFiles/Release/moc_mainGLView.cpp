/****************************************************************************
** Meta object code from reading C++ file 'mainGLView.h'
**
** Created: Mon Apr 8 13:57:03 2013
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../mainGLView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainGLView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_mainGLView[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      17,   11,   12,   11, 0x0a,
      50,   11,   12,   11, 0x0a,
      82,   11,   12,   11, 0x0a,
     112,   11,   12,   11, 0x0a,
     142,   11,   12,   11, 0x0a,
     171,   11,   12,   11, 0x0a,
     195,   11,   12,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_mainGLView[] = {
    "mainGLView\0\0bool\0openPointCloudsGeometry(QString)\0"
    "openPointCloudsTexture(QString)\0"
    "openTriangleGeometry(QString)\0"
    "openTriangleTopology(QString)\0"
    "openTriangleTexture(QString)\0"
    "PositionGaussianNoise()\0"
    "gradientFieldFilterRBFilter()\0"
};

void mainGLView::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        mainGLView *_t = static_cast<mainGLView *>(_o);
        switch (_id) {
        case 0: { bool _r = _t->openPointCloudsGeometry((*reinterpret_cast< QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 1: { bool _r = _t->openPointCloudsTexture((*reinterpret_cast< QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 2: { bool _r = _t->openTriangleGeometry((*reinterpret_cast< QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 3: { bool _r = _t->openTriangleTopology((*reinterpret_cast< QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 4: { bool _r = _t->openTriangleTexture((*reinterpret_cast< QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 5: { bool _r = _t->PositionGaussianNoise();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 6: { bool _r = _t->gradientFieldFilterRBFilter();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        default: ;
        }
    }
}

const QMetaObjectExtraData mainGLView::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject mainGLView::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_mainGLView,
      qt_meta_data_mainGLView, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &mainGLView::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *mainGLView::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *mainGLView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_mainGLView))
        return static_cast<void*>(const_cast< mainGLView*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int mainGLView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE


//////////////////////////////////////////////////////////////////////////



#if !defined (NorstructNormstruct)
#define NorstructNormstruct
#include <map>
#include <vector>
struct KnearestField{
	unsigned int m_numOfNearest;//���������������С��k
	bool m_IfBoundary;//�Ƿ��Ǳ߽� �����Ϊ true������Ϊfalse
	std::vector<int> m_nearest;//��������
};// ����������Ľṹ

struct Triangle{
	 int vertexIndex[3];
};


struct Color {
	float pointColor[3];
};


struct POINTVECTOR3D {
	float pointVector[3];
};

struct Point{
	float vertexPosition[3];
};
#endif
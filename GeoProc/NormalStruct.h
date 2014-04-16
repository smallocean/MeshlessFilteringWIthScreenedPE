
//////////////////////////////////////////////////////////////////////////



#if !defined (NorstructNormstruct)
#define NorstructNormstruct
#include <map>
#include <vector>
struct KnearestField{
	unsigned int m_numOfNearest;//邻域的数量，必须小于k
	bool m_IfBoundary;//是否是边界 如果是为 true，否则为false
	std::vector<int> m_nearest;//存放邻域点
};// 保存点的邻域的结构

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
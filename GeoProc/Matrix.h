#include <math.h>

//��4��4���������ʽ
template <class T>  T Matrix4Det(T a[4][4])
{
	return a[0][0]*a[1][1]*a[2][2]*a[3][3]-a[0][0]*a[1][1]*a[2][3]*a[3][2]-a[0][0]*a[1][2]*a[2][1]*a[3][3]
		+a[0][0]*a[1][2]*a[2][3]*a[3][1]+a[0][0]*a[1][3]*a[2][1]*a[3][2]-a[0][0]*a[1][3]*a[2][2]*a[3][1]
		-a[0][1]*a[1][0]*a[2][2]*a[3][3]+a[0][1]*a[1][0]*a[2][3]*a[3][2]+a[0][1]*a[1][2]*a[2][0]*a[3][3]
		-a[0][1]*a[1][2]*a[2][3]*a[3][0]-a[0][1]*a[1][3]*a[2][0]*a[3][2]+a[0][1]*a[1][3]*a[2][2]*a[3][0]
		+a[0][2]*a[1][0]*a[2][1]*a[3][3]-a[0][2]*a[1][0]*a[2][3]*a[3][1]-a[0][2]*a[1][1]*a[2][0]*a[3][3]
		+a[0][2]*a[1][1]*a[2][3]*a[3][0]+a[0][2]*a[1][3]*a[2][0]*a[3][1]-a[0][2]*a[1][3]*a[2][2]*a[3][0]
		-a[0][3]*a[1][0]*a[2][1]*a[3][2]+a[0][3]*a[1][0]*a[2][2]*a[3][1]+a[0][3]*a[1][1]*a[2][0]*a[3][2]
		-a[0][3]*a[1][1]*a[2][2]*a[3][0]+a[0][3]*a[1][2]*a[2][1]*a[4][0]-a[0][3]*a[1][2]*a[2][0]*a[3][1];
}

//��3��3���������ʽ
template <class T> T Matrix3Det(T a[3][3])
{
	return a[0][0]*a[1][1]*a[2][2]-a[0][0]*a[1][2]*a[2][1]+a[0][1]*a[1][2]*a[2][0]
		-a[0][1]*a[1][0]*a[2][2]+a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0];
}




//��2��2���������ʽ
template <class T> T Matrix2Det(T a[2][2])
{
	return a[0][0]*a[1][1]-a[0][1]*a[1][0];
}


//����3��3�������
template <class T> void Matrix3Multi(T a[3][3] , T b[3][3],T c[3][3])
{
   c[0][0]=a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0];
   c[0][1]=a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1];
   c[0][2]=a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2];

   c[1][0]=a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0];
   c[1][1]=a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1];
   c[1][2]=a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2];

   c[2][0]=a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0];
   c[2][1]=a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1];
   c[2][2]=a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2];

   return;
}
//����2*2�������
template <class T> void Matrix2Multi(T a[2][2],T b[2][2],T c[2][2])
{
	c[0][0]=a[0][0]*b[0][0]+a[0][1]*b[1][0];
	c[0][1]=a[0][0]*b[0][1]+a[0][1]*b[1][1];
	c[1][0]=a[1][0]*b[0][0]+a[1][1]*b[0][1];
	c[1][1]=a[1][0]*b[0][1]+a[1][1]*b[1][1];
	return;
}




//��2��2���������
template <class T>  void Matrix2Invert(T a[2][2],T b[2][2])
{
	T c;
	c=Matrix2Det(a);

	b[0][0]=a[1][1]/c;
	b[0][1]=-a[0][1]/c;
	b[1][0]=-a[1][0]/c;
	b[1][1]=a[0][0]/c;

    return;
}

//��һ����λ�����Ĵ�ֱ����
template <class T> void  Perpendicular(T a[3],T b[3])
{
	
	if (fabs(a[0])<fabs(a[1]))
	{
		if (fabs(a[0])<fabs(a[2]))
		{
			b[0]=1-a[0]*a[0];
			b[1]=-a[0]*a[1];
			b[2]=-a[0]*a[2];
			float sum;
			sum=0;
			for(int i=0;i<3;i++)
			{
				sum+=b[i]*b[i];
			}
			sum=sqrt(sum);
			for(int i=0;i<3;i++)
			{
				b[i]=b[i]/sum;
			}
			return;
			
		}
	}
	else{
		if (fabs(a[1])<fabs(a[2]))
		{
			b[0]=-a[0]*a[1];
			b[1]=1-a[1]*a[1];
			b[2]=-a[1]*a[2];
			float sum;
			sum=0;
			for(int i=0;i<3;i++)
			{
				sum+=b[i]*b[i];
			}
			sum=sqrt(sum);
			for(int i=0;i<3;i++)
			{
				b[i]=b[i]/sum;
			}
			return;
			
		}
	}
	   
		
		
            b[0]=-a[0]*a[2];
            b[1]=-a[1]*a[2];
            b[2]=1-a[2]*a[2];

			float sum;
			sum=0;
			for(int i=0;i<3;i++)
			{
				sum+=b[i]*b[i];
			}
			sum=sqrt(sum);
			for(int i=0;i<3;i++)
			{
				b[i]=b[i]/sum;
			}
			return;
	
}

//����1��3�����Ĳ��
template <class T> void Cross3(T a[3], T b[3],T c[3])
{
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
	return;
}


//4*4�����������ĳ˷�
template <class T> void Matrix4vector(T a[4][4], T b[4],T c[4])
{
	
	for(int i=0;i<4;i++)
		c[i]=a[i][0]*b[0]+a[i][1]*b[1]+a[i][2]*b[2]+a[i][3]*b[3];
	return;
}

// ���������ĵ��
template <class T> T Dot4(T a[4],T b[4])
{
	T c;
	c=0;
	for (int i=0;i<4;i++)
	{
		c+=a[i]*b[i];
	}
	return c;
}



//3*3���������
template <class T> void Matrix3Invert(T a[3][3],T b[3][3])
{
	T c;
	c=Matrix3Det(a);

	b[0][0]=a[1][1]*a[2][2]-a[1][2]*a[2][1];
	b[0][1]=-(a[0][1]*a[2][2]-a[0][2]*a[2][1]);
	b[0][2]=a[0][1]*a[1][2]-a[0][2]*a[1][1];
	b[1][0]=-(a[1][0]*a[2][2]-a[1][2]*a[2][0]);
	b[1][1]=a[0][0]*a[2][2]-a[0][2]*a[2][0];
	b[1][2]=-(a[0][0]*a[1][2]-a[0][2]*a[1][0]);
	b[2][0]=a[1][0]*a[2][1]-a[1][1]*a[2][0];
	b[2][1]=-(a[0][0]*a[2][1]-a[0][1]*a[2][0]);
	b[2][2]=a[0][0]*a[1][1]-a[0][1]*a[1][0];
	
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
       {
           b[i][j]=b[i][j]/c;
	}
	return;
}


//3��3�����ת��
template <class T> void Matrix3Transpose(T a[3][3], T b[3][3])
{
	b[0][0]=a[0][0];
	b[0][1]=a[1][0];
	b[0][2]=a[2][0];
	b[1][0]=a[0][1];
	b[1][1]=a[1][1];
	b[1][2]=a[2][1];
	b[2][0]=a[0][2];
	b[2][1]=a[1][2];
	b[2][2]=a[2][2];
	
	return;
}
//1*3�����ĵ��
template <class T> T Vector3Vector(T a[3],T b[3])
{
	float sum;
	sum=0;
	for(int i=0;i<3;i++)
		sum+=a[i]*b[i];
	return sum;
}
// 3*1������1*3�����ĳ�
template <class T> void Vector3TVector(T a[3],T b[3],T c[3][3])
{

	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			c[i][j]=a[i]*b[j];
		}
	}
	return;
}

// 3*n������ n*3 �������
template <class T> void MatrixN3TMatrix(int n,T* a,T* b,T c[3][3])
{
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			c[i][j]=0;
			for(int k=0;k<n;k++){
				c[i][j]+=a[k*3+i]*b[k*3+j];
			}
		}
	}
	return;
}
// ��һ��1*3�����ľ���
template <class T> T ComputeNormOfVector3(T a[3])
{
	T sum=0;
	for(int i=0;i<3;i++){
		sum+=a[i]*a[i];
	}
	sum=sqrt(sum);
	return sum;	
}
// 3*3�����3*3����
template <class T> void Matrix3MinusMatrix3(T a[3][3],T b[3][3],T c[3][3])
{
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			c[i][j]=a[i][j]-b[i][j];
	return;
}

//3*3�����������ĳ˷�
template <class T> void Matrix3MultiVector(T a[3][3],T b[3],T c[3])
{
	for(int i=0;i<3;i++)
		c[i]=a[i][0]*b[0]+a[i][1]*b[1]+a[i][2]*b[2];
	return;
}

// ����3*3����ĳ˷�
template <class T> void Matrix3MultiScalar(T a[3][3],T b,T c[3][3])
{
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			c[i][j]=b*a[i][j];
		}
	}
	return;
}
/*
// 1*3������3*3����ĳ˷�
template <class T> void Vector3Matrix(T a[3], T b[3][3], T c[3])
{
*/
//////////////////////////////////////////////////////////////////////////
//��λ������
template <class T> void UnitVector(T a[3],T b[3])
{
	float sum;
	sum=0;
	for(int i=0;i<3;i++)
	{
		sum+=a[i]*a[i];
	}
	sum=sqrt(sum);
	for(int i=0;i<3;i++)
	{
		b[i]=a[i]/sum;
	}
	return;

}



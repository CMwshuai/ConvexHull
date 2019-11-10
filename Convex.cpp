#include "Convex.h"

#ifndef PI
#define PI			3.14159
#endif // !PI

#define MIN(x,y)	(((x) < (y)) ? (x):(y))
#define MAX(x,y)	(((x) > (y)) ? (x):(y))

#define DEGTORAD(deg)	(((deg)*2.0*PI)/360.0)
#define RADTODEG(rad)	(((rad)*360.0)/(2.0*PI))

#define SEGMENTLEN(x0,y0,x1,y1) (sqrt(pow(((x1)-(x0)), 2.0) + pow(((y1)-(y0)), 2.0)))

typedef struct Point_ {
	double x;
	double y;
	double z;
} Point;

typedef struct SPoint_ {
	double rho;
	double theta;
	double phi;
} SPoint;


Convex::Convex(QObject *parent)
	: QObject(parent)
{
}

Convex::Convex()
{
}

Convex::~Convex()
{
}

qreal Convex::comparePointClock(const QPointF &point_0, const QPointF &point_c, const QPointF &point_i)
{
	return ((point_i.x() - point_0.x())*(point_c.y() - point_0.y()) - (point_i.y() - point_0.y())*(point_c.x() - point_0.x()));
}

quint32 Convex::removeRepeatPoints(QVector<QPointF> &vecPoints)
{
	if (vecPoints.isEmpty())
		return 0;

	QVector<QPointF> tempVecPorint;
	tempVecPorint = vecPoints;
	vecPoints.clear();
	QPointF tempPoint;
	while (tempVecPorint.size())
	{
		tempPoint = tempVecPorint.at(0);
		tempVecPorint.removeAll(tempPoint);
		vecPoints.push_back(tempPoint);
	}

	return vecPoints.size();
}

QPointF Convex::getMinimumPoint(const QVector<QPointF> &vecPoints)
{
	if (vecPoints.isEmpty())
		return QPointF();

	QPointF minPoint = vecPoints.at(0);
	quint16 point_x = vecPoints.at(0).x(), point_y = vecPoints.at(0).y();
	for (QVector<QPointF>::const_iterator it = vecPoints.constBegin(); it != vecPoints.constEnd(); it++)
	{
		//比较Y坐标，找Y坐标最小的
		if (it->y() < minPoint.y())
		{
			minPoint = (*it);
		}
		else
		{
			//Y坐标相同，找X坐标小的
			if (it->y() == minPoint.y() && it->x() < minPoint.x())
			{
				minPoint = (*it);
			}
		}
	}

	return minPoint;
}

//Jarvis's march 算法，O(nH)，H为点的个数。
qint8 Convex::getConvexHullJarvis(const QVector<QPointF> &vecSourPoints, QVector<QPointF> &vecTarPoints)
{
	if (vecSourPoints.isEmpty())
		return -1;

	QPointF minPoint;
	QPointF lowPoint, point_0, point_i, point_c;
	qreal count = 0,z = 0;
	qreal length_1, length_2;
	QVector<QPointF> tempVecPoint(vecSourPoints);

	vecTarPoints.clear();
	//删除重复坐标
	if (removeRepeatPoints(tempVecPoint) <= 0)
		return -1;

	//查找最小坐标
	minPoint = getMinimumPoint(tempVecPoint);
	lowPoint = minPoint;
	
	point_0 = lowPoint;

	do {
		//起始点point_0压入凸包点集中
		vecTarPoints.push_back(point_0);
		
		count = 0;
		for (QVector<QPointF>::iterator it = tempVecPoint.begin(); it != tempVecPoint.end(); it++)
		{
			//跳过起始坐标
			if ((*it) == point_0)
				continue;

			count++;
			if (count == 1) //把第一个便利的点作为point_c
			{
				point_c = (*it);
				continue;
			}
			//如果z>0则point在point_i和point_c连线的下方，z<0则point_i在连线的上方，z=0则point_i共线
			z = comparePointClock(point_0,point_c,(*it));//((it->x() - point_0.x())*(point_c.y() - point_0.y()) - (it->y() - point_0.y())*(point_c.x() - point_0.x()));
			if (z > 0)
			{
				point_c = (*it);
			}
			else if (z == 0)
			{
				//共线情况找出距离point_0较远的那个点作为point_c
				length_1 = SEGMENTLEN(point_0.x(),point_0.y(),it->x(),it->y());
				length_2 = SEGMENTLEN(point_0.x(), point_0.y(), point_c.x(), point_c.y());
				if (length_1 > length_2)
				{
					point_c = (*it);
				}
			}
		}
		point_0 = point_c;

	} while (point_0 != lowPoint);
	vecTarPoints.push_back(lowPoint);
	if (vecTarPoints.isEmpty())
		return -1;
	return 0;
}


QPointF m_point0;
bool comPolarAngle(const QPointF &point_1, const QPointF &point_2)
{
	qreal z = ((point_2.x() - m_point0.x())*(point_1.y() - m_point0.y()) - (point_2.y() - m_point0.y())*(point_1.x() - m_point0.x()));
	if (fabs(z) < 1e-6)
	{
		qreal length_1 = SEGMENTLEN(m_point0.x(), m_point0.y(), point_1.x(), point_1.y());
		qreal length_2 = SEGMENTLEN(m_point0.x(), m_point0.y(), point_2.x(), point_2.y());

		return length_1 > length_2;
	}
	else
	{
		return z < 0;
	}
}

bool Convex::sortByPolarAngle(QVector<QPointF> &vecPoints)
{
	if (vecPoints.isEmpty())
		return false;

	QVector<QPointF> tempVecPoint(vecPoints);
	tempVecPoint.removeOne(m_point0);

	qreal z = 0;
	qSort(tempVecPoint.begin(), tempVecPoint.end(), comPolarAngle);
	tempVecPoint.push_front(m_point0);

	vecPoints = tempVecPoint;

	return true;
}
//Graham 扫描算法，O(nlgn)。
qint8 Convex::getConvecHullGraham(const QVector<QPointF> &vecSourPoints, QVector<QPointF> &vecTarPoints)
{
	if (vecSourPoints.isEmpty())
		return -1;
	
	QVector<QPointF> tempVecPoint(vecSourPoints);

	//删除重复坐标
	if (removeRepeatPoints(tempVecPoint) <= 0)
		return -1;

	//查找最小坐标
	QPointF minPoint;
	minPoint = getMinimumPoint(tempVecPoint);
	m_point0 = minPoint;

	//按极角进行排序
	if(!sortByPolarAngle(tempVecPoint))
		return -1;

	vecTarPoints.clear();
	vecTarPoints.push_back(tempVecPoint.at(0));
	vecTarPoints.push_back(tempVecPoint.at(1));
	vecTarPoints.push_back(tempVecPoint.at(2));
	qint32 vecTop = 2;
	for (int i = 3; i < tempVecPoint.size(); i++)
	{
		while (vecTop > 0
			&& (comparePointClock(vecTarPoints.at(vecTop - 1), vecTarPoints.at(vecTop), tempVecPoint.at(i)) >= 0))
		{
			vecTop--;
			vecTarPoints.pop_back();
		}
		vecTarPoints.push_back(tempVecPoint.at(i));
		vecTop++;
	}
	vecTarPoints.push_back(minPoint);

	if (vecTarPoints.isEmpty())
		return -1;
	return 0;
}
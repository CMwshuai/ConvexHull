#ifndef _CONVEX_H_
#define _CONVEX_H_

#include <QObject>
#include <QVector>
class Convex : public QObject
{
	Q_OBJECT

public:
	Convex(QObject *parent);
	Convex();
	~Convex();

	//Jarvis's march 步进算法，O(nH)，H为点的个数。
	qint8 getConvexHullJarvis(const QVector<QPointF> &vecSourPoints, QVector<QPointF> &vecTarPoints);
	//Graham 扫描算法，O(nlgn)。
	qint8 getConvecHullGraham(const QVector<QPointF> &vecSourPoints, QVector<QPointF> &vecTarPoints);

private:
	qreal comparePointClock(const QPointF &point_0,const QPointF &point_c,const QPointF &point_i);
	quint32 removeRepeatPoints(QVector<QPointF> &vecPoints);
	QPointF getMinimumPoint(const QVector<QPointF> &vecPoints);
	bool	sortByPolarAngle(QVector<QPointF> &vecPoints);

};

#endif
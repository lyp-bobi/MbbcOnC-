
#pragma once

namespace SpatialIndex
{
	class SIDX_DLL STPoint : public Point, public ITimeShape
	{
	public:
		STPoint();
		STPoint(const double* pCoords, const Tools::IInterval& ti, uint32_t dimension);
		STPoint(const double* pCoords, double t,uint32_t dimension);
		STPoint(const Point& p, const Tools::IInterval& ti);
		STPoint(const Point& p, double tStart, double tEnd);
		STPoint(const STPoint& p);
		virtual ~STPoint();

		virtual STPoint& operator=(const STPoint& p);
		virtual bool operator==(const STPoint& p) const;

		//
		// IObject interface
		//
		virtual STPoint* clone();

		//
		// ISerializable interface
		//
		virtual uint32_t getByteArraySize() const;
		virtual void loadFromByteArray(const uint8_t* data);
		virtual void storeToByteArray(uint8_t** data, uint32_t& len);

		//
		// ITimeShape interface
		//
		virtual bool intersectsShapeInTime(const ITimeShape& in) const;
		virtual bool intersectsShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const;
		virtual bool containsShapeInTime(const ITimeShape& in) const;
		virtual bool containsShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const;
		virtual bool touchesShapeInTime(const ITimeShape& in) const;
		virtual bool touchesShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const;
		virtual double getAreaInTime() const;
		virtual double getAreaInTime(const Tools::IInterval& ivI) const;
		virtual double getIntersectingAreaInTime(const ITimeShape& r) const;
		virtual double getIntersectingAreaInTime(const Tools::IInterval& ivI, const ITimeShape& r) const;

		//
		// IInterval interface
		//
		virtual Tools::IInterval& operator=(const Tools::IInterval&);
		virtual double getLowerBound() const;
		virtual double getUpperBound() const;
		virtual void setBounds(double, double);
		virtual bool intersectsInterval(const Tools::IInterval& ti) const;
		virtual bool intersectsInterval(Tools::IntervalType t, const double start, const double end) const;
		virtual bool containsInterval(const Tools::IInterval& ti) const;
		virtual Tools::IntervalType getIntervalType() const;

		virtual void makeInfinite(uint32_t dimension);
		virtual void makeDimension(uint32_t dimension);

	public:
		double m_time;
        static STPoint* makemid(const STPoint &tp1,const STPoint &tp2,double t);
		friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const STPoint& pt);
	}; // STPoint

	SIDX_DLL std::ostream& operator<<(std::ostream& os, const STPoint& pt);
}

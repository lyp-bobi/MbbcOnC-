
#include <cstring>
#include <limits>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

STPoint::STPoint()
	: Point(), m_time(-std::numeric_limits<double>::max())
{
}

STPoint::STPoint(const double* pCoords, const IInterval& ti, uint32_t dimension)
	: Point(pCoords, dimension), m_time(ti.getLowerBound())
{
}

STPoint::STPoint(const double* pCoords, double tStart, uint32_t dimension)
	: Point(pCoords, dimension), m_time(tStart)
{
}

STPoint::STPoint(double x, double y, double t) {
    makeInfinite(2);
    m_pCoords[0]=x;
    m_pCoords[1]=y;
    m_time=t;
}

STPoint::STPoint(const Point& p, const IInterval& ti)
	: Point(p), m_time(ti.getLowerBound())
{
}

STPoint::STPoint(const Point& p, double tStart)
	: Point(p), m_time(tStart)
{
}

STPoint::STPoint(const STPoint& p)
	: m_time(p.m_time)
{
	m_dimension = p.m_dimension;

	m_pCoords = new double[m_dimension];
	memcpy(m_pCoords, p.m_pCoords, m_dimension * sizeof(double));
}

STPoint::~STPoint()
{
}

STPoint& STPoint::operator=(const STPoint& p)
{
	if (this != &p)
	{
		makeDimension(p.m_dimension);
		memcpy(m_pCoords, p.m_pCoords, m_dimension * sizeof(double));
		m_time = p.m_time;
		m_time = p.m_time;
	}

	return *this;
}

bool STPoint::operator==(const STPoint& p) const
{
	if (
		m_time < p.m_time - std::numeric_limits<double>::epsilon() ||
		m_time > p.m_time + std::numeric_limits<double>::epsilon() ||
		m_time < p.m_time - std::numeric_limits<double>::epsilon() ||
		m_time > p.m_time + std::numeric_limits<double>::epsilon())
		return false;

	for (uint32_t cDim = 0; cDim < m_dimension; ++cDim)
	{
		if (
			m_pCoords[cDim] < p.m_pCoords[cDim] - std::numeric_limits<double>::epsilon() ||
			m_pCoords[cDim] > p.m_pCoords[cDim] + std::numeric_limits<double>::epsilon()) 
			return false;
	}

	return true;
}

//
// IObject interface
//
STPoint* STPoint::clone()
{
	return new STPoint(*this);
}

//
// ISerializable interface
//
uint32_t STPoint::getByteArraySize() const
{
	return (sizeof(uint32_t) + sizeof(double) + m_dimension * sizeof(double));
}

void STPoint::loadFromByteArray(const uint8_t* ptr)
{
	uint32_t dimension;
	memcpy(&dimension, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(&m_time, ptr, sizeof(double));
	ptr += sizeof(double);

	makeDimension(dimension);
	memcpy(m_pCoords, ptr, m_dimension * sizeof(double));
	//ptr += m_dimension * sizeof(double);
}

void STPoint::storeToByteArray(uint8_t** data, uint32_t& len)
{
	len = getByteArraySize();
	*data = new uint8_t[len];
	uint8_t* ptr = *data;

	memcpy(ptr, &m_dimension, sizeof(uint32_t));
	ptr += sizeof(uint32_t);
	memcpy(ptr, &m_time, sizeof(double));
	ptr += sizeof(double);
	memcpy(ptr, m_pCoords, m_dimension * sizeof(double));
	//ptr += m_dimension * sizeof(double);
}

//
// ITimeShape interface
//
bool STPoint::intersectsShapeInTime(const ITimeShape& in) const
{
//	const TimeRegion* pr = dynamic_cast<const TimeRegion*>(&in);
//	if (pr != 0) return pr->containsPointInTime(*this);

	throw Tools::IllegalStateException("intersectsShapeInTime: Not implemented yet!");
}

bool STPoint::intersectsShapeInTime(const IInterval&, const ITimeShape&) const
{
	throw Tools::IllegalStateException("intersectsShapeInTime: Not implemented yet!");
}

bool STPoint::containsShapeInTime(const ITimeShape&) const
{
	return false;
}

bool STPoint::containsShapeInTime(const IInterval&, const ITimeShape&) const
{
	return false;
}

bool STPoint::touchesShapeInTime(const ITimeShape&) const
{
	throw Tools::IllegalStateException("touchesShapeInTime: Not implemented yet!");
}

bool STPoint::touchesShapeInTime(const IInterval&, const ITimeShape&) const
{
	throw Tools::IllegalStateException("touchesShapeInTime: Not implemented yet!");
}

double STPoint::getAreaInTime() const
{
	return 0.0;
}

double STPoint::getAreaInTime(const IInterval&) const
{
	return 0.0;
}

double STPoint::getIntersectingAreaInTime(const ITimeShape&) const
{
	return 0.0;
}

double STPoint::getIntersectingAreaInTime(const IInterval&, const ITimeShape&) const
{
	return 0.0;
}

//
// IInterval interface
//
Tools::IInterval& STPoint::operator=(const Tools::IInterval& i)
{
	if (this != &i)
	{
		m_time = i.getLowerBound();
		m_time = i.getUpperBound();
	}

	return *this;
}

double STPoint::getLowerBound() const
{
	return m_time;
}

double STPoint::getUpperBound() const
{
	return m_time;
}

void STPoint::setBounds(double l, double h)
{
	assert(l <= h);

	m_time = l;
	m_time = h;
}

bool STPoint::intersectsInterval(const IInterval& ti) const
{
	return intersectsInterval(ti.getIntervalType(), ti.getLowerBound(), ti.getUpperBound());
}

bool STPoint::intersectsInterval(Tools::IntervalType, const double start, const double end) const
{
	//if (m_time != start &&
	//		(m_time >= end || m_time <= start)) return false;
	if (m_time >= end || m_time <= start) return false;

	return true;
}

bool STPoint::containsInterval(const IInterval& ti) const
{
	if (m_time <= ti.getLowerBound() && m_time >= ti.getUpperBound()) return true;
	return false;
}

Tools::IntervalType STPoint::getIntervalType() const
{
	return Tools::IT_RIGHTOPEN;
}

void STPoint::makeInfinite(uint32_t dimension)
{
	makeDimension(dimension);
	for (uint32_t cIndex = 0; cIndex < m_dimension; ++cIndex)
	{
		m_pCoords[cIndex] = std::numeric_limits<double>::max();
	}

	m_time = std::numeric_limits<double>::max();
}

void STPoint::makeDimension(uint32_t dimension)
{
	if (m_dimension != dimension)
	{
		m_dimension = dimension;

		delete[] m_pCoords;
		m_pCoords = 0;

		m_pCoords = new double[m_dimension];
	}
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const STPoint& pt)
{
	uint32_t i;

	for (i = 0; i < pt.m_dimension; ++i)
	{
		os << pt.m_pCoords[i] << " ";
	}

	os << ", " << pt.m_time;

	return os;
}

STPoint STPoint::makemid(const STPoint &p1, const STPoint &p2, double t){
    assert(p1.m_dimension==p2.m_dimension);
    if(p1.m_time==p2.m_time)
        return STPoint(p1);
    int dim=p1.m_dimension;
    double* p1c=new double[dim];
    double* p2c=new double[dim];
    double t1=p1.m_time;
    double t2=p2.m_time;
    for(int i=0;i<dim;i++){
        p1c[i]=p1.m_pCoords[i];
        p2c[i]=p2.m_pCoords[i];
    }
    double h1= (t-t1)/(t2-t1);
    double h2= (t2-t)/(t2-t1);
    double* p3c=new double[dim];
    for(int i=0;i<dim;i++){
        p3c[i]=h2*p1c[i]+h1*p2c[i];
    }
    delete[](p1c);
    delete[](p2c);
    auto res=STPoint(p3c,t,dim);
    delete[](p3c);
    return res;
}

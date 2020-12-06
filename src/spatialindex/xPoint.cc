
#include <cstring>
#include <limits>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

xPoint::xPoint(prex x, prex y, prex t) {
    makeInfinite(2);
    m_x=x;
    m_y=y;
    m_t=t;
}


xPoint::xPoint(const xPoint& p)
	:m_x(p.m_x),m_y(p.m_y), m_t(p.m_t)
{
}

xPoint::~xPoint()
{
}

xPoint& xPoint::operator=(const xPoint& p)
{
	if (this != &p)
	{
		m_x = p.m_x;
		m_y = p.m_y;
		m_t = p.m_t;
	}

	return *this;
}

bool xPoint::operator==(const xPoint& p) const
{
	if (
		m_t < p.m_t - std::numeric_limits<double>::epsilon() ||
		m_t > p.m_t + std::numeric_limits<double>::epsilon() ||
        m_x < p.m_x - std::numeric_limits<double>::epsilon() ||
        m_x > p.m_x + std::numeric_limits<double>::epsilon()||
        m_y < p.m_y - std::numeric_limits<double>::epsilon() ||
        m_y > p.m_y + std::numeric_limits<double>::epsilon())
		return false;

	return true;
}

xPoint xPoint::operator+(const xPoint &b) const {
    xPoint res;
    res.m_x = m_x + b.m_x;
    res.m_y = m_y + b.m_y;
    res.m_t = m_t + b.m_t;
    return res;
}

xPoint xPoint::operator-(const xPoint &b) const {
    xPoint res;
    res.m_x = m_x - b.m_x;
    res.m_y = m_y - b.m_y;
    res.m_t = m_t - b.m_t;
    return res;
}

//
// IObject interface
//
xPoint* xPoint::clone()
{
	return new xPoint(*this);
}

//
// ISerializable interface
//
uint32_t xPoint::getByteArraySize() const
{
	return 3* sizeof(prex);
}

void xPoint::loadFromByteArray(const uint8_t* ptr)
{
    memcpy(&m_x, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_y, ptr, sizeof(prex));
    ptr += sizeof(prex);
	memcpy(&m_t, ptr, sizeof(prex));
//	ptr += sizeof(prex);
}

void xPoint::storeToByteArray(uint8_t** data, uint32_t& len)
{
	len = getByteArraySize();
	*data = new uint8_t[len];
	uint8_t* ptr = *data;

    memcpy(ptr, &m_x, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_y, sizeof(prex));
    ptr += sizeof(prex);
	memcpy(ptr, &m_t, sizeof(prex));
//	ptr += sizeof(prex);

}

void xPoint::storeToByteArrayE(uint8_t** data, uint32_t& len)
{
    len = getByteArraySize();
    uint8_t* ptr = *data;

    memcpy(ptr, &m_x, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_y, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_t, sizeof(prex));
//	ptr += sizeof(prex);

}

void xPoint::makeInfinite(uint32_t dimension)
{
}

void xPoint::makeDimension(uint32_t dimension)
{
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const xPoint& pt)
{
	uint32_t i;

	os<<pt.m_x<<","<<pt.m_y << "," << pt.m_t;

	return os;
}

xPoint xPoint::makemid(const xPoint &p1, const xPoint &p2, double t){
    if(p1.m_t==p2.m_t)
        return xPoint(p1);
    double p1c[2];
    double p2c[2];
    double t1=p1.m_t;
    double t2=p2.m_t;
    p1c[0]=p1.m_x;
    p2c[0]=p2.m_x;
    p1c[1]=p1.m_y;
    p2c[1]=p2.m_y;
    double h1= (t-t1)/(t2-t1);
    double h2= (t2-t)/(t2-t1);
    double p3c[2];
    for(int i=0;i<2;i++){
        p3c[i]=h2*p1c[i]+h1*p2c[i];
    }
    auto res=xPoint(p3c[0],p3c[1],t);
    return res;
}


//
// IShape interface
//
bool xPoint::intersectsShape(const IShape& s) const
{

    throw Tools::IllegalStateException(
            "Point::intersectsShape: Not implemented yet!"
    );
}

bool xPoint::containsShape(const IShape&) const
{
    return false;
}

bool xPoint::touchesShape(const IShape& s) const
{
    throw Tools::IllegalStateException(
            "Point::touchesShape: Not implemented yet!"
    );
}

void xPoint::getCenter(Point& out) const
{

}

uint32_t xPoint::getDimension() const
{
    return 3;
}

void xPoint::getxMBR(xMBR& out) const
{
    out.m_xmin=m_x;
    out.m_xmax=m_x;
    out.m_ymin=m_y;
    out.m_ymax=m_y;
    out.m_tmin=m_t;
    out.m_tmax=m_t;
}

double xPoint::getArea() const
{
    return 0.0;
}

double xPoint::getMinimumDistance(const IShape& s) const
{
    const xPoint* ppt = dynamic_cast<const xPoint*>(&s);
    if (ppt != 0)
    {
        return getMinimumDistance(*ppt);
    }

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != 0)
    {
        return pr->getMinimumDistance(*this);
    }

    throw Tools::IllegalStateException(
            "Point::getMinimumDistance: Not implemented yet!"
    );
}

double xPoint::getMinimumDistance(const xPoint& p) const
{

    double ret = 0.0;
    ret += std::pow(m_x - p.m_x, 2.0);
    ret += std::pow(m_y - p.m_y, 2.0);
    return std::sqrt(ret);
}

xPoint xPoint::rotate(xPoint &center, double angle) {
    xPoint delta = *this - center;
    xPoint d2;
    d2.m_x = delta.m_x * cos(angle) - delta.m_y * sin(angle);
    d2.m_y = delta.m_x * sin(angle) + delta.m_y * cos(angle);
    return center + d2;
}
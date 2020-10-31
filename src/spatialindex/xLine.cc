//
// Created by Chuang on 2019/5/29.
//
#define _USE_MATH_DEFINES

#include <cstring>
#include <math.h>
#include <limits>
#include <algorithm>
#include <cfloat>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
xLine::xLine()
{
}
xLine::xLine(const SpatialIndex::xLine &in)
:m_ps(in.m_ps),m_pe(in.m_pe){
}
xLine::xLine(xPoint ps, xPoint pe, double rd, double rv)
:m_ps(ps),m_pe(pe){}


xLine& xLine::operator=(const xLine& r)
{
    m_ps = r.m_ps;
    m_pe = r.m_pe;
    return *this;
}

bool xLine::operator==(const SpatialIndex::xLine &r) const {
    throw Tools::IllegalArgumentException(
            "MBC:== not supported"
    );
    return true;
}
//
// IObject interface
//
xLine* xLine::clone() {
    return new xLine(*this);
}
//
// ISerializable interface
//
uint32_t xLine::getByteArraySize() const {
    if(bCompactMBC) return 7*sizeof(prex);
    else return 8*sizeof(prex);
}

void xLine::loadFromByteArray(const uint8_t* ptr) {
    m_ps.loadFromByteArray(ptr);
    ptr+=m_ps.getByteArraySize();
    m_pe.loadFromByteArray(ptr);
//    ptr+=m_pe.getByteArraySize();
}

void xLine::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint32_t tmp;
    m_ps.storeToByteArrayE(&ptr,tmp);
    ptr+=tmp;
    m_pe.storeToByteArrayE(&ptr,tmp);
    len = getByteArraySize();
//    ptr+=tmp;
}

void xLine::storeToByteArrayE(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    uint8_t* ptr = *data;
    uint32_t tmp;
    m_ps.storeToByteArrayE(&ptr,tmp);
    ptr+=tmp;
    m_pe.storeToByteArrayE(&ptr,tmp);
    len = getByteArraySize();
//    ptr+=tmp;
}

std::pair<xPoint,double> xLine::getCenterRdAtTime(double t) const {
    if(t<=m_ps.m_t) return std::make_pair(xPoint(m_ps),0);
    if(t>=m_pe.m_t) return std::make_pair(xPoint(m_pe),0);
    double x=makemidmacro(m_ps.m_x,m_ps.m_t,m_pe.m_x,m_pe.m_t,t);
    double y=makemidmacro(m_ps.m_y,m_ps.m_t,m_pe.m_y,m_pe.m_t,t);
    double coord[2]={x,y};
//    xPoint *tp = xPoint::makemid(xPoint(m_pLow, m_ps.m_t, m_ps.m_t, m_dimension),
//                                      xPoint(m_pHigh, m_pe.m_t, m_pe.m_t, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double r=0;
    xPoint tp(coord[0],coord[1],t);
    return std::make_pair(tp,r);
}

//
// IShape interface
//
bool xLine::intersectsShape(const SpatialIndex::IShape& s) const {

    const xPoint* ptp = dynamic_cast<const xPoint*>(&s);
    if (ptp != 0) return intersectsxPoint(*ptp);

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != 0) return intersectsxMBR(*pr);

    const xCylinder* pcy = dynamic_cast<const xCylinder*>(&s);
    if (pcy != 0) return pcy->checkRel(*this);

//    const xTrajectory* ptra = dynamic_cast<const xTrajectory*>(&s);
//    if (ptra != 0) return ptra->intersectsShape(*this);

//    const xLine* pxLine = dynamic_cast<const xLine*>(&s);
//    if (pxLine != 0) return intersectsxLine(*pxLine);

//    const TimexMBR* ptr = dynamic_cast<const TimexMBR*>(&s);
//    if (ptr != 0) return intersectsTimexMBR(*ptr);
    throw Tools::NotSupportedException("xLine::intersectsShape:Not supported");
}

bool xLine::intersectsxPoint(const SpatialIndex::xPoint &in) const {
    if(in.m_t<m_ps.m_t||in.m_t>m_pe.m_t) return false;
    auto timed=getCenterRdAtTime(in.m_t);
    return timed.first.getMinimumDistance(in)<timed.second;
}

bool xLine::intersectsxMBR(const SpatialIndex::xMBR &in) const {
    //2 dimensions for space, the third for time
    if(in.m_tmin>m_pe.m_t||in.m_tmax<m_ps.m_t) return false;
    double t0 = m_ps.m_t, t3 = m_pe.m_t;
    double tlow = in.m_tmin, thigh = in.m_tmax;
    double ints = std::max(t0, tlow), inte = std::min(t3, thigh);
    auto a = getCenterRdAtTime(ints), b = getCenterRdAtTime(inte);
    double d = xTrajectory::line2MBRMinSED(a.first, b.first, in);
    if (d > 0) return false;
    return true;
}

inline bool xLine::intersectsxLine(const xLine& in) const{throw Tools::NotSupportedException("xLine::intersectsxLine");}

static std::pair<double,double> getIntersectPeriod(
        double xs,double xe,double ts,double te,
        double xlow,double xhigh){
    //assume the time dimension of xlow and xhigh is infinity
    double dt1=ts+(xlow-xs)/(xe-xs)*(te-ts);
    double dt2=te-(xe-xhigh)/(xe-xs)*(te-ts);
    double is=(dt1>xs&&dt1<xe)?dt1:xs,ie=(dt2>xs&&dt2<xe)?dt2:xe;
    if(is<ie)
        return std::make_pair(is,ie);
    else
        return std::make_pair(ie,is);
}


bool xLine::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool xLine::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("xLine:touchesShape");
}
void xLine::getCenter(Point& out) const{
}
uint32_t xLine::getDimension() const{return 3;}
void xLine::getxMBR(xMBR& out) const{
    out.m_xmin=min(m_ps.m_x, m_pe.m_x);
    out.m_xmax=max(m_ps.m_x, m_pe.m_x);
    out.m_ymin=min(m_ps.m_y, m_pe.m_y);
    out.m_ymax=max(m_ps.m_y, m_pe.m_y);
    out.m_tmin=min(m_ps.m_t, m_pe.m_t);
    out.m_tmax=max(m_ps.m_t, m_pe.m_t);
}

double xLine::getArea() const{
    return 0;
}
double xLine::getMinimumDistance(const IShape& in) const{
    const xPoint* ptp = dynamic_cast<const xPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const xMBR* pr = dynamic_cast<const xMBR*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "xLine::getMinimumDistance: Not implemented yet!"
    );
}

double xLine::getMinimumDistance(const SpatialIndex::xMBR &in) const {
    throw Tools::NotSupportedException("");
}
double xLine::getMinimumDistance(const SpatialIndex::xPoint &in) const {
    auto timed=getCenterRdAtTime(in.m_t);
    double d=in.getMinimumDistance(timed.first);
    if(d<timed.second) return 0;
    else return d-timed.second;
}

void xLine::makeInfinite(uint32_t dimension)
{
    m_ps.m_x=m_ps.m_y=m_ps.m_t = 1e30;
    m_pe.m_x=m_pe.m_y=m_pe.m_t = -1e30;
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const xLine& r)
{
    uint32_t i;

    os<<r.m_ps<<","<<r.m_pe;
    return os;
}

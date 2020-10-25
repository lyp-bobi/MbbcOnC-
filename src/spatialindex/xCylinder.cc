//
// Created by Chuang on 2019/10/17.
//
#define _USE_MATH_DEFINES

#include <cstring>
#include <math.h>
#include <limits>
#include <algorithm>
#include <cfloat>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
xCylinder::xCylinder(){}
xCylinder::xCylinder(const SpatialIndex::xCylinder &in) {
    m_p = in.m_p;
    m_r=in.m_r;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
}
xCylinder::xCylinder(const double *p,double r,double sTime,double eTime, uint32_t dimension){
    m_p.m_x=p[0];
    m_p.m_y = p[1];
    m_startTime=sTime;
    m_endTime=eTime;
    m_r=r;
}



xCylinder& xCylinder::operator=(const xCylinder& in)
{
    m_p = in.m_p;
    m_r=in.m_r;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
    return *this;
}
//
// IObject interface
//
xCylinder* xCylinder::clone() {
    return new xCylinder(*this);
}
//
// ISerializable interface
//
uint32_t xCylinder::getByteArraySize() const {
    return 5*sizeof(prex);
}

void xCylinder::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&(m_p.m_x), ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&(m_p.m_y), ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_r, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_startTime, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_endTime, ptr, sizeof(prex));
//    ptr += sizeof(prex);
}

void xCylinder::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;

    memcpy(ptr, &(m_p.m_x), sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &(m_p.m_y), sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_r, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_startTime, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(ptr, &m_endTime, sizeof(prex));
    ptr += sizeof(prex);
}

//
// IShape interface
//
bool xCylinder::intersectsShape(const SpatialIndex::IShape& s) const {

    const xPoint* ptp = dynamic_cast<const xPoint*>(&s);
    if (ptp != 0) return intersectsxPoint(*ptp);

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != 0) return intersectsxMBR(*pr);

//    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
//    if (ptra != 0) return ptra->intersectsShape(*this);

//    const xCylinder* pxCylinder = dynamic_cast<const xCylinder*>(&s);
//    if (pxCylinder != 0) return intersectsxCylinder(*pxCylinder);

//    const TimexMBR* ptr = dynamic_cast<const TimexMBR*>(&s);
//    if (ptr != 0) return intersectsTimexMBR(*ptr);
    throw Tools::NotSupportedException("xCylinder::intersectsShape:Not supported");
}

bool xCylinder::intersectsxPoint(const SpatialIndex::xPoint &in) const {
    if(in.m_t<m_startTime||in.m_t>m_endTime) return false;
    return m_p.getMinimumDistance(in)<m_r;
}

bool xCylinder::intersectsxMBR(const SpatialIndex::xMBR &in) const {
    //2 dimensions for space, the third for time
    if(in.m_tmin>m_endTime||in.m_tmax<m_startTime) return false;
    return m_p.getMinimumDistance(in)<=m_r+1e-7;
}
inline bool xCylinder::intersectsxCylinder(const xCylinder& in) const{throw Tools::NotSupportedException("xCylinder::intersectsxCylinder");}


bool xCylinder::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool xCylinder::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("xCylinder:touchesShape");
}
void xCylinder::getCenter(Point& out) const{
    double t=(m_startTime+m_endTime)/2;
    double p3d[3];
    p3d[0]=(m_p.m_x);
    p3d[1]=(m_p.m_y);
    p3d[2]=t;
    out=Point(p3d,3);
}
uint32_t xCylinder::getDimension() const{return 2;}
void xCylinder::getxMBR(xMBR& out) const{
    out.m_xmin=(m_p.m_x)-m_r;
    out.m_xmax=(m_p.m_x)+m_r;
    out.m_ymin=(m_p.m_y)-m_r;
    out.m_ymax=(m_p.m_y)+m_r;
    out.m_tmin=m_startTime;
    out.m_tmax=m_endTime;
}

double xCylinder::getArea() const{
    double res=M_PI*sq(m_r)*(m_endTime-m_startTime);
    return res;
}
double xCylinder::getMinimumDistance(const IShape& in) const{
    const xPoint* ptp = dynamic_cast<const xPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const xMBR* pr = dynamic_cast<const xMBR*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "xCylinder::getMinimumDistance: Not implemented yet!"
    );
}

double xCylinder::getMinimumDistance(const SpatialIndex::xMBR &in) const {
    throw Tools::NotSupportedException("");
}
double xCylinder::getMinimumDistance(const SpatialIndex::xPoint &in) const {
    double d=in.getMinimumDistance(xPoint((m_p.m_x),(m_p.m_y),0));
    if(d<m_r) return 0;
    else return d-m_r;
}

void xCylinder::makeInfinite(uint32_t dimension)
{
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}



std::ostream& SpatialIndex::operator<<(std::ostream& os, const xCylinder& r)
{
    uint32_t i;

    os<<r.m_p.m_x<<r.m_p.m_y<<r.m_r<<r.m_startTime<<r.m_endTime;
    return os;
}

/*
 * checkRel, 0 for not intersects, 1 for intersects, and 2 for timely contains
 */

int xCylinder::checkRel(const xMBR &br) const {
    if(br.m_tmin>m_endTime||br.m_tmax<m_startTime) return 0;
    bool b = intersectsxMBR(br);
    if (!b) return 0;
    double x1 = br.m_xmin, x2 =  br.m_xmax, y1 = br.m_ymin, y2 =  br.m_ymax;
    double x= (m_p.m_x), y =(m_p.m_y);
    double r2 = m_r*m_r;
    double d1 = sq(x-x1) + sq(y-y1), d2 = sq(x-x1) + sq(y-y2), d3 = sq(x-x2) + sq(y-y1), d4 = sq(x-x2) + sq(y-y2);
    if (d1<r2 && d2<r2&&d3<r2&&d4<r2)
        return 2;
    return 1;
}

static bool below0(double a,double b,double c,double ts,double te){
    if(a*ts*ts+b*ts+c<=0||a*te*te+b*te+c<=0) return true;
    double m = -b/2/a;
    if(ts<m && m<te && a*m*m +b*m +c<=0) return true;
    return false;
}

int xCylinder::checkRel(const xLine &l) const {
    double t0 = l.m_ps.m_t,  t3 = l.m_pe.m_t;
    double tlow = m_startTime, thigh = m_endTime;
    double ints = std::max(t0, tlow), inte = std::min(t3, thigh);
    auto a = l.getCenterRdAtTime(ints), b = l.getCenterRdAtTime(inte);
    double dxs=(m_p.m_x)-a.first.m_x;
    double dys=(m_p.m_y)-a.first.m_y;
    double dxe=(m_p.m_x)-b.first.m_x;
    double dye=(m_p.m_y)-b.first.m_y;
    double ts = ints, te = inte;
    double c1=sq(dxs-dxe)+sq(dys-dye),
            c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
            c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
            c4=te-ts;
    double d;
    if(c1<1e-7){
        d= std::sqrt(sq(dxs)+sq(dys));
    }else{
        double middle=-c2/c1/2;
        if(middle>ints&&middle<inte){
            d= sqrt((4*c1*c3-c2*c2)/4/c1)/c4;
        }
        else{
            d= std::min(std::sqrt(sq(dxs)+sq(dys)),std::sqrt(sq(dxe)+sq(dye)));
        }
    }
    if(d<=m_r){
        return 2;
    }
    return 0;
}



int xCylinder::checkRel(const xMBC &bc) const {
    if (bc.m_ps.m_t > m_endTime || bc.m_pe.m_t < m_startTime) return 0;
    if (m_startTime == m_endTime) {
        auto timed = bc.getCenterRdAtTime(m_startTime);
        double d = timed.first.getMinimumDistance(m_p);
        if (d + timed.second<=m_r){
            return 2;
        }else if (d <= timed.second+m_r){
            return 1;
        }else{
            return 0;
        }
    } else {
        double t0 = bc.m_ps.m_t, t1 = bc.m_ps.m_t + bc.m_rd / bc.m_rv, t2 =
                bc.m_pe.m_t - bc.m_rd / bc.m_rv, t3 = bc.m_pe.m_t;
        if (bc.m_rv < 1e-7) {
            t1 = bc.m_ps.m_t;
            t2 = bc.m_pe.m_t;
        }
        double tlow = m_startTime, thigh = m_endTime;
        double ints = std::max(t0, tlow), inte = std::min(t3, thigh);
        double ts = ints, te = inte;
        auto a = bc.getCenterRdAtTime(ints), b = bc.getCenterRdAtTime(inte);
        double d,dmid;
        double dxs=(m_p.m_x)-a.first.m_x;
        double dys=(m_p.m_y)-a.first.m_y;
        double dxe=(m_p.m_x)-b.first.m_x;
        double dye=(m_p.m_y)-b.first.m_y;


        double c1=sq(dxs-dxe)+sq(dys-dye),
                c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                c4=te-ts;
        if(c1<1e-7){
            d= std::sqrt(sq(dxs)+sq(dys));
        }else{
            double middle=-c2/c1/2;
            if(middle>ints&&middle<inte){
                d= sqrt((4*c1*c3-c2*c2)/4/c1)/c4;
            }
            else{
                d= std::min(std::sqrt(sq(dxs)+sq(dys)),std::sqrt(sq(dxe)+sq(dye)));
            }
        }
        int tmpRes=0;
        ints=std::max(t1,tlow);inte=std::min(t2,thigh);
        if(c1<1e-7){
            dmid = std::sqrt(sq(dxs)+sq(dys));
        }else{
            double middle=-c2/c1/2;
            if(middle>ints&&middle<inte){
                dmid= sqrt((4*c1*c3-c2*c2)/4/c1)/c4;
            }
            else{
                dmid= std::min(std::sqrt(sq(dxs)+sq(dys)),std::sqrt(sq(dxe)+sq(dye)));
            }
        }

        if (d > bc.m_rd + m_r) return 0;
        if (d <= m_r-bc.m_rd) return 2;

        if (dmid <= bc.m_rd + m_r) tmpRes = 1;
        //else tmpRes=0;

        double d1,d2,d3;
        //lower cone
        ints=std::max(t0,tlow);inte=std::min(t1,thigh);
        if(ints<inte){
            double c1=sq(dxs-dxe)+sq(dys-dye),
                    c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                    c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                    c4=te-ts;
            d1=c1-sq(c4*bc.m_rv);
            d2=c2+sq(c4)*2*(m_r+bc.m_rv*t0)*bc.m_rv;
            d3=c3-sq(c4*(m_r+bc.m_rv*t0));
            if(below0(d1,d2,d3,ints,inte)) return 2;
            d2 = c2 - 2*sq(c4)*(m_r - bc.m_rv*t0)*bc.m_rv;
            d3 = c3 - sq(c4*(m_r-bc.m_rv*t0));
            if(below0(d1,d2,d3,ints,inte)) tmpRes = max(1, tmpRes);
            // else 0;
        }
        //higher cone
        ints=std::max(t2,tlow);inte=std::min(t3,thigh);
        if(ints<inte){
            double c1=sq(dxs-dxe)+sq(dys-dye),
                    c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                    c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                    c4=te-ts;
            d1=c1-sq(c4*bc.m_rv);
            d2=c2-sq(c4)*2*(m_r-bc.m_rv*t3)*bc.m_rv;
            d3=c3-sq(c4*(m_r-bc.m_rv*t3));
            if(below0(d1,d2,d3,ints,inte)) return 2;
            d2=c2+sq(c4)*2*(m_r+bc.m_rv*t3)*bc.m_rv;
            d3=c3-sq(c4*(m_r+bc.m_rv*t3));
            if(below0(d1,d2,d3,ints,inte)) tmpRes = max(1, tmpRes);
            //else 0
        }
        return tmpRes;
    }
}

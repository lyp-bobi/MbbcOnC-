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

extern bool bCompactMBC = false;

using namespace SpatialIndex;
xMBC::xMBC()
{
}
xMBC::xMBC(const SpatialIndex::xMBC &in)
:m_ps(in.m_ps),m_pe(in.m_pe),m_rd(in.m_rd),m_rv(in.m_rv){
}
xMBC::xMBC(xPoint ps, xPoint pe, double rd, double rv)
:m_ps(ps),m_pe(pe),m_rd(rd),m_rv(rv){}


xMBC& xMBC::operator=(const xMBC& r)
{
    m_rd=r.m_rd;
    m_rv=r.m_rv;
    m_ps = r.m_ps;
    m_pe = r.m_pe;
    return *this;
}

bool xMBC::operator==(const SpatialIndex::xMBC &r) const {
    throw Tools::IllegalArgumentException(
            "MBC:== not supported"
    );
    return true;
}
//
// IObject interface
//
xMBC* xMBC::clone() {
    return new xMBC(*this);
}
//
// ISerializable interface
//
uint32_t xMBC::getByteArraySize() const {
    if(bCompactMBC) return 7*sizeof(prex);
    else return 8*sizeof(prex);
}

void xMBC::loadFromByteArray(const uint8_t* ptr) {
    m_ps.loadFromByteArray(ptr);
    ptr+=m_ps.getByteArraySize();
    m_pe.loadFromByteArray(ptr);
    ptr+=m_pe.getByteArraySize();
    memcpy(&m_rd, ptr, sizeof(double));
    if(!bCompactMBC) {
        ptr += sizeof(double);
        memcpy(&m_rv, ptr, sizeof(double));
//    ptr += sizeof(double);
    }
}

void xMBC::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint32_t tmp;
    m_ps.storeToByteArray(&ptr,tmp);
    ptr+=tmp;
    m_pe.storeToByteArray(&ptr,tmp);
    ptr+=tmp;
    memcpy(ptr, &m_rd, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_rv, sizeof(double));
//    ptr += sizeof(double);
    //ptr += m_dimension * sizeof(double);
}

//
// IEvolvingShape interface
//
void xMBC::getVMBR(xMBR& out) const{
}
void xMBC::getMBRAtTime(double t, SpatialIndex::xMBR &out) const {
    auto c = getCenterRdAtTime(t);
    out.m_xmin = c.first.m_x - c.second;
    out.m_xmax = c.first.m_x + c.second;
    out.m_ymin = c.first.m_y - c.second;
    out.m_ymax = c.first.m_y + c.second;
    out.m_tmin = c.first.m_t;
    out.m_tmax = c.first.m_t;
}


std::pair<xPoint,double> xMBC::getCenterRdAtTime(double t) const {
    if(t<=m_ps.m_t) return std::make_pair(xPoint(m_ps),0);
    if(t>=m_pe.m_t) return std::make_pair(xPoint(m_pe),0);
    double x=makemidmacro(m_ps.m_x,m_ps.m_t,m_pe.m_x,m_pe.m_t,t);
    double y=makemidmacro(m_ps.m_y,m_ps.m_t,m_pe.m_y,m_pe.m_t,t);
    double coord[2]={x,y};
//    xPoint *tp = xPoint::makemid(xPoint(m_pLow, m_ps.m_t, m_ps.m_t, m_dimension),
//                                      xPoint(m_pHigh, m_pe.m_t, m_pe.m_t, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double r;
    if(std::isfinite(m_rv)) {
        r = std::min(std::min((t - m_ps.m_t) * m_rv, m_rd),(m_pe.m_t-t) * m_rv);
    }else{
        r=m_rd;
    }
    xPoint tp(coord[0],coord[1],t);
    return std::make_pair(tp,r);
}

//
// IShape interface
//
bool xMBC::intersectsShape(const SpatialIndex::IShape& s) const {

    const xPoint* ptp = dynamic_cast<const xPoint*>(&s);
    if (ptp != 0) return intersectsxPoint(*ptp);

    const xMBR* pr = dynamic_cast<const xMBR*>(&s);
    if (pr != 0) return intersectsxMBR(*pr);

    const xCylinder* pcy = dynamic_cast<const xCylinder*>(&s);
    if (pcy != 0) return pcy->checkRel(*this);

//    const xTrajectory* ptra = dynamic_cast<const xTrajectory*>(&s);
//    if (ptra != 0) return ptra->intersectsShape(*this);

//    const xMBC* pxMBC = dynamic_cast<const xMBC*>(&s);
//    if (pxMBC != 0) return intersectsxMBC(*pxMBC);

//    const TimexMBR* ptr = dynamic_cast<const TimexMBR*>(&s);
//    if (ptr != 0) return intersectsTimexMBR(*ptr);
    throw Tools::NotSupportedException("xMBC::intersectsShape:Not supported");
}

bool xMBC::intersectsxPoint(const SpatialIndex::xPoint &in) const {
    if(in.m_t<m_ps.m_t||in.m_t>m_pe.m_t) return false;
    auto timed=getCenterRdAtTime(in.m_t);
    return timed.first.getMinimumDistance(in)<timed.second;
}

bool xMBC::intersectsxMBR(const SpatialIndex::xMBR &in) const {
    //2 dimensions for space, the third for time
    if(in.m_tmin>m_pe.m_t||in.m_tmax<m_ps.m_t) return false;
    if(in.m_tmin==in.m_tmax) {
        auto timed=getCenterRdAtTime(in.m_tmin);
        return timed.first.getMinimumDistance(in)<=timed.second+1e-7;
    }else{
        double t0=m_ps.m_t,t1=m_ps.m_t+m_rd/m_rv,t2=m_pe.m_t-m_rd/m_rv,t3=m_pe.m_t;
        if(m_rv<1e-7){
            t1 = t0;
            t2=t3;
        }
        double tlow=in.m_tmin,thigh=in.m_tmax;
        double ints=std::max(t0,tlow),inte=std::min(t3,thigh);
        auto a=getCenterRdAtTime(ints),b=getCenterRdAtTime(inte);
        double d=xTrajectory::line2MBRMinSED(a.first,b.first,in);
        if(d>m_rd) return false;
//        return true;
        //cylinder
        double ts=std::max(t1,tlow),te=std::min(t2,thigh);
        if(ts<=te) {
            a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            d=xTrajectory::line2MBRMinSED(a.first,b.first,in);
            if (d <= m_rd) return true;
        }

        //top
        ts=std::max(t2,tlow),te=std::min(t3,thigh);
        if(ts<=te) {
            a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            if(a.first.getMinimumDistance(in)<=a.second) return true;
            if(b.first.getMinimumDistance(in)<=b.second) return true;
            auto part = xTrajectory::cutByPhase(a.first, b.first, in);
            for (const auto &p:part) {
                int tmpsr = xTrajectory::getPhase(in, p.first, p.second);
                if(tmpsr==5) return true;
                if(tmpsr%2==0) {
                    a = getCenterRdAtTime(p.first.m_t), b = getCenterRdAtTime(p.second.m_t);
                    if(a.first.getMinimumDistance(in)<=a.second) return true;
                    if(b.first.getMinimumDistance(in)<=b.second) return true;
                }
                if(tmpsr%2==1){
                    double px, py;
                    if (tmpsr == 1 || tmpsr == 7) px = in.m_xmin;
                    else px = in.m_xmax;
                    if (tmpsr == 1 || tmpsr == 3) py = in.m_ymin;
                    else py = in.m_ymax;
                    const xPoint *s,*e;
                    if(p.first.m_t<p.second.m_t){
                        s=&p.first;
                        e=&p.second;
                    }else{
                        s=&p.second;
                        e=&p.first;
                    }
                    double _ts=s->m_t,_te=e->m_t;
                    double dxs=s->m_x-px;
                    double dys=s->m_y-py;
                    double dxe=e->m_x-px;
                    double dye=e->m_y-py;
                    double c1=sq(dxs-dxe)+sq(dys-dye),
                            c2=2*((dxe*_ts-dxs*_te)*(dxs-dxe)+(dye*_ts-dys*_te)*(dys-dye)),
                            c3=sq(dxe*_ts-dxs*te)+sq(dye*_ts-dys*te),
                            c4=_te-_ts;
                    double _a=c1*c1-c4*c4*m_rv*m_rv,
                            _b=c2,_c=c3;
                    double delta=_b*_b-4*_a*_c;
                    if(_a*_ts*_ts+_b*_ts+_c<=0||_a*_te*_te+_b*_te+_c<=0) return true;
                    if(delta>0&&_ts<-_b/2/_a&&-_b/2/_a<_te) return true;
                }
            }
        }
        //bottom
        ts=std::max(t0,tlow),te=std::min(t1,thigh);
        if(ts<=te) {
            a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            if(a.first.getMinimumDistance(in)<=a.second) return true;
            if(b.first.getMinimumDistance(in)<=b.second) return true;
            auto part = xTrajectory::cutByPhase(a.first, b.first, in);
            for (const auto &p:part) {
                int tmpsr = xTrajectory::getPhase(in, p.first, p.second);
                if(tmpsr==5) return true;
                if(tmpsr%2==0) {
                    a = getCenterRdAtTime(p.first.m_t), b = getCenterRdAtTime(p.second.m_t);
                    if(a.first.getMinimumDistance(in)<=a.second) return true;
                    if(b.first.getMinimumDistance(in)<=b.second) return true;
                }
                if(tmpsr%2==1){
                    double px, py;
                    if (tmpsr == 1 || tmpsr == 7) px = in.m_xmin;
                    else px = in.m_xmax;
                    if (tmpsr == 1 || tmpsr == 3) py = in.m_ymin;
                    else py = in.m_ymax;
                    const xPoint *s,*e;
                    if(p.first.m_t<p.second.m_t){
                        s=&p.first;
                        e=&p.second;
                    }else{
                        s=&p.second;
                        e=&p.first;
                    }
                    double _ts=s->m_t,_te=e->m_t;
                    double dxs=s->m_x-px;
                    double dys=s->m_y-py;
                    double dxe=e->m_x-px;
                    double dye=e->m_y-py;
                    double c1=sq(dxs-dxe)+sq(dys-dye),
                            c2=2*((dxe*_ts-dxs*_te)*(dxs-dxe)+(dye*_ts-dys*_te)*(dys-dye)),
                            c3=sq(dxe*_ts-dxs*te)+sq(dye*_ts-dys*te),
                            c4=_te-_ts;
                    double _a=c1*c1-c4*c4*m_rv*m_rv,
                            _b=c2,_c=c3;
                    double delta=_b*_b-4*_a*_c;
                    if(_a*_ts*_ts+_b*_ts+_c<=0||_a*_te*_te+_b*_te+_c<=0) return true;
                    if(delta>0&&_ts<-_b/2/_a&&-_b/2/_a<_te) return true;
                }
            }
        }
        return false;
//        if(t0>te||t3<ts) return false;
//        if(ts>t0&&ts<t3){
//            auto a =getCenterRdAtTime(ts);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        if(te>t0&&te<t3){
//            auto a =getCenterRdAtTime(te);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        if(t0>ts&&t0<te){
//            auto a =getCenterRdAtTime(t0);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        if(t1>ts&&t1<te){
//            auto a =getCenterRdAtTime(t1);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        if(t2>ts&&t2<te){
//            auto a =getCenterRdAtTime(t2);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        if(t3>ts&&t3<te){
//            auto a =getCenterRdAtTime(t3);
//            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
//        }
//        return false;
    }
}

inline bool xMBC::intersectsxMBC(const xMBC& in) const{throw Tools::NotSupportedException("xMBC::intersectsxMBC");}

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

bool xMBC::prevalidate(const SpatialIndex::xMBR &in) const {
    //simple
    if(2*m_rd>in.m_xmax-in.m_xmin||2*m_rd>in.m_xmax-in.m_xmin) return false;
    xMBR shrank=in;
    shrank.m_xmax-=m_rd;shrank.m_ymax-=m_rd;
    shrank.m_xmin+=m_rd;shrank.m_ymin+=m_rd;
    double ts=std::max(in.m_ymin,m_ps.m_t),te=std::min(in.m_ymax,m_pe.m_t);
    auto ps=getCenterRdAtTime(ts),pe=getCenterRdAtTime(te);
    double d=xTrajectory::line2MBRMinSED(ps.first,pe.first,shrank);
    return d==0;


//    xMBR mbr2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
//    double t0=m_ps.m_t,t1=m_ps.m_t+m_rd/m_rv,t2=m_pe.m_t-m_rd/m_rv,t3=m_pe.m_t;
//    if(m_rv<1e-7){
//        t1=t2=(m_pe.m_t+m_ps.m_t)/2;
//    }
//    xMBR rs,re,rc;
//    double xlow=in.m_pLow[0],ylow=in.m_pLow[1],tlow=in.m_pLow[2],
//            xhigh=in.m_pHigh[0],yhigh=in.m_pHigh[1],thigh=in.m_pHigh[2];
//    double ts,te;
//    //cylinder part
//    ts=std::max(t1,tlow);te=std::min(t2,thigh);
//    if(ts<=te&&xhigh-xlow>=2*m_rd&&yhigh-ylow>=2*m_rd){
//        getMBRAtTime(ts,rs);
//        getMBRAtTime(te,re);
//        rs.combinexMBR(re);
//        if(rs.intersectsxMBR(mbr2d)) {
//            double plow[3] = {xlow + m_rd, ylow + m_rd, ts},
//                    phigh[3] = {xhigh - m_rd, yhigh - m_rd, te};
//            xMBR shrank(plow, phigh, 3);
//            auto a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
//            if (xTrajectory::line2MBRMinSED(a.first, b.first, shrank) == 0) return true;
//        }
//    }
//    //top part
//    ts=std::max(t2,tlow);te=std::min(t3,thigh);
//    if(ts<=te){
//        getMBRAtTime(ts,rs);
//        getMBRAtTime(te,re);
//        if(mbr2d.containsxMBR(rs)) return true;
//        if(mbr2d.containsxMBR(re)) return true;
//        rs.getCombinedxMBR(rc,re);
//        if(rc.intersectsxMBR(mbr2d)){
//            double x1s=rs.m_pLow[0],x2s=rs.m_pHigh[0],
//                x1e=re.m_pLow[0],x2e=re.m_pHigh[0],
//                y1s=rs.m_pLow[1],y2s=rs.m_pHigh[1],
//                y1e=re.m_pLow[1],y2e=re.m_pHigh[1];
//            if((x1s<xlow&&x1e<xlow)||(x2s>xhigh&&x2e>xhigh)||(y1s<ylow&&y1e<ylow)||(y2s>yhigh&&y2e>yhigh)){}
//            else{
//                auto prd1=getIntersectPeriod(x1s,x1e,ts,te,xlow,xhigh);
//                auto prd2=getIntersectPeriod(x2s,x2e,ts,te,xlow,xhigh);
//                auto prd3=getIntersectPeriod(y1s,y1e,ts,te,ylow,yhigh);
//                auto prd4=getIntersectPeriod(y2s,y2e,ts,te,ylow,yhigh);
//                double low=std::max(std::max(prd1.first,prd2.first),std::max(prd3.first,prd4.first)),
//                        high=std::min(std::min(prd1.second,prd2.second),std::min(prd3.second,prd4.second));
//                if(low<=high) return true;
//            }
//        }
//    }
//    //bottom part
//    ts=std::max(t0,tlow);te=std::min(t1,thigh);
//    if(ts<=te){
//        getMBRAtTime(ts,rs);
//        getMBRAtTime(te,re);
//        if(mbr2d.containsxMBR(rs)) return true;
//        if(mbr2d.containsxMBR(re)) return true;
//        rs.getCombinedxMBR(rc,re);
//        if(rc.intersectsxMBR(mbr2d)){
//            double x1s=rs.m_pLow[0],x2s=rs.m_pHigh[0],
//                    x1e=re.m_pLow[0],x2e=re.m_pHigh[0],
//                    y1s=rs.m_pLow[1],y2s=rs.m_pHigh[1],
//                    y1e=re.m_pLow[1],y2e=re.m_pHigh[1];
//            if((x1s<xlow&&x1e<xlow)||(x2s>xhigh&&x2e>xhigh)||(y1s<ylow&&y1e<ylow)||(y2s>yhigh&&y2e>yhigh)){}
//            else{
//                auto prd1=getIntersectPeriod(x1s,x1e,ts,te,xlow,xhigh);
//                auto prd2=getIntersectPeriod(x2s,x2e,ts,te,xlow,xhigh);
//                auto prd3=getIntersectPeriod(y1s,y1e,ts,te,ylow,yhigh);
//                auto prd4=getIntersectPeriod(y2s,y2e,ts,te,ylow,yhigh);
//                double low=std::max(std::max(prd1.first,prd2.first),std::max(prd3.first,prd4.first)),
//                        high=std::min(std::min(prd1.second,prd2.second),std::min(prd3.second,prd4.second));
//                if(low<=high) return true;
//            }
//        }
//    }
//    return false;
}


bool xMBC::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool xMBC::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("xMBC:touchesShape");
}
void xMBC::getCenter(Point& out) const{
    double t=(m_ps.m_t+m_pe.m_t)/2;
    xMBR br;
    getMBRAtTime(t,br);
    Point p;
    br.getCenter(p);
    double p3d[3];
    p3d[0]=p.m_pCoords[0];
    p3d[1]=p.m_pCoords[1];
    p3d[2]=t;
    out=Point(p3d,3);
}
uint32_t xMBC::getDimension() const{return 3;}
void xMBC::getxMBR(xMBR& out) const{
    out.makeInfinite(2);
    if(m_rv>1e-9) {
        double ts = m_rd / m_rv;
        xMBR tmpbr;
        getMBRAtTime(m_ps.m_t + ts, tmpbr);
        out.combinexMBR(tmpbr);
        getMBRAtTime(m_pe.m_t - ts, tmpbr);
        out.combinexMBR(tmpbr);
    }
    out.combinexPoint(m_ps);
    out.combinexPoint(m_pe);
}

double xMBC::getArea() const{
    if(m_rv==0) return 0;
    double dt=m_rd/m_rv;
    if(m_rv<1e-7){
        dt=0;
    }
    double res=M_PI*sq(m_rd)*(m_pe.m_t-m_ps.m_t-2*dt)+2.0/3*M_PI*sq(m_rd);
    return res;
}
double xMBC::getMinimumDistance(const IShape& in) const{
    const xPoint* ptp = dynamic_cast<const xPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const xMBR* pr = dynamic_cast<const xMBR*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "xMBC::getMinimumDistance: Not implemented yet!"
    );
}

double xMBC::getMinimumDistance(const SpatialIndex::xMBR &in) const {
    throw Tools::NotSupportedException("");
}
double xMBC::getMinimumDistance(const SpatialIndex::xPoint &in) const {
    auto timed=getCenterRdAtTime(in.m_t);
    double d=in.getMinimumDistance(timed.first);
    if(d<timed.second) return 0;
    else return d-timed.second;
}

void xMBC::makeInfinite(uint32_t dimension)
{
    m_ps.m_x=m_ps.m_y=m_ps.m_t = 1e30;
    m_pe.m_x=m_pe.m_y=m_pe.m_t = -1e30;
    m_rd=m_rv=0;
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const xMBC& r)
{
    uint32_t i;

    os<<r.m_ps<<","<<r.m_pe<<","<<r.m_rd<<","<<r.m_rv;
    return os;
}


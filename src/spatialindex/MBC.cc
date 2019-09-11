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
//todo: solve the problem of 2 and 3 dimensions

MBC::MBC():m_dimension(0),m_pLow(0), m_pHigh(0)
{
}
MBC::MBC(const SpatialIndex::MBC &in) {
    m_dimension=in.m_dimension;
    m_rd=in.m_rd;
    m_rv=in.m_rv;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
    try
    {
        m_pLow = new double[m_dimension];
        m_pHigh = new double[m_dimension];
    }
    catch (...)
    {
        delete[] m_pLow;
        throw;
    }

    memcpy(m_pLow, in.m_pLow, (m_dimension) * sizeof(double));
    memcpy(m_pHigh, in.m_pHigh, (m_dimension) * sizeof(double));
}
MBC::MBC(const double *pLow, const double *pHigh,double sTime,double eTime, uint32_t dimension, double rd, double rv){
    m_dimension=dimension;
    try
    {
        m_pLow = new double[m_dimension];
        m_pHigh = new double[m_dimension];
    }
    catch (...)
    {
        delete[] m_pLow;
        throw;
    }

    memcpy(m_pLow, pLow, (m_dimension)* sizeof(double));
    memcpy(m_pHigh, pHigh, (m_dimension)* sizeof(double));
    m_startTime=sTime;
    m_endTime=eTime;
    m_rd=rd;
    m_rv=rv;
}

MBC::~MBC(){
    delete[] m_pLow;
    delete[] m_pHigh;
}



MBC& MBC::operator=(const MBC& r)
{
    if(this != &r)
    {
        makeInfinite(r.m_dimension);
        memcpy(m_pLow, r.m_pLow, (m_dimension) * sizeof(double));
        memcpy(m_pHigh, r.m_pHigh, (m_dimension) * sizeof(double));
    }
    m_startTime=r.m_startTime;
    m_endTime=r.m_endTime;
    m_rv=r.m_rv;
    m_rd=r.m_rd;
    return *this;
}

bool MBC::operator==(const SpatialIndex::MBC &r) const {
    if (m_dimension != r.m_dimension)
        throw Tools::IllegalArgumentException(
                "Region::operator==: Regions have different number of dimensions."
        );
    if(m_rd!=r.m_rd||m_rv!=r.m_rv) return false;
    for (uint32_t i = 0; i < m_dimension; ++i)
    {
        if (
                m_pLow[i] < r.m_pLow[i] - std::numeric_limits<double>::epsilon() ||
                m_pLow[i] > r.m_pLow[i] + std::numeric_limits<double>::epsilon() ||
                m_pHigh[i] < r.m_pHigh[i] - std::numeric_limits<double>::epsilon() ||
                m_pHigh[i] > r.m_pHigh[i] + std::numeric_limits<double>::epsilon())
            return false;
    }
    return true;
}
//
// IObject interface
//
MBC* MBC::clone() {
    return new MBC(*this);
}
//
// ISerializable interface
//
uint32_t MBC::getByteArraySize() const {
    return (sizeof(uint32_t) + 2 * (m_dimension+1) * sizeof(double)+2* sizeof(double));
}

void MBC::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    makeInfinite(m_dimension);
    memcpy(m_pLow, ptr, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(m_pHigh, ptr, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_rd, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_rv, ptr, sizeof(double));
    ptr += sizeof(double);
    //ptr += m_dimension * sizeof(double);
}

void MBC::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, m_pLow, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(ptr, m_pHigh, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_rd, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_rv, sizeof(double));
//    ptr += sizeof(double);
    //ptr += m_dimension * sizeof(double);
}

//
// IEvolvingShape interface
//
void MBC::getVMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    double v;
    for(int i=0;i<m_dimension;i++){
        v=(m_pHigh[i]-m_pLow[i])/(m_endTime-m_startTime);
        out.m_pLow[i]=v-m_rv;
        out.m_pHigh[i]=v+m_rv;
    }
}
void MBC::getMBRAtTime(double t, SpatialIndex::Region &out) const {
//    STPoint *tp = STPoint::makemid(STPoint(m_pLow, m_startTime, m_startTime, m_dimension),
//                                      STPoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double x=makemidmacro(m_pLow[0],m_startTime,m_pHigh[0],m_endTime,t);
    double y=makemidmacro(m_pLow[1],m_startTime,m_pHigh[1],m_endTime,t);
    double coord[2]={x,y};
    Point plow(coord,2),phigh(coord,2);
    if(std::isfinite(m_rv)) {
        double r = std::min(std::min((t - m_startTime) * m_rv, m_rd),(m_endTime-t) * m_rv);
        for (int i = 0; i < m_dimension; i++) {
            plow.m_pCoords[i] = coord[i] - r;
            phigh.m_pCoords[i] = coord[i] + r;
        }
    }else{
        for (int i = 0; i < m_dimension; i++) {
            plow.m_pCoords[i] = std::max(coord[i],std::min(m_pLow[i],m_pHigh[i]));
            phigh.m_pCoords[i] = std::min(coord[i] ,std::max(m_pLow[i],m_pHigh[i]));
        }
    }
    out = Region(plow, phigh);
}


std::pair<STPoint,double> MBC::getCenterRdAtTime(double t) const {
    double x=makemidmacro(m_pLow[0],m_startTime,m_pHigh[0],m_endTime,t);
    double y=makemidmacro(m_pLow[1],m_startTime,m_pHigh[1],m_endTime,t);
    double coord[2]={x,y};
//    STPoint *tp = STPoint::makemid(STPoint(m_pLow, m_startTime, m_startTime, m_dimension),
//                                      STPoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double r;
    if(std::isfinite(m_rv)) {
        r = std::min(std::min((t - m_startTime) * m_rv, m_rd),(m_endTime-t) * m_rv);
    }else{
        r=m_rd;
    }
    STPoint tp(coord,t,2);
    return std::make_pair(tp,r);
}

//
// IShape interface
//
bool MBC::intersectsShape(const SpatialIndex::IShape& s) const {

    const STPoint* ptp = dynamic_cast<const STPoint*>(&s);
    if (ptp != 0) return intersectsSTPoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

//    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
//    if (ptra != 0) return ptra->intersectsShape(*this);

//    const MBC* pMBC = dynamic_cast<const MBC*>(&s);
//    if (pMBC != 0) return intersectsMBC(*pMBC);

//    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
//    if (ptr != 0) return intersectsTimeRegion(*ptr);
    throw Tools::NotSupportedException("MBC::intersectsShape:Not supported");
}

bool MBC::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    auto timed=getCenterRdAtTime(in.m_startTime);
    return timed.first.getMinimumDistance(in)<timed.second;
}
bool MBC::intersectsSTPoint(const SpatialIndex::STPoint &in) const {
    auto timed=getCenterRdAtTime(in.m_time);
    return timed.first.getMinimumDistance(in)<timed.second;
}

bool MBC::intersectsRegion(const SpatialIndex::Region &in) const {
    //2 dimensions for space, the third for time
    if(in.m_pLow[m_dimension]>m_endTime||in.m_pHigh[m_dimension]<m_startTime) return false;
    if(in.m_pLow[m_dimension]==in.m_pHigh[m_dimension]) {
        auto timed=getCenterRdAtTime(in.m_pLow[m_dimension]);
        return timed.first.getMinimumDistance(Region(in.m_pLow,in.m_pHigh,m_dimension))<=timed.second+1e-10;
    }else{
        double t0=m_startTime,t1=m_startTime+m_rd/m_rv,t2=m_endTime-m_rd/m_rv,t3=m_endTime;
        if(m_rv<1e-7){
            t1=t2=(m_endTime+m_startTime)/2;
        }
        double tlow=in.m_pLow[in.m_dimension-1],thigh=in.m_pHigh[in.m_dimension-1];
        double ints=std::max(t0,tlow),inte=std::min(t3,thigh);
        auto a=getCenterRdAtTime(ints),b=getCenterRdAtTime(inte);
        double d=Trajectory::line2MBRMinSED(a.first,b.first,in);
        if(d>m_rd) return false;
//        return true;
        //cylinder
        double ts=std::max(t1,tlow),te=std::min(t2,thigh);
        if(ts<=te) {
            a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            d=Trajectory::line2MBRMinSED(a.first,b.first,in);
            if (d <= m_rd) return true;
        }
        Region mbr2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
        //top
        ts=std::max(t2,tlow),te=std::min(t3,thigh);
        if(ts<=te) {
            a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
            if(b.first.getMinimumDistance(mbr2d)<=b.second) return true;
            auto part = Trajectory::cutByPhase(a.first, b.first, in);
            for (const auto &p:part) {
                int tmpsr = Trajectory::getPhase(in, p.first, p.second);
                if(tmpsr==5) return true;
                if(tmpsr%2==0) {
                    a = getCenterRdAtTime(p.first.m_time), b = getCenterRdAtTime(p.second.m_time);
                    if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
                    if(b.first.getMinimumDistance(mbr2d)<=b.second) return true;
                }
                if(tmpsr%2==1){
                    double px, py;
                    if (tmpsr == 1 || tmpsr == 7) px = in.m_pLow[0];
                    else px = in.m_pHigh[0];
                    if (tmpsr == 1 || tmpsr == 3) py = in.m_pLow[1];
                    else py = in.m_pHigh[1];
                    const STPoint *s,*e;
                    if(p.first.m_time<p.second.m_time){
                        s=&p.first;
                        e=&p.second;
                    }else{
                        s=&p.second;
                        e=&p.first;
                    }
                    double _ts=s->m_time,_te=e->m_time;
                    double dxs=s->m_pCoords[0]-px;
                    double dys=s->m_pCoords[1]-py;
                    double dxe=e->m_pCoords[0]-px;
                    double dye=e->m_pCoords[1]-py;
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
            if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
            if(b.first.getMinimumDistance(mbr2d)<=b.second) return true;
            auto part = Trajectory::cutByPhase(a.first, b.first, in);
            for (const auto &p:part) {
                int tmpsr = Trajectory::getPhase(in, p.first, p.second);
                if(tmpsr==5) return true;
                if(tmpsr%2==0) {
                    a = getCenterRdAtTime(p.first.m_time), b = getCenterRdAtTime(p.second.m_time);
                    if(a.first.getMinimumDistance(mbr2d)<=a.second) return true;
                    if(b.first.getMinimumDistance(mbr2d)<=b.second) return true;
                }
                if(tmpsr%2==1){
                    double px, py;
                    if (tmpsr == 1 || tmpsr == 7) px = in.m_pLow[0];
                    else px = in.m_pHigh[0];
                    if (tmpsr == 1 || tmpsr == 3) py = in.m_pLow[1];
                    else py = in.m_pHigh[1];
                    const STPoint *s,*e;
                    if(p.first.m_time<p.second.m_time){
                        s=&p.first;
                        e=&p.second;
                    }else{
                        s=&p.second;
                        e=&p.first;
                    }
                    double _ts=s->m_time,_te=e->m_time;
                    double dxs=s->m_pCoords[0]-px;
                    double dys=s->m_pCoords[1]-py;
                    double dxe=e->m_pCoords[0]-px;
                    double dye=e->m_pCoords[1]-py;
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
inline bool MBC::intersectsMBC(const MBC& in) const{throw Tools::NotSupportedException("MBC::intersectsMBC");}

std::pair<double,double> getIntersectPeriod(
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

bool MBC::prevalidate(const SpatialIndex::Region &in) const {
    Region mbr2d(in.m_pLow,in.m_pHigh,in.m_dimension-1);
    double t0=m_startTime,t1=m_startTime+m_rd/m_rv,t2=m_endTime-m_rd/m_rv,t3=m_endTime;
    if(m_rv<1e-7){
        t1=t2=(m_endTime+m_startTime)/2;
    }
    Region rs,re,rc;
    double xlow=in.m_pLow[0],ylow=in.m_pLow[1],tlow=in.m_pLow[2],
            xhigh=in.m_pHigh[0],yhigh=in.m_pHigh[1],thigh=in.m_pHigh[2];
    double ts,te;
    //cylinder part
    ts=std::max(t1,tlow);te=std::min(t2,thigh);
    if(ts<=te&&xhigh-xlow>=2*m_rd&&yhigh-ylow>=2*m_rd){
        getMBRAtTime(ts,rs);
        getMBRAtTime(te,re);
        rs.combineRegion(re);
        if(rs.intersectsRegion(mbr2d)) {
            double plow[3] = {xlow + m_rd, ylow + m_rd, ts},
                    phigh[3] = {xhigh - m_rd, yhigh - m_rd, te};
            Region shrank(plow, phigh, 3);
            auto a = getCenterRdAtTime(ts), b = getCenterRdAtTime(te);
            if (Trajectory::line2MBRMinSED(a.first, b.first, shrank) == 0) return true;
        }
    }
    //top part
    ts=std::max(t2,tlow);te=std::min(t3,thigh);
    if(ts<=te){
        getMBRAtTime(ts,rs);
        getMBRAtTime(te,re);
        if(mbr2d.containsRegion(rs)) return true;
        if(mbr2d.containsRegion(re)) return true;
        rs.getCombinedRegion(rc,re);
        if(rc.intersectsRegion(mbr2d)){
            double x1s=rs.m_pLow[0],x2s=rs.m_pHigh[0],
                x1e=re.m_pLow[0],x2e=re.m_pHigh[0],
                y1s=rs.m_pLow[1],y2s=rs.m_pHigh[1],
                y1e=re.m_pLow[1],y2e=re.m_pHigh[1];
            if((x1s<xlow&&x1e<xlow)||(x2s>xhigh&&x2e>xhigh)||(y1s<ylow&&y1e<ylow)||(y2s>yhigh&&y2e>yhigh)){}
            else{
                auto prd1=getIntersectPeriod(x1s,x1e,ts,te,xlow,xhigh);
                auto prd2=getIntersectPeriod(x2s,x2e,ts,te,xlow,xhigh);
                auto prd3=getIntersectPeriod(y1s,y1e,ts,te,ylow,yhigh);
                auto prd4=getIntersectPeriod(y2s,y2e,ts,te,ylow,yhigh);
                double low=std::max(std::max(prd1.first,prd2.first),std::max(prd3.first,prd4.first)),
                        high=std::min(std::min(prd1.second,prd2.second),std::min(prd3.second,prd4.second));
                if(low<=high) return true;
            }
        }
    }
    //bottom part
    ts=std::max(t0,tlow);te=std::min(t1,thigh);
    if(ts<=te){
        getMBRAtTime(ts,rs);
        getMBRAtTime(te,re);
        if(mbr2d.containsRegion(rs)) return true;
        if(mbr2d.containsRegion(re)) return true;
        rs.getCombinedRegion(rc,re);
        if(rc.intersectsRegion(mbr2d)){
            double x1s=rs.m_pLow[0],x2s=rs.m_pHigh[0],
                    x1e=re.m_pLow[0],x2e=re.m_pHigh[0],
                    y1s=rs.m_pLow[1],y2s=rs.m_pHigh[1],
                    y1e=re.m_pLow[1],y2e=re.m_pHigh[1];
            if((x1s<xlow&&x1e<xlow)||(x2s>xhigh&&x2e>xhigh)||(y1s<ylow&&y1e<ylow)||(y2s>yhigh&&y2e>yhigh)){}
            else{
                auto prd1=getIntersectPeriod(x1s,x1e,ts,te,xlow,xhigh);
                auto prd2=getIntersectPeriod(x2s,x2e,ts,te,xlow,xhigh);
                auto prd3=getIntersectPeriod(y1s,y1e,ts,te,ylow,yhigh);
                auto prd4=getIntersectPeriod(y2s,y2e,ts,te,ylow,yhigh);
                double low=std::max(std::max(prd1.first,prd2.first),std::max(prd3.first,prd4.first)),
                        high=std::min(std::min(prd1.second,prd2.second),std::min(prd3.second,prd4.second));
                if(low<=high) return true;
            }
        }
    }
    return false;
}

bool MBC::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool MBC::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("MBC:touchesShape");
}
void MBC::getCenter(Point& out) const{
    double t=(m_startTime+m_endTime)/2;
    Region br;
    getMBRAtTime(t,br);
    Point p;
    br.getCenter(p);
    double p3d[3];
    p3d[0]=p.m_pCoords[0];
    p3d[1]=p.m_pCoords[1];
    p3d[2]=t;
    out=Point(p3d,3);
}
uint32_t MBC::getDimension() const{return m_dimension;}
void MBC::getMBR(Region& out) const{
    //return a 3d mbr
    Region br;
    br.makeInfinite(m_dimension);
    if(std::isfinite(m_rv)) {
        double ts = m_rd / m_rv;
        Region tmpbr;
        getMBRAtTime(m_startTime + ts, tmpbr);
        br.combineRegion(tmpbr);
        getMBRAtTime(m_endTime - ts, tmpbr);
        br.combineRegion(tmpbr);
    }
    br.combinePoint(Point(m_pLow,m_dimension));
    br.combinePoint(Point(m_pHigh,m_dimension));
    out.makeInfinite(m_dimension+1);
    double *pLow=new double[m_dimension+1];
    double *pHigh=new double[m_dimension+1];
    for(int i=0;i<m_dimension;i++){
        pLow[i]=br.m_pLow[i];
        pHigh[i]=br.m_pHigh[i];
    }
    pLow[m_dimension]=m_startTime;
    pHigh[m_dimension]=m_endTime;
    out=Region(pLow,pHigh,3);
    delete[](pLow);
    delete[](pHigh);
}

double MBC::getArea() const{
    if(m_rv==0) return 0;
    double dt=m_rd/m_rv;
    if(m_rv<1e-7){
        dt=(m_endTime-m_startTime)/2;
    }
    double res=M_PI*sq(m_rd)*(m_endTime-m_startTime-2*dt)+2.0/3*M_PI*sq(m_rd);
    return res;
}
double MBC::getMinimumDistance(const IShape& in) const{
    const STPoint* ptp = dynamic_cast<const STPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "MBC::getMinimumDistance: Not implemented yet!"
    );
}

double MBC::getMinimumDistance(const SpatialIndex::Region &in) const {
    throw Tools::NotSupportedException("");
}
double MBC::getMinimumDistance(const SpatialIndex::STPoint &in) const {
    auto timed=getCenterRdAtTime(in.m_time);
    double d=in.getMinimumDistance(timed.first);
    if(d<timed.second) return 0;
    else return d-timed.second;
}

void MBC::makeInfinite(uint32_t dimension)
{
    if (m_dimension != dimension)
    {
        delete[] m_pLow;
        delete[] m_pHigh;

        // remember that this is not a constructor. The object will be destructed normally if
        // something goes wrong (bad_alloc), so we must take care not to leave the object at an intermediate state.
        m_pLow = 0; m_pHigh = 0;

        m_dimension = dimension;
        m_pLow = new double[m_dimension];
        m_pHigh = new double[m_dimension];
    }
    for (uint32_t cIndex = 0; cIndex < m_dimension; ++cIndex)
    {
        m_pLow[cIndex] = std::numeric_limits<double>::max();
        m_pHigh[cIndex] = -std::numeric_limits<double>::max();
    }
//    m_startTime=0;
//    m_endTime=200;
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}

int MBC::getOrient() const {
    int res=0;
    for(int i=0;i<m_dimension;i++){
        if(m_pHigh[i]>m_pLow[i]) res+=int(std::pow(2,i));
    }
    return res;
}


std::ostream& SpatialIndex::operator<<(std::ostream& os, const MBC& r)
{
    uint32_t i;

    os<<"Time: "<<r.m_startTime<<" "<<r.m_endTime<<std::endl;

    os << "Low: ";
    for (i = 0; i < r.m_dimension; ++i)
    {
        os << r.m_pLow[i] << " ";
    }

    os << ", High: ";

    for (i = 0; i < r.m_dimension; ++i)
    {
        os << r.m_pHigh[i] << " ";
    }
    os<<std::endl;
    os<<"rd: "<<r.m_rd<<",rv "<<r.m_rv<<std::endl;
    return os;
}

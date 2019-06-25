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
    ptr += sizeof(double);
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
//    TimePoint *tp = TimePoint::makemid(TimePoint(m_pLow, m_startTime, m_startTime, m_dimension),
//                                      TimePoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
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


std::pair<TimePoint,double> MBC::getCenterRdAtTime(double t) const {
    double x=makemidmacro(m_pLow[0],m_startTime,m_pHigh[0],m_endTime,t);
    double y=makemidmacro(m_pLow[1],m_startTime,m_pHigh[1],m_endTime,t);
    double coord[2]={x,y};
//    TimePoint *tp = TimePoint::makemid(TimePoint(m_pLow, m_startTime, m_startTime, m_dimension),
//                                      TimePoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double r;
    if(std::isfinite(m_rv)) {
        r = std::min(std::min((t - m_startTime) * m_rv, m_rd),(m_endTime-t) * m_rv);
    }else{
        r=m_rd;
    }
    TimePoint tp(coord,t,t,2);
    return std::make_pair(tp,r);
}

//
// IShape interface
//
bool MBC::intersectsShape(const SpatialIndex::IShape& s) const {

    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&s);
    if (ptp != 0) return intersectsTimePoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
    if (ptra != 0) return ptra->intersectsShape(*this);

//    const MBC* pMBC = dynamic_cast<const MBC*>(&s);
//    if (pMBC != 0) return intersectsMBC(*pMBC);

//    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
//    if (ptr != 0) return intersectsTimeRegion(*ptr);
}

bool MBC::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    auto timed=getCenterRdAtTime(in.m_startTime);
    return timed.first.getMinimumDistance(in)<timed.second;
}
bool MBC::intersectsTimePoint(const SpatialIndex::TimePoint &in) const {
    auto timed=getCenterRdAtTime(in.m_startTime);
    return timed.first.getMinimumDistance(in)<timed.second;
}
bool MBC::intersectsRegion(const SpatialIndex::Region &in) const {
    //2 dimensions for space, the third for time
    if(in.m_pLow[m_dimension]>m_endTime||in.m_pHigh[m_dimension]<m_startTime) return false;
    if(in.m_pLow[m_dimension]==in.m_pHigh[m_dimension]) {
        auto timed=getCenterRdAtTime(in.m_pLow[m_dimension]);
        return timed.first.getMinimumDistance(Region(in.m_pLow,in.m_pHigh,m_dimension))<=timed.second+1e-10;
    }else{
        throw Tools::NotSupportedException("MBC:intersect 3DMBR");
    }
}
bool MBC::intersectsMBC(const MBC& in) const{return false;}

bool MBC::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool MBC::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("MBC:touchesShape");
}
void MBC::getCenter(Point& out) const{
    double t=(m_startTime+m_endTime)/2;
    Region br;
    getMBRAtTime(t,br);
    Point p;
    br.getCenter(out);
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
//    delete[](pLow);
//    delete[](pHigh);
}

void MBC::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    out.makeInfinite(m_dimension);
    if(std::isfinite(m_rv)) {
        double ts = m_rd / m_rv;
        Region tmpbr;
        getMBRAtTime(m_startTime + ts, tmpbr);
        out.combineRegion(tmpbr);
        getMBRAtTime(m_endTime - ts, tmpbr);
        out.combineRegion(tmpbr);
    }
    out.combinePoint(Point(m_pLow,m_dimension));
    out.combinePoint(Point(m_pHigh,m_dimension));
    out.m_startTime=m_startTime;
    out.m_endTime=m_endTime;
}

double MBC::getArea() const{
    double dt=m_rd/m_rv;
    double res=M_PI*sq(m_rd)*(m_endTime-m_startTime-2*dt)+2.0/3*M_PI*sq(m_rd);
    return res;
}
double MBC::getMinimumDistance(const IShape& in) const{
    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "Region::getMinimumDistance: Not implemented yet!"
    );
}

double MBC::getMinimumDistance(const SpatialIndex::Region &in) const {
    throw Tools::NotSupportedException("");
}
double MBC::getMinimumDistance(const SpatialIndex::TimePoint &in) const {
    auto timed=getCenterRdAtTime(in.m_startTime);
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
//[[deprecated]]
//void MBC::combineMBC(const MBC& r)
//{
////    throw Tools::NotSupportedException("don't combine MBC here,load MBCs directly to Node Level1");
//    TimeRegion br1, br2;
//    getTimeMBR(br1);
//    r.getTimeMBR(br2);
//    br1.combineRegion(br2);
//    int ori = getOrient();
//    int mod = 2;
//    double *pLow = new double[m_dimension];
//    double *pHigh = new double[m_dimension];
//    for (int i = 0; i < m_dimension; i++) {
//        if (ori % mod == 1) {
//            pLow[i] = br1.m_pLow[i];
//            pHigh[i] = br1.m_pHigh[i];
//        } else {
//            pHigh[i] = br1.m_pLow[i];
//            pLow[i] = br1.m_pHigh[i];
//        }
//        ori /= 2;
//    }
//    double stime = std::min(m_startTime, r.m_startTime),
//            etime = std::min(m_endTime, r.m_endTime);
//    TimePoint p1(pLow, stime, stime, m_dimension);
//    TimePoint p2(pHigh, etime, etime, m_dimension);
//    Point p;
//    Region tmpbr;
//    double d;
//
//    double dt=0;
//    double t1=0,t2=0;
//    if(std::isfinite(m_rv))
//        dt=m_rd/m_rv;
//    else
//        dt=m_rd/(Point(pLow,m_dimension).getMinimumDistance(Point(pHigh,m_dimension))/(m_endTime-m_startTime));
//    t1 = m_startTime + dt, t2 = m_endTime - dt;
//    double newrd = 0;
//    getMBRAtTime(m_startTime, tmpbr);
//    d = tmpbr.getMinimumDistance(*TimePoint::makemid(p1, p2, m_startTime));
//    if (d > newrd) newrd = d;
//    getMBRAtTime(t1, tmpbr);
//    tmpbr.getCenter(p);
//    d = p.getMinimumDistance(*TimePoint::makemid(p1, p2, t1)) + m_rd;
//    if (d > newrd) newrd = d;
//    getMBRAtTime(t2, tmpbr);
//    tmpbr.getCenter(p);
//    d = p.getMinimumDistance(*TimePoint::makemid(p1, p2, t2)) + m_rd;
//    if (d > newrd) newrd = d;
//    getMBRAtTime(m_endTime, tmpbr);
//    d = tmpbr.getMinimumDistance(*TimePoint::makemid(p1, p2, m_endTime));
//    if (d > newrd) newrd = d;
//
//    if(std::isfinite(r.m_rv))
//        dt=r.m_rd/r.m_rv;
//    else
//        dt=r.m_rd/(Point(pLow,r.m_dimension).getMinimumDistance(Point(pHigh,r.m_dimension))/(r.m_endTime-r.m_startTime));
//    t1=r.m_startTime+dt;t2=r.m_endTime-dt;
//    r.getMBRAtTime(r.m_startTime,tmpbr);
//    d=tmpbr.getMinimumDistance(*TimePoint::makemid(p1,p2,r.m_startTime));
//    if(d>newrd) newrd=d;
//    r.getMBRAtTime(t1,tmpbr);
//    tmpbr.getCenter(p);
//    d=p.getMinimumDistance(*TimePoint::makemid(p1,p2,t1))+r.m_rd;
//    if(d>newrd) newrd=d;
//    r.getMBRAtTime(t2,tmpbr);
//    tmpbr.getCenter(p);
//    d=p.getMinimumDistance(*TimePoint::makemid(p1,p2,t2))+r.m_rd;
//    if(d>newrd) newrd=d;
//    getMBRAtTime(r.m_endTime,tmpbr);
//    d=tmpbr.getMinimumDistance(*TimePoint::makemid(p1,p2,r.m_endTime));
//    if(d>newrd) newrd=d;
//    newrd=std::min(newrd,(Point(pLow,r.m_dimension).getMinimumDistance(Point(pHigh,r.m_dimension))));
//    m_startTime=stime;
//    m_endTime=etime;
//    m_pLow=pLow;
//    m_pHigh=pHigh;
//    m_rv=std::numeric_limits<double>::infinity();
//    m_rd=newrd;
//}
//
//bool MBC::containsMBC(const SpatialIndex::MBC &r) {
//    throw Tools::NotSupportedException("containsMBC");
//}
//
//void MBC::getCombinedMBC(MBC& out, const MBC& in) const
//{
//    out = *this;
//    out.combineMBC(in);
//}

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

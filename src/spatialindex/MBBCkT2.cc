//
// Created by Chuang on 2019/5/17.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

MBBCk::MBBCk() {
    m_k=2;
    m_startTime=std::numeric_limits<double>::max();
    m_endTime=-std::numeric_limits<double>::max();
}
MBBCk::MBBCk(int k) {
    m_k=k;
    m_startTime=std::numeric_limits<double>::max();
    m_endTime=-std::numeric_limits<double>::max();
}

MBBCk::MBBCk(const std::vector<TimePoint> tps, double tStart, double tEnd) {
    m_k=tps.size();
    m_points=tps;
    m_startTime=tStart;
    m_endTime=tEnd;
}
MBBCk::MBBCk(const SpatialIndex::MBBCk &in) {
    m_k=in.m_k;
    m_points=in.m_points;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
}

MBBCk& MBBCk::operator=(const MBBCk& r)
{
    if(this != &r)
    {
        m_k=r.m_k;
        m_points=r.m_points;
        m_startTime=r.m_startTime;
        m_endTime=r.m_endTime;
    }

    return *this;
}

bool MBBCk::operator==(const SpatialIndex::MBBCk &r) const {
    if (m_startTime < r.m_startTime - std::numeric_limits<double>::epsilon() ||
        m_startTime > r.m_startTime + std::numeric_limits<double>::epsilon() ||
        m_endTime < r.m_endTime - std::numeric_limits<double>::epsilon() ||
        m_endTime > r.m_endTime + std::numeric_limits<double>::epsilon())
        return false;
    if (m_k==r.m_k&&m_points==r.m_points)
        return false;
    return true;
}
//
// IObject interface
//
MBBCk* MBBCk::clone() {
    return new MBBCk(*this);
}
//
// ISerializable interface
//
uint32_t MBBCk::getByteArraySize() const {
    return sizeof(int)+m_points[0].getByteArraySize()*m_k+2 * sizeof(double);
}

void MBBCk::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_k, ptr, sizeof(int));
    ptr += sizeof(int);
    m_points.resize(m_k);
    for(int i=0;i<m_k;i++){
        m_points[i].loadFromByteArray(ptr);
        ptr+=m_points[i].getByteArraySize();
    }
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    //ptr += sizeof(double);
}

void MBBCk::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_k, sizeof(int));
    ptr += sizeof(int);
    for(int i=0;i<m_k;i++){
        m_points[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        ptr += tmplen;
    }
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    //ptr += sizeof(double);
    assert(len==(ptr - *data)+sizeof(double));
}

inline int MBBCk::getPhase(double t) const{
    for(int i=0;i<m_k;i++){
        if(t>m_points[i].m_startTime) return i;
    }
}
//
// IEvolvingShape interface
//
void MBBCk::getVMBR(Region& out) const{out= Region();}
void MBBCk::getMBRAtTime(double t, SpatialIndex::Region &out) const {
    TimePoint tp= m_points[getPhase(t)];
    out=Region(tp,tp);
}




//
// IShape interface
//
bool MBBCk::intersectsShape(const SpatialIndex::IShape& s) const {
    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&s);
    if (ptp != 0) return intersectsTimePoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return pr->intersectsShape(*this);

    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
    if (ptra != 0) return ptra->intersectsShape(*this);

    const MBBCk* pMBBCk = dynamic_cast<const MBBCk*>(&s);
    if (pMBBCk != 0) return intersectsMBBCk(*pMBBCk);
}

bool MBBCk::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    int p=getPhase(in.m_startTime);
    Region br;
    getMBR(br);
    if(!br.intersectsShape(in)&&m_points[p].intersectsShape(in))
        system("pause");
//    std::cout<<"Query is\n"<<in<<"\nIntersect\n"<<*this<<"\nand result is "<<m_points[p].intersectsShape(in)<<"\n\n\n\n\n";
    return m_points[p].intersectsShape(in);
}
bool MBBCk::intersectsTimePoint(const SpatialIndex::TimePoint &in) const {
    int p=getPhase(in.m_startTime);
    return m_points[p].intersectsShape(in);
}
bool MBBCk::intersectsRegion(const SpatialIndex::Region &in) const {
    for(auto mp:m_points){
        if(in.intersectsShape(mp)) return true;
    }
    return false;
}
bool MBBCk::intersectsMBBCk(const MBBCk& in) const{return false;}

bool MBBCk::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool MBBCk::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("MBBCk:touchesShape");
}
void MBBCk::getCenter(Point& out) const{
    int size=m_points.size();
    out.makeDimension(m_dimension);
    for(int i=0;i<m_dimension;i++) out.m_pCoords[i]=0;
    for(auto mbr:m_points){
        Point center;
        mbr.getCenter(center);
        for(int i=0;i<m_dimension;i++) out.m_pCoords[i]+=center.m_pCoords[i]/size;
    }
}
uint32_t MBBCk::getDimension() const{return m_dimension;}
void MBBCk::getMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    out.combinePoint(m_points[0]);
    for(auto mp:m_points) out.combinePoint(mp);
}
double MBBCk::getArea() const{
    double timeSeg=PeriodLen/m_k;
    double sum=0;
    for(auto mbr:m_points){
        sum+=mbr.getArea();
    }
    sum*=timeSeg;
    return sum;
}
double MBBCk::getMinimumDistance(const IShape& in) const{
    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "Region::getMinimumDistance: Not implemented yet!"
    );
}

double MBBCk::getMinimumDistance(const SpatialIndex::Region &in) const {
    double sum=0;
    double timeSeg=PeriodLen/m_k;
    for(auto mbr:m_points){
        sum+=mbr.getMinimumDistance(in);
    }
    sum*=timeSeg;
    return sum;
}
double MBBCk::getMinimumDistance(const SpatialIndex::TimePoint &in) const {
    Region tmp;
    getMBRAtTime(in.m_startTime,tmp);
    return tmp.getMinimumDistance(in);
}

void MBBCk::makeInfinite(uint32_t dimension,int k)
{
//    m_dimension=dimension;
    m_k=k;
    m_points.resize(m_k);
    m_vmbr.resize(m_k-1);
    for(int i=0;i<m_k;i++){
        m_points[i].makeInfinite(m_dimension);
    }
    for(int i=0;i<m_k-1;i++){
        m_vmbr[i].makeInfinite(m_dimension);
    }
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}


void MBBCk::combineMBBCk(const MBBCk& r)
{
    assert(m_k==r.m_k);
    assert(m_startTime ==r.m_startTime);
    assert(m_endTime ==r.m_endTime);
    throw Tools::NotSupportedException("combine MBBC error");
//    for(int i=0;i<m_k;i++){
//        m_points[i].combineRegion(r.m_points[i]);
//    }
}

bool MBBCk::containsMBBCk(const SpatialIndex::MBBCk &r) {
    return false;
//    if(m_k!=r.m_k||m_startTime !=r.m_startTime||m_endTime ==r.m_endTime)
//        return false;
//    for(int i=0;i<m_k;i++){
//        if(!m_points[i].containsRegion(r.m_points[i])) return  false;
//    }
//    if(m_startTime>r.m_startTime) return false;
//    if(m_endTime<r.m_endTime) return false;
//    return true;
}

void MBBCk::getCombinedMBBCk(MBBCk& out, const MBBCk& in) const
{
    out = *this;
    out.combineMBBCk(in);
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const MBBCk& r)
{
    os<<"MBR-"<<r.m_k<<":\n";

    for (auto mbr:r.m_points)
    {
        os << mbr.toString()<<"\n";
    }

    return os;
}
//
// Created by Chuang on 2019/5/15.
//

//
// Created by chuang on 4/1/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

MBRk::MBRk() {
    m_k=2;
    m_startTime=std::numeric_limits<double>::max();
    m_endTime=-std::numeric_limits<double>::max();
}
MBRk::MBRk(int k) {
    m_k=k;
    m_startTime=std::numeric_limits<double>::max();
    m_endTime=-std::numeric_limits<double>::max();
}

MBRk::MBRk(const std::vector<Region> mbrs, double tStart, double tEnd) {
    m_k=mbrs.size();
    m_mbrs=mbrs;
    m_startTime=tStart;
    m_endTime=tEnd;
}
MBRk::MBRk(const SpatialIndex::MBRk &in) {
    m_k=in.m_k;
    m_mbrs=in.m_mbrs;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
}

MBRk& MBRk::operator=(const MBRk& r)
{
    if(this != &r)
    {
        m_k=r.m_k;
        m_mbrs=r.m_mbrs;
        m_startTime=r.m_startTime;
        m_endTime=r.m_endTime;
    }

    return *this;
}

bool MBRk::operator==(const SpatialIndex::MBRk &r) const {
    if (m_startTime < r.m_startTime - std::numeric_limits<double>::epsilon() ||
        m_startTime > r.m_startTime + std::numeric_limits<double>::epsilon() ||
        m_endTime < r.m_endTime - std::numeric_limits<double>::epsilon() ||
        m_endTime > r.m_endTime + std::numeric_limits<double>::epsilon())
        return false;
    if (m_k==r.m_k&&m_mbrs==r.m_mbrs)
        return false;
    return true;
}
//
// IObject interface
//
MBRk* MBRk::clone() {
    return new MBRk(*this);
}
//
// ISerializable interface
//
uint32_t MBRk::getByteArraySize() {
    return sizeof(int)+m_mbrs[0].getByteArraySize()*m_k+2 * sizeof(double);
}

void MBRk::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_k, ptr, sizeof(int));
    ptr += sizeof(int);
    m_mbrs.resize(m_k);
    for(int i=0;i<m_k;i++){
        m_mbrs[i].loadFromByteArray(ptr);
        ptr+=m_mbrs[i].getByteArraySize();
    }
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    //ptr += sizeof(double);
}

void MBRk::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_k, sizeof(int));
    ptr += sizeof(int);
    for(int i=0;i<m_k;i++){
        m_mbrs[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        ptr += tmplen;
    }
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    //ptr += sizeof(double);
    assert(len==(ptr - *data)+sizeof(double));
}

inline int MBRk::getPhase(double t) const{
    double d=t-getPeriodStart(t);
    int p=floor(d/(PeriodLen/m_k));
    return p;
}
//
// IEvolvingShape interface
//
void MBRk::getVMBR(Region& out) const{out= Region();}
void MBRk::getMBRAtTime(double t, SpatialIndex::Region &out) const {
    out= m_mbrs[getPhase(t)];
}




//
// IShape interface
//
bool MBRk::intersectsShape(const SpatialIndex::IShape& s) const {
    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&s);
    if (ptp != 0) return intersectsTimePoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return pr->intersectsShape(*this);

    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
    if (ptra != 0) return ptra->intersectsShape(*this);

    const MBRk* pMBRk = dynamic_cast<const MBRk*>(&s);
    if (pMBRk != 0) return intersectsMBRk(*pMBRk);
}

bool MBRk::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    int p=getPhase(in.m_startTime);
    Region br;
    getMBR(br);
    if(!br.intersectsShape(in)&&m_mbrs[p].intersectsShape(in))
        system("pause");
//    std::cout<<"Query is\n"<<in<<"\nIntersect\n"<<*this<<"\nand result is "<<m_mbrs[p].intersectsShape(in)<<"\n\n\n\n\n";
    return m_mbrs[p].intersectsShape(in);
}
bool MBRk::intersectsTimePoint(const SpatialIndex::TimePoint &in) const {
    int p=getPhase(in.m_startTime);
    return m_mbrs[p].intersectsShape(in);
}
bool MBRk::intersectsRegion(const SpatialIndex::Region &in) const {
    for(auto mbr:m_mbrs){
        if(mbr.intersectsRegion(in)) return true;
    }
    return false;
}
bool MBRk::intersectsMBRk(const MBRk& in) const{return false;}

bool MBRk::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool MBRk::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("MBRk:touchesShape");
}
void MBRk::getCenter(Point& out) const{
    int size=m_mbrs.size();
    out.makeDimension(m_dimension);
    for(int i=0;i<m_dimension;i++) out.m_pCoords[i]=0;
    for(auto mbr:m_mbrs){
        Point center;
        mbr.getCenter(center);
        for(int i=0;i<m_dimension;i++) out.m_pCoords[i]+=center.m_pCoords[i]/size;
    }
}
uint32_t MBRk::getDimension() const{return m_dimension;}
void MBRk::getMBR(Region& out) const{
    out=m_mbrs[0];
    for(auto mbr:m_mbrs) out.combineRegion(mbr);
}
double MBRk::getArea() const{
    double timeSeg=PeriodLen/m_k;
    double sum=0;
    for(auto mbr:m_mbrs){
        sum+=mbr.getArea();
    }
    sum*=timeSeg;
    return sum;
}
double MBRk::getMinimumDistance(const IShape& in) const{
    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "Region::getMinimumDistance: Not implemented yet!"
    );
}

double MBRk::getMinimumDistance(const SpatialIndex::Region &in) const {
    double sum=0;
    double timeSeg=PeriodLen/m_k;
    for(auto mbr:m_mbrs){
        sum+=mbr.getMinimumDistance(in);
    }
    sum*=timeSeg;
    return sum;
}
double MBRk::getMinimumDistance(const SpatialIndex::TimePoint &in) const {
    Region tmp;
    getMBRAtTime(in.m_startTime,tmp);
    return tmp.getMinimumDistance(in);
}

void MBRk::makeInfinite(uint32_t dimension,int k)
{
//    m_dimension=dimension;
    m_k=k;
    m_mbrs.resize(m_k);
    for(int i=0;i<m_k;i++){
        m_mbrs[i].makeInfinite(m_dimension);
    }
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}


void MBRk::combineMBRk(const MBRk& r)
{
    assert(m_k==r.m_k);
    assert(m_startTime ==r.m_startTime);
    assert(m_endTime ==r.m_endTime);
    for(int i=0;i<m_k;i++){
        m_mbrs[i].combineRegion(r.m_mbrs[i]);
    }
}

bool MBRk::containsMBRk(const SpatialIndex::MBRk &r) {
    if(m_k!=r.m_k||m_startTime !=r.m_startTime||m_endTime ==r.m_endTime)
        return false;
    for(int i=0;i<m_k;i++){
        if(!m_mbrs[i].containsRegion(r.m_mbrs[i])) return  false;
    }
    if(m_startTime>r.m_startTime) return false;
    if(m_endTime<r.m_endTime) return false;
    return true;
}

void MBRk::getCombinedMBRk(MBRk& out, const MBRk& in) const
{
    out = *this;
    out.combineMBRk(in);
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const MBRk& r)
{
    os<<"MBR-"<<r.m_k<<":\n";

    for (auto mbr:r.m_mbrs)
    {
        os << mbr.toString()<<"\n";
    }

    return os;
}
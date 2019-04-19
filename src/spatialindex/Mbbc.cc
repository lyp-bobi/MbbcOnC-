//
// Created by chuang on 4/1/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

Mbbc::Mbbc() {
    m_smbr=*new Region;
    m_embr=*new Region;
    m_vmbr=*new Region;
    m_wmbr=*new Region;
    m_startTime=std::numeric_limits<double>::max();;
    m_endTime=-std::numeric_limits<double>::max();;
}

Mbbc::Mbbc(const SpatialIndex::Region &smbr, const SpatialIndex::Region &embr, const SpatialIndex::Region &vbr,
           const SpatialIndex::Region &wmbr, double tStart, double tEnd) {
    m_smbr=smbr;
    m_embr=embr;
    m_vmbr=vbr;
    m_wmbr=wmbr;
    m_startTime=tStart;
    m_endTime=tEnd;
}
Mbbc::Mbbc(const SpatialIndex::Mbbc &in) {
    m_smbr=in.m_smbr;
    m_embr=in.m_embr;
    m_vmbr=in.m_vmbr;
    m_wmbr=in.m_wmbr;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
}

Mbbc& Mbbc::operator=(const Mbbc& r)
{
    if(this != &r)
    {
        m_smbr=r.m_smbr;
        m_embr=r.m_embr;
        m_vmbr=r.m_vmbr;
        m_wmbr=r.m_wmbr;
        m_startTime=r.m_startTime;
        m_endTime=r.m_endTime;
    }

    return *this;
}

bool Mbbc::operator==(const SpatialIndex::Mbbc &r) const {
    if (m_startTime < r.m_startTime - std::numeric_limits<double>::epsilon() ||
        m_startTime > r.m_startTime + std::numeric_limits<double>::epsilon() ||
        m_endTime < r.m_endTime - std::numeric_limits<double>::epsilon() ||
        m_endTime > r.m_endTime + std::numeric_limits<double>::epsilon())
        return false;
    if (!(m_smbr==r.m_smbr)||!(m_embr==r.m_embr)||!(m_vmbr==r.m_vmbr)||!(m_wmbr==r.m_wmbr))
        return false;
    return true;
}
//
// IObject interface
//
Mbbc* Mbbc::clone() {
    return new Mbbc(*this);
}
//
// ISerializable interface
//
uint32_t Mbbc::getByteArraySize() {
    return m_smbr.getByteArraySize()+m_embr.getByteArraySize()+m_vmbr.getByteArraySize()+
        m_wmbr.getByteArraySize()+2 * sizeof(double);
}

void Mbbc::loadFromByteArray(const byte* ptr) {
    m_smbr.loadFromByteArray(ptr);
    ptr+=m_smbr.getByteArraySize();
    m_embr.loadFromByteArray(ptr);
    ptr+=m_embr.getByteArraySize();
    m_vmbr.loadFromByteArray(ptr);
    ptr+=m_vmbr.getByteArraySize();
    m_wmbr.loadFromByteArray(ptr);
    ptr+=m_wmbr.getByteArraySize();
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    //ptr += sizeof(double);
}

void Mbbc::storeToByteArray(byte **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new byte[len];
    byte* ptr = *data;
    byte* tmpb;
    u_int32_t tmplen;
    m_smbr.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    ptr += tmplen;
    m_embr.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    ptr += tmplen;
    m_vmbr.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    ptr += tmplen;
    m_wmbr.storeToByteArray(&tmpb,tmplen);
    memcpy(ptr, tmpb, tmplen);
    ptr += tmplen;
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    //ptr += sizeof(double);
    assert(len==(ptr - *data)+sizeof(double));
}
//
// IEvolvingShape interface
//
void Mbbc::getVMBR(Region& out) const{out= m_vmbr;}
void Mbbc::getMBRAtTime(double t, SpatialIndex::Region &out) const {
    if(t<m_startTime||t>m_endTime){
        throw Tools::IllegalStateException(
                "Mbbc::getMBRAtTime: time not in interval!"
                );
    }
    out.makeDimension(2);
    double xlow=std::max(m_smbr.m_pLow[0]+(t-m_startTime)*m_vmbr.m_pLow[0],
            m_embr.m_pLow[0]-(m_endTime-t)*m_vmbr.m_pHigh[0]);
    double xhigh=std::min(m_smbr.m_pHigh[0]+(t-m_startTime)*m_vmbr.m_pHigh[0],
                         m_embr.m_pHigh[0]-(m_endTime-t)*m_vmbr.m_pLow[0]);
    double ylow=std::max(m_smbr.m_pLow[1]+(t-m_startTime)*m_vmbr.m_pLow[1],
                         m_embr.m_pLow[1]-(m_endTime-t)*m_vmbr.m_pHigh[1]);
    double yhigh=std::min(m_smbr.m_pHigh[1]+(t-m_startTime)*m_vmbr.m_pHigh[1],
                          m_embr.m_pHigh[1]-(m_endTime-t)*m_vmbr.m_pLow[1]);
    out.m_pLow[0]=xlow;
    out.m_pLow[1]=ylow;
    out.m_pHigh[0]=xhigh;
    out.m_pHigh[1]=yhigh;

}




//
// IShape interface
//
bool Mbbc::intersectsShape(const SpatialIndex::IShape& s) const {
    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    const LineSegment* pls = dynamic_cast<const LineSegment*>(&s);
    if (pls != 0) return intersectsLineSegment(*pls);

    const Point* ppt = dynamic_cast<const Point*>(&s);
    if (ppt != 0) return containsPoint(*ppt);

    const Mbbc* pmbbc = dynamic_cast<const Mbbc*>(&s);
    if (pmbbc != 0) return intersectsMbbc(*pmbbc);
}

bool Mbbc::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    if(!m_wmbr.intersectsRegion(in)) return false;
    Region timed;
    getMBRAtTime(in.m_startTime,timed);
    return timed.intersectsRegion(in);
}
bool Mbbc::intersectsRegion(const Region& in) const{
    return m_wmbr.intersectsShape(in);
}
bool Mbbc::intersectsLineSegment(const LineSegment& in) const{return false;}
bool Mbbc::containsPoint(const Point& in) const{return false;}
bool Mbbc::intersectsMbbc(const Mbbc& in) const{return false;}

bool Mbbc::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool Mbbc::touchesShape(const SpatialIndex::IShape& in) const{return false;}
void Mbbc::getCenter(Point& out) const{
    m_wmbr.getCenter(out);
}
uint32_t Mbbc::getDimension() const{return 3;}
void Mbbc::getMBR(Region& out) const{out= m_wmbr;}
double Mbbc::getArea() const{ return 0;}
double Mbbc::getMinimumDistance(const IShape& in) const{return 0;}

void Mbbc::makeInfinite()
{
    m_smbr.makeInfinite(2);
    m_embr.makeInfinite(2);
    m_vmbr.makeInfinite(2);
    m_wmbr.makeInfinite(2);
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}


void Mbbc::combineMbbc(const Mbbc& r)
{
    //assert(m_startTime=r.m_startTime);
    m_smbr.combineRegion(r.m_smbr);
    m_embr.combineRegion(r.m_embr);
    m_vmbr.combineRegion(r.m_vmbr);
    m_wmbr.combineRegion(r.m_wmbr);
    m_startTime = std::min(m_startTime, r.m_startTime);
    m_endTime = std::max(m_endTime, r.m_endTime);
}

bool Mbbc::containsMbbc(const SpatialIndex::Mbbc &r) {
    if(!m_smbr.containsRegion(r.m_smbr)) return false;
    if(!m_embr.containsRegion(r.m_embr)) return false;
    if(!m_vmbr.containsRegion(r.m_vmbr)) return false;
    if(!m_wmbr.containsRegion(r.m_wmbr)) return false;
    if(m_startTime>r.m_startTime) return false;
    if(m_endTime<r.m_endTime) return false;
    return true;
}

void Mbbc::getCombinedMbbc(Mbbc& out, const Mbbc& in) const
{
    out = *this;
    out.combineMbbc(in);
}
const std::string Mbbc::toString() const{
    std::string s ="smbr:"+ std::to_string(m_smbr.m_pLow[0])+" "+
            std::to_string(m_smbr.m_pHigh[0])+" "+
            std::to_string(m_smbr.m_pLow[1])+" "+
            std::to_string(m_smbr.m_pHigh[1])+"\n"+
            "embr:"+ std::to_string(m_embr.m_pLow[0])+" "+
            std::to_string(m_embr.m_pHigh[0])+" "+
            std::to_string(m_embr.m_pLow[1])+" "+
            std::to_string(m_embr.m_pHigh[1])+"\n"+
            "vmbr:"+ std::to_string(m_vmbr.m_pLow[0])+" "+
            std::to_string(m_vmbr.m_pHigh[0])+" "+
            std::to_string(m_vmbr.m_pLow[1])+" "+
            std::to_string(m_vmbr.m_pHigh[1])+"\n"+
            "wmbr:"+ std::to_string(m_wmbr.m_pLow[0])+" "+
            std::to_string(m_wmbr.m_pHigh[0])+" "+
            std::to_string(m_wmbr.m_pLow[1])+" "+
            std::to_string(m_wmbr.m_pHigh[1])+"\n"+
            "time:"+ std::to_string(m_startTime)+" "+
            std::to_string(m_endTime);
    return s;
}
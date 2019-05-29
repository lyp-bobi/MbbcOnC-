//
// Created by Chuang on 2019/5/29.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

MBC::MBC() {
    m_dimension=2;
    makeInfinite(3);
}
MBC::MBC(const double *pLow, const double *pHigh, uint32_t dimension, double rd, double rv) {
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
    m_dimension=dimension;
    m_rd=rd;
    m_rv=rv;
}

MBC::~MBC(){
    delete[] m_pLow;
    delete[] m_pHigh;
}

MBC::MBC(const SpatialIndex::MBC &in) {
    m_dimension=in.m_dimension;
    m_rd=in.m_rd;
    m_rv=in.m_rv;
    m_startTime=in.m_dimension;
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

MBC& MBC::operator=(const MBC& r)
{
    if(this != &r)
    {
        makeInfinite(r.m_dimension);
        memcpy(m_pLow, r.m_pLow, (m_dimension) * sizeof(double));
        memcpy(m_pHigh, r.m_pHigh, (m_dimension) * sizeof(double));
    }
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
    return (sizeof(uint32_t) + 2 * m_dimension * sizeof(double)+2* sizeof(double));
}

void MBC::loadFromByteArray(const uint8_t* ptr) {
    uint32_t dimension;
    memcpy(&dimension, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(&m_rd, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_rv, ptr, sizeof(double));
    ptr += sizeof(double);
    makeInfinite(dimension);
    memcpy(m_pLow, ptr, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(m_pHigh, ptr, m_dimension * sizeof(double));
    //ptr += m_dimension * sizeof(double);
}

void MBC::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;

    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, &m_rd, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_rv, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, m_pLow, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(ptr, m_pHigh, m_dimension * sizeof(double));
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
    TimePoint tp=TimePoint::makemid(TimePoint(m_pLow,m_startTime,m_startTime,m_dimension),
                       TimePoint(m_pHigh,m_endTime,m_endTime,m_dimension),t);
    double r=std::min((t-m_startTime)*m_rv,m_rd);
    Point plow=tp,phigh=tp;
    for(int i=0;i<m_dimension;i++){
        plow.m_pCoords[i]=tp.m_pCoords[i]-r;
        phigh.m_pCoords[i]=tp.m_pCoords[i]+r;
    }
    out=Region(plow,phigh);
}


//
// IShape interface
//
bool MBC::intersectsShape(const SpatialIndex::IShape& s) const {
    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const TimePoint* ptp = dynamic_cast<const TimePoint*>(&s);
    if (ptp != 0) return intersectsTimePoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return pr->intersectsShape(*this);

    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
    if (ptra != 0) return ptra->intersectsShape(*this);

    const MBC* pMBC = dynamic_cast<const MBC*>(&s);
    if (pMBC != 0) return intersectsMBC(*pMBC);
}

bool MBC::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    Region br;
    getMBRAtTime(in.m_startTime,br);
    Point p;
    br.getCenter(p);
    double r=br.m_pHigh[0]-p.m_pCoords[0];
    return p.getMinimumDistance(in)<r;
}
bool MBC::intersectsTimePoint(const SpatialIndex::TimePoint &in) const {
    Region br;
    getMBRAtTime(in.m_startTime,br);
    Point p;
    br.getCenter(p);
    double r=br.m_pHigh[0]-p.m_pCoords[0];
    return p.getMinimumDistance(in)<r;
}
bool MBC::intersectsRegion(const SpatialIndex::Region &in) const {
    throw Tools::NotSupportedException("MBC::intersectsRegion");
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
    Region br;
    br.makeInfinite(m_dimension);
    double ts=m_rd/m_rv;
    Region tmpbr;
    getMBRAtTime(m_startTime+ts,tmpbr);
    br.combineRegion(tmpbr);
    getMBRAtTime(m_endTime-ts,tmpbr);
    br.combineRegion(tmpbr);
    br.combinePoint(Point(m_pLow,m_dimension));
    br.combinePoint(Point(m_pHigh,m_dimension));
    out.makeInfinite(m_dimension+1);
    double *pLow=new double(m_dimension+1);
    double *pHigh=new double(m_dimension+1);
    for(int i=0;i<m_dimension;i++){
        pLow[i]=br.m_pLow[i];
        pHigh[i]=br.m_pHigh[i];
    }
    pLow[m_dimension]=m_startTime;
    pHigh[m_dimension]=m_endTime;
    delete[](pLow);
    delete[](pHigh);
    out=Region(pLow,pHigh,2);
}

void MBC::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    double ts=m_rd/m_rv;
    Region tmpbr;
    getMBRAtTime(m_startTime+ts,tmpbr);
    out.combineRegion(tmpbr);
    getMBRAtTime(m_endTime-ts,tmpbr);
    out.combineRegion(tmpbr);
    out.combinePoint(Point(m_pLow,m_dimension));
    out.combinePoint(Point(m_pHigh,m_dimension));
    out.m_startTime=m_startTime;
    out.m_endTime=m_endTime;
}

double MBC::getArea() const{
    throw Tools::NotSupportedException("");
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
    Region tmp;
    getMBRAtTime(in.m_startTime,tmp);
    return tmp.getMinimumDistance(in);
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
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}


void MBC::combineMBC(const MBC& r)
{
    throw Tools::NotSupportedException("DON'T try to combine MBBC!");
}

bool MBC::containsMBC(const SpatialIndex::MBC &r) {
    throw Tools::NotSupportedException("containsMBC");
}

void MBC::getCombinedMBC(MBC& out, const MBC& in) const
{
    out = *this;
    out.combineMBC(in);
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
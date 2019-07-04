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

SBR::SBR():m_dimension(0),m_pLow(0), m_pHigh(0)
{
}
SBR::SBR(const SpatialIndex::SBR &in) {
    m_dimension=in.m_dimension;
    m_rd=in.m_rd;
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
SBR::SBR(const double *pLow, const double *pHigh,double sTime,double eTime, uint32_t dimension, double rd){
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
}

SBR::~SBR(){
    delete[] m_pLow;
    delete[] m_pHigh;
}



SBR& SBR::operator=(const SBR& r)
{
    if(this != &r)
    {
        makeInfinite(r.m_dimension);
        memcpy(m_pLow, r.m_pLow, (m_dimension) * sizeof(double));
        memcpy(m_pHigh, r.m_pHigh, (m_dimension) * sizeof(double));
    }
    m_startTime=r.m_startTime;
    m_endTime=r.m_endTime;
    m_rd=r.m_rd;
    return *this;
}

bool SBR::operator==(const SpatialIndex::SBR &r) const {
    if (m_dimension != r.m_dimension)
        throw Tools::IllegalArgumentException(
                "Region::operator==: Regions have different number of dimensions."
        );
    if(m_rd!=r.m_rd) return false;
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
SBR* SBR::clone() {
    return new SBR(*this);
}
//
// ISerializable interface
//
uint32_t SBR::getByteArraySize() const {
    return (sizeof(uint32_t) + 2 * (m_dimension+1) * sizeof(double)+ sizeof(double));
}

void SBR::loadFromByteArray(const uint8_t* ptr) {
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
    //ptr += m_dimension * sizeof(double);
}

void SBR::storeToByteArray(uint8_t **data, uint32_t &len) {
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
    //ptr += m_dimension * sizeof(double);
}

//
// IEvolvingShape interface
//
void SBR::getVMBR(Region& out) const{
    throw Tools::NotSupportedException("SBR::getVMBR:not supported");
}
void SBR::getMBRAtTime(double t, SpatialIndex::Region &out) const {
//    STPoint *tp = STPoint::makemid(STPoint(m_pLow, m_startTime, m_startTime, m_dimension),
//                                      STPoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
//    Point plow = *tp, phigh = *tp;
    double x=makemidmacro(m_pLow[0],m_startTime,m_pHigh[0],m_endTime,t);
    double y=makemidmacro(m_pLow[1],m_startTime,m_pHigh[1],m_endTime,t);
    double coord[2]={x,y};
    Point plow(coord,2),phigh(coord,2);
    for (int i = 0; i < m_dimension; i++) {
        plow.m_pCoords[i] -=m_rd;
        phigh.m_pCoords[i] +=m_rd;
    }

    out = Region(plow, phigh);
}


//
// IShape interface
//
bool SBR::intersectsShape(const SpatialIndex::IShape& s) const {

    const STPoint* ptp = dynamic_cast<const STPoint*>(&s);
    if (ptp != 0) return intersectsSTPoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
    if (ptra != 0) return ptra->intersectsShape(*this);

//    const SBR* pSBR = dynamic_cast<const SBR*>(&s);
//    if (pSBR != 0) return intersectsSBR(*pSBR);

//    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
//    if (ptr != 0) return intersectsTimeRegion(*ptr);
}

bool SBR::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    Region r;
    getMBRAtTime(in.m_startTime,r);
    return r.intersectsRegion(in);
}
bool SBR::intersectsSTPoint(const SpatialIndex::STPoint &in) const {
    Region r;
    getMBRAtTime(in.m_time,r);
    return r.containsPoint(in);
}
bool SBR::intersectsRegion(const SpatialIndex::Region &in) const {
    //2 dimensions for space, the third for time
    if(in.m_pLow[m_dimension]>m_endTime||in.m_pHigh[m_dimension]<m_startTime) return false;
    if(in.m_pLow[m_dimension]==in.m_pHigh[m_dimension]) {
        Region r;
        getMBRAtTime(in.m_pLow[m_dimension],r);
        Region in2d(in.m_pLow,in.m_pHigh,2);
        return r.intersectsRegion(in2d);
    }else{
        throw Tools::NotSupportedException("SBR:intersect 3DMBR");
    }
}

bool SBR::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool SBR::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("SBR:touchesShape");
}
bool SBR::touchesSBR(const SBR& in) const{
    throw Tools::NotSupportedException("SBR:touchesSBR");
}
void SBR::getCenter(Point& out) const{
    double t=(m_startTime+m_endTime)/2;
    Region br;
    getMBRAtTime(t,br);
    Point p;
    br.getCenter(out);
}
uint32_t SBR::getDimension() const{return m_dimension;}
void SBR::getMBR(Region& out) const{
    //return a 3d mbr
    Region br;
    br.makeInfinite(m_dimension);
    br.combinePoint(Point(m_pLow,m_dimension));
    br.combinePoint(Point(m_pHigh,m_dimension));
    out.makeInfinite(m_dimension+1);
    double *pLow=new double[m_dimension+1];
    double *pHigh=new double[m_dimension+1];
    for(int i=0;i<m_dimension;i++){
        pLow[i]=br.m_pLow[i]-m_rd;
        pHigh[i]=br.m_pHigh[i]+m_rd;
    }
    pLow[m_dimension]=m_startTime;
    pHigh[m_dimension]=m_endTime;
    out=Region(pLow,pHigh,3);
    delete[](pLow);
    delete[](pHigh);
}

double SBR::getArea() const{

    double res=4*m_rd*m_rd*(m_endTime-m_startTime);
    return res;
}
double SBR::getMinimumDistance(const IShape& in) const{
    const STPoint* ptp = dynamic_cast<const STPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);

    throw Tools::IllegalStateException(
            "Region::getMinimumDistance: Not implemented yet!"
    );
}

double SBR::getMinimumDistance(const SpatialIndex::Region &in) const {
    throw Tools::NotSupportedException("");
}
double SBR::getMinimumDistance(const SpatialIndex::STPoint &in) const {
    Region r;
    getMBRAtTime(in.m_time,r);
    return r.getMinimumDistance(in);
}

void SBR::makeInfinite(uint32_t dimension)
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

double SBR::getIntersectingArea(const SpatialIndex::SBR &in) const {
    if(m_endTime<in.m_startTime||m_startTime>in.m_endTime){
        return 0;
    }
    double ts,te;
    ts=std::max(m_startTime,in.m_startTime);
    te=std::min(m_endTime,in.m_endTime);
    Region sbr,ebr;
    Region othersbr,otherebr;
    getMBRAtTime(ts,sbr);
    in.getMBRAtTime(ts,othersbr);
    getMBRAtTime(te,ebr);
    in.getMBRAtTime(te,otherebr);
    return sbr.getIntersectingArea(othersbr)+ebr.getIntersectingArea(otherebr);
}


void SBR::getCombinedSBR(SBR& out, const SBR& in) const
{
    if (m_dimension != in.m_dimension)
        throw Tools::IllegalArgumentException(
                "Region::getCombinedRegion: Regions have different number of dimensions."
        );

    out = *this;
    out.combineSBR(in);
}

bool SBR::containsSBR(const SpatialIndex::SBR &in) {
    if(m_startTime>in.m_startTime||m_endTime<in.m_endTime)
        return false;
    Region sbr,ebr;
    Region othersbr,otherebr;
    getMBRAtTime(in.m_startTime,sbr);
    in.getMBRAtTime(in.m_startTime,othersbr);
    if(!sbr.containsRegion(othersbr)) return false;
    getMBRAtTime(in.m_endTime,ebr);
    in.getMBRAtTime(in.m_startTime,otherebr);
    if(!ebr.containsRegion(otherebr)) return false;
    return true;
}

double SBR::getMargin() const {
    Region sbr,ebr;
    getMBRAtTime(m_startTime,sbr);
    getMBRAtTime(m_endTime,ebr);
    return sbr.getMargin()+ebr.getMargin();
}

int SBR::getOrient() const {
    int res=0;
    for(int i=0;i<m_dimension;i++){
        if(m_pHigh[i]>m_pLow[i]) res+=int(std::pow(2,i));
    }
    return res;
}

void SBR::combineSBR(const SpatialIndex::SBR &in) {
//    throw Tools::NotSupportedException("this function have proved its inefficiency");
    std::cerr<<"SBR::combineSBR:this function have proved its inefficiency\n";
    if(m_dimension!=in.m_dimension){
        throw Tools::IllegalStateException("combineSBR:different dimension");
    }
    assert(m_startTime<=in.m_startTime&&m_endTime>=in.m_endTime);
    Region sbr,ebr,osbr,oebr;
    getMBRAtTime(m_startTime,sbr);
    getMBRAtTime(m_endTime,ebr);
    in.getMBRAtTime(m_startTime,osbr);
    in.getMBRAtTime(m_endTime,oebr);
    sbr.combineRegion(osbr);
    ebr.combineRegion(oebr);
    Point sct,ect;
    sbr.getCenter(sct);
    ebr.getCenter(ect);
    double sr=std::max(sbr.m_pHigh[0]-sbr.m_pLow[0],sbr.m_pHigh[1]-sbr.m_pLow[1])/2;
    double er=std::max(ebr.m_pHigh[0]-ebr.m_pLow[0],ebr.m_pHigh[1]-ebr.m_pLow[1])/2;
    m_pLow[0]=sct.m_pCoords[0];
    m_pLow[1]=sct.m_pCoords[1];
    m_pHigh[0]=ect.m_pCoords[0];
    m_pHigh[1]=ect.m_pCoords[1];
    m_rd=std::max(sr,er);
}

SBR SBR::getSBR(std::vector<SpatialIndex::SBR> &in) {
    assert(in.size()>0);
    double ts=in.front().m_startTime,te=in.front().m_endTime;
    Region smbr,embr;
    Region tr1,tr2;
    smbr.makeInfinite(in.front().m_dimension);
    embr.makeInfinite(in.front().m_dimension);
    for(const auto &sbr:in){
        sbr.getMBRAtTime(ts,tr1);
        sbr.getMBRAtTime(te,tr2);
        smbr.combineRegion(tr1);
        embr.combineRegion(tr2);
    }
    std::cerr<<"s,e MBR is"<<smbr<<"\n"<<embr<<"\n";
    Point sct,ect;
    smbr.getCenter(sct);
    embr.getCenter(ect);
    double sr=std::max(smbr.m_pHigh[0]-smbr.m_pLow[0],smbr.m_pHigh[1]-smbr.m_pLow[1])/2;
    double er=std::max(embr.m_pHigh[0]-embr.m_pLow[0],embr.m_pHigh[1]-embr.m_pLow[1])/2;
    SBR ret;
    ret.makeInfinite(in.front().m_dimension);
    ret.m_startTime=ts;
    ret.m_endTime=te;
    ret.m_pLow[0]=sct.m_pCoords[0];
    ret.m_pLow[1]=sct.m_pCoords[1];
    ret.m_pHigh[0]=ect.m_pCoords[0];
    ret.m_pHigh[1]=ect.m_pCoords[1];
    ret.m_rd=std::max(sr,er);
    return ret;
}


void SBR::combineMBC(const SpatialIndex::MBC &in) {
    assert(m_dimension==in.m_dimension);
    SBR ins;
    ins.getFromMBC(in,m_startTime,in.m_endTime);
    combineSBR(ins);
}

void SBR::getFromMBC(const SpatialIndex::MBC &in,double tstart, double tend) {
    if(in.m_startTime<tstart||in.m_endTime>tend)
        throw Tools::IllegalStateException("SBR:getFromMBC: bad time period");
    makeInfinite(in.m_dimension);
    double x1=makemidmacro(in.m_pLow[0],in.m_startTime,in.m_pHigh[0],in.m_endTime,tstart);
    double y1=makemidmacro(in.m_pLow[1],in.m_startTime,in.m_pHigh[1],in.m_endTime,tstart);
    double coord1[2]={x1,y1};
    double x2=makemidmacro(in.m_pLow[0],in.m_startTime,in.m_pHigh[0],in.m_endTime,tend);
    double y2=makemidmacro(in.m_pLow[1],in.m_startTime,in.m_pHigh[1],in.m_endTime,tend);
    double coord2[2]={x2,y2};
    for(int i=0;i<m_dimension;i++){
        m_pLow[i]=coord1[i];
        m_pHigh[i]=coord2[i];
    }
    m_rd=in.m_rd;
    m_startTime=tstart;
    m_endTime=tend;
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const SBR& r)
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
    os<<"rd: "<<r.m_rd<<std::endl;
    return os;
}

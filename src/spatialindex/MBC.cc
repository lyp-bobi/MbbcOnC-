//
// Created by Chuang on 2019/5/29.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cfloat>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

MBC::MBC(){
    m_dimension=2;
    m_pLow = new double[m_dimension];
    m_pHigh = new double[m_dimension];
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
    if(m_startTime<-1)
        std::cout<<m_startTime;
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
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(&m_rd, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_rv, ptr, sizeof(double));
    ptr += sizeof(double);
    makeInfinite(m_dimension);
    memcpy(m_pLow, ptr, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(m_pHigh, ptr, m_dimension * sizeof(double));
    //ptr += m_dimension * sizeof(double);
}

void MBC::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    ptr += sizeof(double);
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
    TimePoint tp = TimePoint::makemid(TimePoint(m_pLow, m_startTime, m_startTime, m_dimension),
                                      TimePoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
    Point plow = tp, phigh = tp;
    if(std::isfinite(m_rv)) {
        double r = std::min((t - m_startTime) * m_rv, m_rd);
        for (int i = 0; i < m_dimension; i++) {
            plow.m_pCoords[i] = tp.m_pCoords[i] - r;
            phigh.m_pCoords[i] = tp.m_pCoords[i] + r;
        }
    }else{
        for (int i = 0; i < m_dimension; i++) {
            plow.m_pCoords[i] = std::max(tp.m_pCoords[i],std::min(m_pLow[i],m_pHigh[i]));
            phigh.m_pCoords[i] = std::min(tp.m_pCoords[i] ,std::max(m_pLow[i],m_pHigh[i]));
        }
    }
    out = Region(plow, phigh);
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
    if (pr != 0) return intersectsRegion(*pr);

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
    //2 dimensions for space, the third for time
    if(in.m_pLow[m_dimension]==in.m_pHigh[m_dimension]) {
        Region br;
        getMBRAtTime(in.m_pLow[m_dimension], br);
        Point p;
        br.getCenter(p);
        double r = br.m_pHigh[0] - p.m_pCoords[0];
        br = Region(in.m_pLow, in.m_pHigh, 2);
        return p.getMinimumDistance(in) < r;
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
    double *pLow=new double(m_dimension+1);
    double *pHigh=new double(m_dimension+1);
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

int MBC::getOrient() const {
    int res=0;
    for(int i=0;i<m_dimension;i++){
        if(m_pHigh[i]>m_pLow[i]) res+=int(std::pow(2,i));
    }
    return res;
}

void MBC::combineMBC(const MBC& r)
{
//    throw Tools::NotSupportedException("don't combine MBC here,load MBCs directly to Node Level1");
    TimeRegion br1, br2;
    getTimeMBR(br1);
    r.getTimeMBR(br2);
    br1.combineRegion(br2);
    int ori = getOrient();
    int mod = 2;
    double *pLow = new double(m_dimension);
    double *pHigh = new double(m_dimension);
    for (int i = 0; i < m_dimension; i++) {
        if (ori % mod == 1) {
            pLow[i] = br1.m_pLow[i];
            pHigh[i] = br1.m_pHigh[i];
        } else {
            pHigh[i] = br1.m_pLow[i];
            pLow[i] = br1.m_pHigh[i];
        }
        ori /= 2;
    }
    double stime = std::min(m_startTime, r.m_startTime),
            etime = std::min(m_endTime, r.m_endTime);
    TimePoint p1(pLow, stime, stime, m_dimension);
    TimePoint p2(pHigh, etime, etime, m_dimension);
    Point p;
    Region tmpbr;
    double d;

    double dt=0;
    double t1=0,t2=0;
    if(std::isfinite(m_rv))
        dt=m_rd/m_rv;
    else
        dt=m_rd/(Point(pLow,m_dimension).getMinimumDistance(Point(pHigh,m_dimension))/(m_endTime-m_startTime));
    t1 = m_startTime + dt, t2 = m_endTime - dt;
    double newrd = 0;
    getMBRAtTime(m_startTime, tmpbr);
    d = tmpbr.getMinimumDistance(TimePoint::makemid(p1, p2, m_startTime));
    if (d > newrd) newrd = d;
    getMBRAtTime(t1, tmpbr);
    tmpbr.getCenter(p);
    d = p.getMinimumDistance(TimePoint::makemid(p1, p2, t1)) + m_rd;
    if (d > newrd) newrd = d;
    getMBRAtTime(t2, tmpbr);
    tmpbr.getCenter(p);
    d = p.getMinimumDistance(TimePoint::makemid(p1, p2, t2)) + m_rd;
    if (d > newrd) newrd = d;
    getMBRAtTime(m_endTime, tmpbr);
    d = tmpbr.getMinimumDistance(TimePoint::makemid(p1, p2, m_endTime));
    if (d > newrd) newrd = d;

    if(std::isfinite(r.m_rv))
        dt=r.m_rd/r.m_rv;
    else
        dt=r.m_rd/(Point(pLow,r.m_dimension).getMinimumDistance(Point(pHigh,r.m_dimension))/(r.m_endTime-r.m_startTime));
    t1=r.m_startTime+dt;t2=r.m_endTime-dt;
    r.getMBRAtTime(r.m_startTime,tmpbr);
    d=tmpbr.getMinimumDistance(TimePoint::makemid(p1,p2,r.m_startTime));
    if(d>newrd) newrd=d;
    r.getMBRAtTime(t1,tmpbr);
    tmpbr.getCenter(p);
    d=p.getMinimumDistance(TimePoint::makemid(p1,p2,t1))+r.m_rd;
    if(d>newrd) newrd=d;
    r.getMBRAtTime(t2,tmpbr);
    tmpbr.getCenter(p);
    d=p.getMinimumDistance(TimePoint::makemid(p1,p2,t2))+r.m_rd;
    if(d>newrd) newrd=d;
    getMBRAtTime(r.m_endTime,tmpbr);
    d=tmpbr.getMinimumDistance(TimePoint::makemid(p1,p2,r.m_endTime));
    if(d>newrd) newrd=d;
    newrd=std::min(newrd,(Point(pLow,r.m_dimension).getMinimumDistance(Point(pHigh,r.m_dimension))));
    m_startTime=stime;
    m_endTime=etime;
    m_pLow=pLow;
    m_pHigh=pHigh;
    m_rv=std::numeric_limits<double>::infinity();
    m_rd=newrd;
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

MBCs* MBCs::clone() {
    throw Tools::NotSupportedException("clone");
}

MBCs::MBCs(const SpatialIndex::MBCs &in) {
    m_ids=in.m_ids;
    m_mbcs=in.m_mbcs;
}
uint32_t MBCs::getByteArraySize() const {
    return sizeof(id_type)+(m_mbcs[0].getByteArraySize()+ sizeof(id_type))*m_mbcs.size();
}

void MBCs::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
    uint32_t tmplen;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    id_type size=m_mbcs.size();
    memcpy(ptr, &size, sizeof(id_type));
    ptr += sizeof(id_type);
    for(int i=0;i<size;i++){
        memcpy(ptr, &m_ids[i], sizeof(id_type));
        m_mbcs[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        if(i!=size-1)
            ptr += tmplen;
    }
}
void MBCs::loadFromByteArray(const uint8_t *ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    id_type size;
    memcpy(&size, ptr, sizeof(id_type));
    ptr += sizeof(id_type);
    m_ids.resize(size);
    m_mbcs.resize(size);
    for(int i=0;i<size;i++){
        memcpy(&m_ids[i], ptr, sizeof(id_type));
        ptr += sizeof(id_type);
        m_mbcs[i].loadFromByteArray(ptr);
        if(i!=size-1)
            ptr+=m_mbcs[i].getByteArraySize();
    }
}

//
// IEvolvingShape interface
//
void MBCs::getVMBR(Region& out) const{
    throw Tools::NotSupportedException("clone");
}
void MBCs::getMBRAtTime(double t, Region& out) const{
    throw Tools::NotSupportedException("clone");
}


//
// IShape interface
//
bool MBCs::intersectsShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
bool MBCs::containsShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
bool MBCs::touchesShape(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
void MBCs::getCenter(Point& out) const{
    double* aver=new double[m_dimension];
    for(int i=0;i<m_dimension;++i) aver[i]=0;
    Point tmpp;
    for(auto mbc:m_mbcs){
        mbc.getCenter(tmpp);
        for(int i=0;i<m_dimension;++i){
            aver[i]+=tmpp.m_pCoords[i];
        }
    }
    out.makeInfinite(m_dimension);
    for(int i=0;i<m_dimension;++i){
        out.m_pCoords[i]=aver[i]/m_ids.size();
    }
}
uint32_t MBCs::getDimension() const{
    throw Tools::NotSupportedException("clone");
}
void MBCs::getMBR(Region& out) const{
    throw Tools::NotSupportedException("clone");
}
void MBCs::getTimeMBR(SpatialIndex::TimeRegion &out) const {
    out.makeInfinite(m_dimension);
    for(auto bc:m_mbcs){
        TimeRegion tmpbr;
        bc.getTimeMBR(tmpbr);
        out.combineRegionInTime(tmpbr);
    }
}
double MBCs::getArea() const{
    throw Tools::NotSupportedException("clone");
}
double MBCs::getMinimumDistance(const IShape& in) const{
    throw Tools::NotSupportedException("clone");
}
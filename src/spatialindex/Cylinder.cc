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
//todo: solve the problem of 2 and 3 dimensions

Cylinder::Cylinder():m_dimension(0),m_p(0)
{
}
Cylinder::Cylinder(const SpatialIndex::Cylinder &in) {
    m_dimension=in.m_dimension;
    m_r=in.m_r;
    m_startTime=in.m_startTime;
    m_endTime=in.m_endTime;
    try
    {
        m_p = new double[m_dimension];
    }
    catch (...)
    {
        delete[] m_p;
        throw;
    }

    memcpy(m_p, in.m_p, (m_dimension) * sizeof(double));
}
Cylinder::Cylinder(const double *p,double r,double sTime,double eTime, uint32_t dimension){
    m_dimension=dimension;
    try
    {
        m_p = new double[m_dimension];
    }
    catch (...)
    {
        delete[] m_p;
        throw;
    }

    memcpy(m_p, p, (m_dimension)* sizeof(double));
    m_startTime=sTime;
    m_endTime=eTime;
    m_r=r;
}

Cylinder::~Cylinder(){
    delete[] m_p;
}



Cylinder& Cylinder::operator=(const Cylinder& r)
{
    if(this != &r)
    {
        makeInfinite(r.m_dimension);
        memcpy(m_p, r.m_p, (m_dimension) * sizeof(double));
    }
    m_startTime=r.m_startTime;
    m_endTime=r.m_endTime;
    m_r=r.m_r;
    return *this;
}

bool Cylinder::operator==(const SpatialIndex::Cylinder &r) const {
    if (m_dimension != r.m_dimension)
        throw Tools::IllegalArgumentException(
                "Region::operator==: Regions have different number of dimensions."
        );
    if(m_r!=r.m_r) return false;
    for (uint32_t i = 0; i < m_dimension; ++i)
    {
        if (
                m_p[i] < r.m_p[i] - std::numeric_limits<double>::epsilon() ||
                m_p[i] > r.m_p[i] + std::numeric_limits<double>::epsilon())
            return false;
    }
    return true;
}
//
// IObject interface
//
Cylinder* Cylinder::clone() {
    return new Cylinder(*this);
}
//
// ISerializable interface
//
uint32_t Cylinder::getByteArraySize() const {
    return sizeof(uint32_t) + (m_dimension+1) * sizeof(double)+sizeof(double);
}

void Cylinder::loadFromByteArray(const uint8_t* ptr) {
    memcpy(&m_dimension, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    makeInfinite(m_dimension);
    memcpy(m_p, ptr, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(&m_startTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_endTime, ptr, sizeof(double));
    ptr += sizeof(double);
    memcpy(&m_r, ptr, sizeof(double));
//    ptr += sizeof(double);
}

void Cylinder::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    memcpy(ptr, &m_dimension, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    memcpy(ptr, m_p, m_dimension * sizeof(double));
    ptr += m_dimension * sizeof(double);
    memcpy(ptr, &m_startTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_endTime, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_r, sizeof(double));
//    ptr += sizeof(double);
}

//
// IEvolvingShape interface
//
void Cylinder::getVMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    double v;
    for(int i=0;i<m_dimension;i++){
        out.m_pLow[i]=0;
        out.m_pHigh[i]=0;
    }
}
void Cylinder::getMBRAtTime(double t, SpatialIndex::Region &out) const {
//    STPoint *tp = STPoint::makemid(STPoint(m_p, m_startTime, m_startTime, m_dimension),
//                                      STPoint(m_pHigh, m_endTime, m_endTime, m_dimension), t);
//    Point plow = *tp, phigh = *tp;

//only for two dimensions
    double x=m_p[0];
    double y=m_p[1];
    double coord[2]={x,y};
    Point plow(coord,2),phigh(coord,2);
    for (int i = 0; i < m_dimension; i++) {
        plow.m_pCoords[i] = coord[i] - m_r;
        phigh.m_pCoords[i] = coord[i] + m_r;
    }
    out = Region(plow, phigh);
}


std::pair<STPoint,double> Cylinder::getCenterRdAtTime(double t) const {
    STPoint tp(m_p,t,2);
    return std::make_pair(tp,m_r);
}

//
// IShape interface
//
bool Cylinder::intersectsShape(const SpatialIndex::IShape& s) const {

    const STPoint* ptp = dynamic_cast<const STPoint*>(&s);
    if (ptp != 0) return intersectsSTPoint(*ptp);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

//    const Trajectory* ptra = dynamic_cast<const Trajectory*>(&s);
//    if (ptra != 0) return ptra->intersectsShape(*this);

//    const Cylinder* pCylinder = dynamic_cast<const Cylinder*>(&s);
//    if (pCylinder != 0) return intersectsCylinder(*pCylinder);

//    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
//    if (ptr != 0) return intersectsTimeRegion(*ptr);
    throw Tools::NotSupportedException("Cylinder::intersectsShape:Not supported");
}

bool Cylinder::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    auto timed=getCenterRdAtTime(in.m_startTime);
    return timed.first.getMinimumDistance(in)<timed.second;
}
bool Cylinder::intersectsSTPoint(const SpatialIndex::STPoint &in) const {
    if(in.m_time<m_startTime||in.m_time>m_endTime) return false;
    return Point(m_p,2).getMinimumDistance(in)<m_r;
}

bool Cylinder::intersectsRegion(const SpatialIndex::Region &in) const {
    //2 dimensions for space, the third for time
    if(in.m_pLow[m_dimension]>m_endTime||in.m_pHigh[m_dimension]<m_startTime) return false;
    return Point(m_p,2).getMinimumDistance(Region(in.m_pLow,in.m_pHigh,m_dimension))<=m_r+1e-7;
}
inline bool Cylinder::intersectsCylinder(const Cylinder& in) const{throw Tools::NotSupportedException("Cylinder::intersectsCylinder");}


bool Cylinder::containsShape(const SpatialIndex::IShape& in) const{return false;}
bool Cylinder::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException("Cylinder:touchesShape");
}
void Cylinder::getCenter(Point& out) const{
    double t=(m_startTime+m_endTime)/2;
    double p3d[3];
    p3d[0]=m_p[0];
    p3d[1]=m_p[1];
    p3d[2]=t;
    out=Point(p3d,3);
}
uint32_t Cylinder::getDimension() const{return m_dimension;}
void Cylinder::getMBR(Region& out) const{
    //return a 3d mbr
    double *pLow=new double[m_dimension+1];
    double *pHigh=new double[m_dimension+1];
    for(int i=0;i<m_dimension;i++){
        pLow[i]=m_p[i]-m_r;
        pHigh[i]=m_p[i]+m_r;
    }
    pLow[m_dimension]=m_startTime;
    pHigh[m_dimension]=m_endTime;
    out=Region(pLow,pHigh,3);
    delete[](pLow);
    delete[](pHigh);
}

double Cylinder::getArea() const{
    double res=M_PI*sq(m_r)*(m_endTime-m_startTime);
    return res;
}
double Cylinder::getMinimumDistance(const IShape& in) const{
    const STPoint* ptp = dynamic_cast<const STPoint*>(&in);
    if (ptp != 0) return getMinimumDistance(*ptp);
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::IllegalStateException(
            "Cylinder::getMinimumDistance: Not implemented yet!"
    );
}

double Cylinder::getMinimumDistance(const SpatialIndex::Region &in) const {
    throw Tools::NotSupportedException("");
}
double Cylinder::getMinimumDistance(const SpatialIndex::Point &in) const {
    double d=in.getMinimumDistance(Point(m_p,2));
    if(d<m_r) return 0;
    else return d-m_r;
}

void Cylinder::makeInfinite(uint32_t dimension)
{
    if (m_dimension != dimension)
    {
        delete[] m_p;

        // remember that this is not a constructor. The object will be destructed normally if
        // something goes wrong (bad_alloc), so we must take care not to leave the object at an intermediate state.
        m_p = 0;

        m_dimension = dimension;
        m_p = new double[m_dimension];
    }
    for (uint32_t cIndex = 0; cIndex < m_dimension; ++cIndex)
    {
        m_p[cIndex] = 0;
    }
    m_startTime = std::numeric_limits<double>::max();
    m_endTime = -std::numeric_limits<double>::max();
}



std::ostream& SpatialIndex::operator<<(std::ostream& os, const Cylinder& r)
{
    uint32_t i;

    os<<"Time: "<<r.m_startTime<<" "<<r.m_endTime<<std::endl;

    os << "Center ";
    for (i = 0; i < r.m_dimension; ++i)
    {
        os << r.m_p[i] << " ";
    }
    os<<std::endl;
    os<<"rd: "<<r.m_r<<std::endl;
    return os;
}

/*
 * checkRel, 0 for not intersects, 1 for intersects, and 2 for timely contains
 */

int Cylinder::checkRel(const Region &br) const {
    if(br.m_pLow[m_dimension]>m_endTime||br.m_pHigh[m_dimension]<m_startTime) return 0;
    bool b = intersectsRegion(br);
    if (!b) return 0;
    double x1 = br.m_pLow[0], x2 =  br.m_pHigh[0], y1 = br.m_pLow[1], y2 =  br.m_pHigh[1];
    double x= m_p[0], y =m_p[1];
    double r2 = m_r*m_r;
    double d1 = sq(x-x1) + sq(y-y1), d2 = sq(x-x1) + sq(y-y2), d3 = sq(x-x2) + sq(y-y1), d4 = sq(x-x2) + sq(y-y2);
    if (d1<r2 && d2<r2&&d3<r2&&d4<r2)
        return 2;
    return 1;
}

int Cylinder::checkRel(const MBC &bc) const {
    if (bc.m_startTime > m_endTime || bc.m_endTime < m_startTime) return 0;
    if (m_startTime == m_endTime) {
        auto timed = bc.getCenterRdAtTime(m_startTime);
        double d = timed.first.getMinimumDistance(Point(m_p,m_dimension));
        if (d + timed.second<=m_r){
            return 2;
        }else if (d <= timed.second+m_r){
            return 1;
        }else{
            return 0;
        }
    } else {
        double t0 = bc.m_startTime, t1 = bc.m_startTime + bc.m_rd / bc.m_rv, t2 =
                m_endTime - bc.m_rd / bc.m_rv, t3 = m_endTime;
        if (bc.m_rv < 1e-7) {
            t1 = bc.m_startTime;
            t2 = bc.m_endTime;
        }
        double tlow = m_startTime, thigh = m_endTime;
        double ints = std::max(t0, tlow), inte = std::min(t3, thigh);
        auto a = bc.getCenterRdAtTime(ints), b = bc.getCenterRdAtTime(inte);
        double d = Trajectory::line2lineMinSED(a.first, b.first, STPoint(m_p, ints, 2), STPoint(m_p, inte, 2));
        if (d > bc.m_rd + m_r) return 0;
        if (d <= m_r-bc.m_rd) return 2;
        double ts = ints, te = inte;
        double dxs=m_p[0]-a.first.m_pCoords[0];
        double dys=m_p[1]-a.first.m_pCoords[1];
        double dxe=m_p[0]-b.first.m_pCoords[0];
        double dye=m_p[1]-b.first.m_pCoords[1];
        //lower cone
        ints=std::max(t0,tlow);inte=std::min(t1,thigh);
        if(ints<inte){
            double c1=sq(dxs-dxe)+sq(dys-dye),
                    c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                    c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                    c4=te-ts;
            c1=c1-sq(c4*bc.m_rv);
            c2=c2+sq(c4)*2*(m_r+bc.m_rv*t0)*bc.m_rv;
            c3=c3-sq(c4*(m_r+bc.m_rv*t0));
            double middle=-c2/c1/2;
            if(middle>ints&&middle<inte){
                if(c1*middle*middle+c2*middle+c3<0){
                    return 2;
                }
            }
            else{
                if(c1*ints*ints+c2*ints+c3<=0||c1*inte*inte+c2*inte+c3<=0) return 2;
            }
        }
        //higher cone
        ints=std::max(t2,tlow);inte=std::min(t3,thigh);
        if(ints<inte){
            double c1=sq(dxs-dxe)+sq(dys-dye),
                    c2=2*((dxe*ts-dxs*te)*(dxs-dxe)+(dye*ts-dys*te)*(dys-dye)),
                    c3=sq(dxe*ts-dxs*te)+sq(dye*ts-dys*te),
                    c4=te-ts;
            c1=c1-sq(c4*bc.m_rv);
            c2=c2-sq(c4)*2*(m_r-bc.m_rv*t3)*bc.m_rv;
            c3=c3-sq(c4*(m_r-bc.m_rv*t3));
            double middle=-c2/c1/2;
            if(middle>ints&&middle<inte){
                if(c1*middle*middle+c2*middle+c3<0){
                    return 2;
                }
            }
            else{
                if(c1*ints*ints+c2*ints+c3<=0||c1*inte*inte+c2*inte+c3<=0) return 2;
            }
        }
        return 0;
    }
}


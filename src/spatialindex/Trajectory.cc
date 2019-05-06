//
// Created by chuang on 4/23/19.
//

#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;

Trajectory::Trajectory() {
}
Trajectory::Trajectory(std::vector<SpatialIndex::TimePoint>& in) {
    points=in;
}

Trajectory::Trajectory(const SpatialIndex::Trajectory &in) {
    points=in.points;
}

Trajectory& Trajectory::operator=(const Trajectory& r)
{
    if(this != &r)
    {
        points=r.points;
    }

    return *this;
}

bool Trajectory::operator==(const SpatialIndex::Trajectory &r) const {
    if (points==r.points)
        return true;
    return false;
}
//
// IObject interface
//
Trajectory* Trajectory::clone() {
    return new Trajectory(*this);
}
//
// ISerializable interface
//
uint32_t Trajectory::getByteArraySize() {
    return sizeof(unsigned long)+points[0].getByteArraySize()*points.size();
}

void Trajectory::loadFromByteArray(const byte* ptr) {
    unsigned long size;
    memcpy(&size, ptr, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    std::vector<TimePoint> p(size);
    for(int i=0;i<size;i++){
        p[i].loadFromByteArray(ptr);
        if(i!=size-1){
            ptr+=p[i].getByteArraySize();
        }
    }

}

void Trajectory::storeToByteArray(byte **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new byte[len];
    byte* ptr = *data;
    byte* tmpb;
    uint32_t tmplen;
    unsigned long size=points.size();
    memcpy(ptr, &size, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    for(int i=0;i<size;i++){
        points[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        if(i!=size-1){
            ptr += tmplen;
        }
    }

    assert(len==(ptr - *data)+tmplen);
}


//
// IShape interface
//
bool Trajectory::intersectsShape(const SpatialIndex::IShape& s) const {
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);

    const Mbbc* pbc = dynamic_cast<const Mbbc*>(&s);
    if (pbc != 0) return intersectsMbbc(*pbc);

    const TimeRegion* ptr = dynamic_cast<const TimeRegion*>(&s);
    if (ptr != 0) return intersectsTimeRegion(*ptr);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return intersectsRegion(*pr);

    const LineSegment* pls = dynamic_cast<const LineSegment*>(&s);
    if (pls != 0) return intersectsLineSegment(*pls);

    const Point* ppt = dynamic_cast<const Point*>(&s);
    if (ppt != 0) return containsPoint(*ppt);

}

//todo: check whether the link of points intersects Regions
bool Trajectory::intersectsMbbc(const SpatialIndex::Mbbc &in) const {
    for(int i=0;i<points.size()-1;i++){
        if(in.intersectsTimePoint(points[i])){
            return true;
        }
    }
    return false;
}

bool Trajectory::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    for(int i=0;i<points.size()-1;i++){
        if(points[i].intersectsShapeInTime(in)){
            return true;
        }
    }
    return false;
}
bool Trajectory::intersectsRegion(const Region& in) const{
    for(int i=0;i<points.size();i++){
        if(points[i].intersectsShape(in)){
            return true;
        }
    }
    return false;
}
bool Trajectory::intersectsLineSegment(const LineSegment& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
bool Trajectory::containsPoint(const Point& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
bool Trajectory::intersectsTrajectory(const Trajectory& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

bool Trajectory::containsShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
bool Trajectory::touchesShape(const SpatialIndex::IShape& in) const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
void Trajectory::getCenter(Point& out) const{
    throw Tools::NotSupportedException("not supported now");
}
uint32_t Trajectory::getDimension() const{return 2;}
void Trajectory::getMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    for(int i=0;i<points.size();i++){
        out.combinePoint(points[i]);
    }
}
void Trajectory::getMbbc(Mbbc& out) const{
    out.makeInfinite();

    double startx=points.begin()->m_pCoords[0],starty=points.begin()->m_pCoords[1],startt=points.begin()->m_startTime;
    double endx=points.back().m_pCoords[0],endy=points.back().m_pCoords[1],endt=points.back().m_startTime;
    double maxvxP=0,maxvxN=0,maxvyP=0,maxvyN=0;
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<points.size();i++){
        if(points[i].m_startTime!=startt){
            double vx=(points[i].m_pCoords[0]-startx)/(points[i].m_startTime-startt);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;

            double vy=(points[i].m_pCoords[1]-starty)/(points[i].m_startTime-startt);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(points[i].m_startTime!=endt){
            double vx=(endx-points[i].m_pCoords[0])/(endt-points[i].m_startTime);
            if(vx>maxvxP) maxvxP=vx;
            if(vx>1)
                std::cout<<"WARNING"<<points[i].toString()<<std::endl<<endx
                <<" "<<endt<<std::endl;
            if(vx<maxvxN) maxvxN=vx;
            double vy=(endy-points[i].m_pCoords[1])/(endt-points[i].m_startTime);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }

        if(points[i].m_pCoords[0]<minx) minx=points[i].m_pCoords[0];
        if(points[i].m_pCoords[0]>maxx) maxx=points[i].m_pCoords[0];
        if(points[i].m_pCoords[1]<miny) miny=points[i].m_pCoords[1];
        if(points[i].m_pCoords[1]>maxy) maxy=points[i].m_pCoords[1];
    }
    double sLow[2]={startx,starty};
    double sHigh[2]={startx,starty};
    double eLow[2]={endx,endy};
    double eHigh[2]={endx,endy};
    double vLow[2]={maxvxN,maxvyN};
    double vHigh[2]={maxvxP,maxvyP};
    double wLow[2]={minx,miny};
    double wHigh[2]={maxx,maxy};
    double stime=int(startt/PeriodLen)*PeriodLen;
    out= Mbbc(Region(sLow,sHigh,2),Region(eLow,eHigh,2),
                Region(vLow,vHigh,2),Region(wLow,wHigh,2),stime,stime+PeriodLen);

}
double Trajectory::getArea() const{ return 0;}
double Trajectory::getMinimumDistance(const IShape& s) const{
    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);

    const Mbbc* pbc = dynamic_cast<const Mbbc*>(&s);
    if (pbc != 0) return getMinimumDistance(*pbc);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return getMinimumDistance(*pr);

    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

double Trajectory::getMinimumDistance(const SpatialIndex::Region &in) const {
    double sum=0;
    sum+=points[0].getMinimumDistance(in)*(points[1].m_startTime-points[0].m_startTime);
    for(int i=1;i<points.size()-1;i++){
        sum+=points[i].getMinimumDistance(in)*(points[i+1].m_startTime-points[i-1].m_startTime);
    }
    sum+=points[points.size()-1].getMinimumDistance(in)*(points[points.size()-1].m_startTime-points[points.size()-2].m_startTime);
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::Mbbc &in) const {
    double sum=0;
    sum+=in.getMinimumDistance(points[0])*(points[1].m_startTime-points[0].m_startTime);
    for(int i=1;i<points.size()-1;i++){
        sum+=in.getMinimumDistance(points[i])*(points[i+1].m_startTime-points[i-1].m_startTime);
    }
    sum+=in.getMinimumDistance(points[points.size()-1])*(points[points.size()-1].m_startTime-points[points.size()-2].m_startTime);
    sum=sum/2;
    return sum;
}

double Trajectory::getMinimumDistance(const SpatialIndex::Trajectory &in) const {
    //NOTICE: this function is not symmetry!
    int cursor=0;
    double time1,time2;
    time1=points[0].m_startTime;
    time2=in.points[0].m_startTime;
    double dist;
    double sum=0;
    while(cursor<in.points.size()-1&&time1>in.points[cursor+1].m_startTime){
        cursor++;
        time2=in.points[cursor].m_startTime;
    }
    if(time1<time2||cursor==in.points.size()-1){//if no corresponding data
        dist=points[0].getMinimumDistance(in.points[cursor]);
    }else{
        dist=points[0].getMinimumDistance(TimePoint::makemid(in.points[cursor],in.points[cursor+1],time1));
    }
    sum+=dist*(points[1].m_startTime-points[0].m_startTime);
    for(int i=1;i<points.size()-1;i++){
        time1=points[i].m_startTime;
        while(cursor<in.points.size()-1&&time1>in.points[cursor+1].m_startTime){
            cursor++;
            time2=in.points[cursor].m_startTime;
        }
        if(time1<time2||cursor==in.points.size()-1){//if no corresponding data
            dist=points[i].getMinimumDistance(in.points[cursor]);
        }else{
            dist=points[i].getMinimumDistance(TimePoint::makemid(in.points[cursor],in.points[cursor+1],time1));
        }
        sum+=dist*(points[i+1].m_startTime-points[i-1].m_startTime);
    }
    time1=points[points.size()-1].m_startTime;
    while(cursor<in.points.size()-1&&time1>in.points[cursor+1].m_startTime){
        cursor++;
        time2=in.points[cursor].m_startTime;
    }
    if(time1<time2||cursor==in.points.size()-1){//if no corresponding data
        dist=points[points.size()-1].getMinimumDistance(in.points[cursor]);
    }else{
        dist=points[points.size()-1].getMinimumDistance(TimePoint::makemid(in.points[cursor],in.points[cursor+1],time1));
    }
    sum+=dist*(points[points.size()-1].m_startTime-points[points.size()-2].m_startTime);
    sum=sum/2;
    return sum;
}

void Trajectory::makeInfinite()
{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}


void Trajectory::combineTrajectory(const Trajectory& r)
{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

bool Trajectory::containsTrajectory(const SpatialIndex::Trajectory &r) {
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

void Trajectory::getCombinedTrajectory(Trajectory& out, const Trajectory& in) const
{
    out = *this;
    out.combineTrajectory(in);
}
const std::string Trajectory::toString() const{
    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}
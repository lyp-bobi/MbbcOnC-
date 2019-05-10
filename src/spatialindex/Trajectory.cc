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

void Trajectory::loadFromByteArray(const uint8_t* ptr) {
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
    points.clear();
    points=p;
}

void Trajectory::storeToByteArray(uint8_t **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new uint8_t[len];
    uint8_t* ptr = *data;
    uint8_t* tmpb;
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

TimePoint Trajectory::getPointAtTime(const double time) const {
    if(time<points.front().m_startTime||time>points.back().m_endTime){
        throw Tools::IllegalArgumentException(
                "Trajectory::getPointAtTime: time"+std::to_string(time)+"is illegal."
                );}
    if(points.size()==1){
        return points[0];
    }
    auto pre =points.begin(),next=points.begin();
    next++;
    while(next->m_startTime<time&&next!=points.end()){
        pre++;next++;
    }
    double h1= (time-pre->m_startTime)/(next->m_startTime-pre->m_startTime);
    double h2= (next->m_startTime-time)/(next->m_startTime-pre->m_startTime);
    double *coords= new double(m_dimension);
    for (int i = 0; i < m_dimension; ++i) {
        coords[i]=h2*pre->m_pCoords[i]+h1*next->m_pCoords[i];
    }
    return TimePoint(coords,time,time,m_dimension);
}

bool Trajectory::intersectsMbbc(const SpatialIndex::Mbbc &in) const {
    for(int i=0;i<points.size()-1;i++){
        if(in.intersectsTimePoint(points[i])){
            return true;
        }
    }
    return false;
}

bool Trajectory::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    if(in.m_startTime==in.m_endTime){//time slice
        if(in.m_startTime<points.front().m_startTime||in.m_startTime>points.back().m_endTime){
            return false;
        }
        TimePoint tp=getPointAtTime(in.m_startTime);
        return tp.intersectsShape(in);
    }else{
        throw Tools::NotSupportedException("time interval range not supported");
    }
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
//todo: should have implement a time divivsion class!
void Trajectory::getMbbc(Mbbc& out) const{
    out.makeInfinite();

    double startx=points.begin()->m_pCoords[0],starty=points.begin()->m_pCoords[1],startt=points.begin()->m_startTime;
    double endx=points.back().m_pCoords[0],endy=points.back().m_pCoords[1],endt=points.back().m_startTime;
    double maxvxP=-std::numeric_limits<double>::max(),
        maxvxN=std::numeric_limits<double>::max(),
        maxvyP=-std::numeric_limits<double>::max(),
        maxvyN=std::numeric_limits<double>::max();
    double minx=startx,maxx=startx,miny=starty,maxy=starty;
    for(int i=0;i<points.size();i++){
        if(points[i].m_startTime-startt>0){
            double vx=(points[i].m_pCoords[0]-startx)/(points[i].m_startTime-startt);
            if(vx>maxvxP) maxvxP=vx;
            if(vx<maxvxN) maxvxN=vx;

            double vy=(points[i].m_pCoords[1]-starty)/(points[i].m_startTime-startt);
            if(vy>maxvyP) maxvyP=vy;
            if(vy<maxvyN) maxvyN=vy;
        }
        if(endt-points[i].m_startTime>0){
            double vx=(endx-points[i].m_pCoords[0])/(endt-points[i].m_startTime);
            if(vx>maxvxP) maxvxP=vx;
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
    double nstartx,nendx,nstarty,nendy;
    nstartx=startx-(endx-startx)/(endt-startt)*startt;
    nstarty=starty-(endy-starty)/(endt-startt)*startt;
    nendx=endx+(endx-startx)/(endt-startt)*(PeriodLen-endt);
    nendy=endy+(endy-starty)/(endt-startt)*(PeriodLen-endt);
    double sLow[2]={nstartx,nstarty};
    double sHigh[2]={nstartx,nstarty};
    double eLow[2]={nendx,nendy};
    double eHigh[2]={nendx,nendy};
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
    std::string s;
    s="Trajectory length:"+std::to_string(points.size())+"\n"+
            "points are"+"\n";
    for(auto p:points){
        s+=std::to_string(p.m_pCoords[0])+","+std::to_string(p.m_pCoords[1])+
                ","+std::to_string(p.m_startTime)+"\t";
    }
    s+="\n";
    return s;
}
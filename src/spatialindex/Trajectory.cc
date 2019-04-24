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
    unsigned long len;
    memcpy(&len, ptr, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    std::vector<TimePoint> p(len);
    for(int i=0;i<len;i++){
        p[i].loadFromByteArray(ptr);
        if(i!=len-1){
            ptr+=p[i].getByteArraySize();
        }
    }

}

void Trajectory::storeToByteArray(byte **data, uint32_t &len) {
    len = getByteArraySize();
    *data = new byte[len];
    byte* ptr = *data;
    byte* tmpb;
    u_int32_t tmplen;
    unsigned long size=points.size();
    memcpy(ptr, &size, sizeof(unsigned long));
    ptr += sizeof(unsigned long);
    for(int i=0;i<len;i++){
        points[i].storeToByteArray(&tmpb,tmplen);
        memcpy(ptr, tmpb, tmplen);
        if(i!=len-1){
            ptr += tmplen;
        }
    }

    assert(len==(ptr - *data)+tmplen);
}


//
// IShape interface
//
bool Trajectory::intersectsShape(const SpatialIndex::IShape& s) const {
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

    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);
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
double Trajectory::getArea() const{ return 0;}
double Trajectory::getMinimumDistance(const IShape& s) const{
    const Mbbc* pbc = dynamic_cast<const Mbbc*>(&s);
    if (pbc != 0) return getMinimumDistance(*pbc);

    const Region* pr = dynamic_cast<const Region*>(&s);
    if (pr != 0) return getMinimumDistance(*pr);

    const Trajectory* pTrajectory = dynamic_cast<const Trajectory*>(&s);
    if (pTrajectory != 0) return intersectsTrajectory(*pTrajectory);

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
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

bool Trajectory::intersectsTimeRegion(const SpatialIndex::TimeRegion &in) const {
    for(int i=0;i<points.size()-1;i++){
        if(MovingPoint(points[i],points.at(i+1),points[i].m_startTime,points.at(i+1).m_startTime).intersectsShapeInTime(in)){
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
uint32_t Trajectory::getDimension() const{return 3;}
void Trajectory::getMBR(Region& out) const{
    out.makeInfinite(m_dimension);
    for(int i=0;i<points.size();i++){
        out.combinePoint(points[i]);
    }
}
double Trajectory::getArea() const{ return 0;}
double Trajectory::getMinimumDistance(const IShape& in) const{
    const Region* pr = dynamic_cast<const Region*>(&in);
    if (pr != 0) return getMinimumDistance(*pr);


    throw Tools::NotSupportedException(
            "Trajectory::getMinimumDistance: Not implemented yet!"
    );
}

double Trajectory::getMinimumDistance(const SpatialIndex::Region &in) const {
    double sum=0;
    for(int i=0;i<points.size();i++){
        sum+=points[i].getMinimumDistance(in);
    }
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
//
// Created by chuang on 4/23/19.
//


#pragma once

namespace SpatialIndex
{
    class SIDX_DLL Trajectory: public Tools::IObject, public virtual IShape{

    public:
    Trajectory();
    Trajectory(std::vector<TimePoint> &in);
    Trajectory(const Trajectory& in);
    virtual Trajectory &operator=(const Trajectory &r);

    virtual bool operator==(const Trajectory &r) const;


    //
    // IObject interface
    //
    virtual Trajectory *clone();

    //
    // ISerializable interface
    //
    virtual uint32_t getByteArraySize();

    virtual void loadFromByteArray(const byte *data);

    virtual void storeToByteArray(byte **data, uint32_t &len);



    //
    // IShape interface
    //
    virtual bool intersectsShape(const IShape& in) const;
    virtual bool containsShape(const IShape& in) const;
    virtual bool touchesShape(const IShape& in) const;
    virtual void getCenter(Point& out) const;
    virtual uint32_t getDimension() const;
    virtual void getMBR(Region& out) const;
    virtual double getArea() const;
    virtual double getMinimumDistance(const IShape& in) const;
    virtual double getMinimumDistance(const Region& in) const;
    virtual double getMinimumDistance(const Mbbc& in) const;
    virtual double getMinimumDistance(const Trajectory& in) const;

    virtual bool intersectsTimeRegion(const TimeRegion& in) const;
    virtual bool intersectsRegion(const Region& in) const;
    virtual bool intersectsMbbc(const Mbbc& in) const;
    virtual bool intersectsLineSegment(const LineSegment& in) const;
    virtual bool containsPoint(const Point& in) const;
    virtual bool intersectsTrajectory(const Trajectory& in) const;

    virtual void combineTrajectory(const Trajectory& in);
    virtual bool containsTrajectory(const Trajectory& in);
    virtual void getCombinedTrajectory(Trajectory& out, const Trajectory& in) const;

    const std::string toString() const;


    std::vector<TimePoint> points;

    friend SIDX_DLL std::ostream operator<<(std::ostream os,const MovingRegion &r);

    static const uint32_t m_dimension=2;
    virtual void makeInfinite();
    private:
};
typedef Tools::PoolPointer<Trajectory> TrajectoryPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Trajectory& r);
}

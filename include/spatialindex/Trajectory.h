//
// Created by chuang on 4/23/19.
//


#pragma once

#include "ShapeList.h"

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
    virtual uint32_t getByteArraySize() const;

    virtual void loadFromByteArray(const uint8_t *data);

    virtual void storeToByteArray(uint8_t **data, uint32_t &len);



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
    virtual double getMinimumDistance(const TimePoint& in) const;
    virtual double getMinimumDistance(const Region& in) const;
    virtual double getMinimumDistance(const MBC& in) const;
    virtual double getMinimumDistance(const Trajectory& in) const;
    virtual double getMinimumDistance(const ShapeList& in) const;

    virtual bool intersectsTimeRegion(const TimeRegion& in) const;
    virtual bool intersectsRegion(const Region& in) const;
    virtual bool intersectsTrajectory(const Trajectory& in) const;

    virtual void combineTrajectory(const Trajectory& in);
    virtual bool containsTrajectory(const Trajectory& in);
    virtual void getCombinedTrajectory(Trajectory& out, const Trajectory& in) const;

    virtual void getMBC(MBC& out) const;
    virtual void getMBRfull(Region& out) const;
    virtual void getTimeMBR(TimeRegion& out) const;
    TimePoint getPointAtTime(double time) const;
    static std::vector<SpatialIndex::TimePoint> simplifyWithRDP(std::vector<SpatialIndex::TimePoint>& Points, double threshold);
    std::vector<Trajectory> cuttraj(std::vector<SpatialIndex::TimePoint>);
    std::vector<Trajectory> getSegments(double threshold);
    void linkTrajectory(Trajectory other);
    virtual void getPartialTrajectory(double tstart,double tend,Trajectory &out) const;

    void loadFromString(std::string s);

    static double line2lineDistance(const TimePoint &p1s,const TimePoint &p1e,const TimePoint &p2s,const TimePoint &p2e);
    static double line2MBRDistance(const TimePoint &ps,const TimePoint &pe,const Region &r);
    static double line2MBCDistance(const TimePoint &ps,const TimePoint &pe,const MBC &r);

    std::vector<TimePoint> m_points;

    friend SIDX_DLL std::ostream operator<<(std::ostream os,const Trajectory &r);

    static const uint32_t m_dimension=2;
    virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<Trajectory> TrajectoryPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Trajectory& r);
}

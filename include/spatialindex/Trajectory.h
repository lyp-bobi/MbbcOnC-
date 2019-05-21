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
    virtual double getMinimumDistance(const Region& in) const;
    virtual double getMinimumDistance(const Mbbc& in) const;
    virtual double getMinimumDistance(const MBRk& in) const;
    virtual double getMinimumDistance(const MBBCk& in) const;
    virtual double getMinimumDistance(const Trajectory& in) const;

    virtual bool intersectsTimeRegion(const TimeRegion& in) const;
    virtual bool intersectsRegion(const Region& in) const;
    virtual bool intersectsMbbc(const Mbbc& in) const;
    virtual bool intersectsMBRk(const MBRk& in) const;
    virtual bool intersectsTrajectory(const Trajectory& in) const;

    virtual void combineTrajectory(const Trajectory& in);
    virtual bool containsTrajectory(const Trajectory& in);
    virtual void getCombinedTrajectory(Trajectory& out, const Trajectory& in) const;

    const std::string toString() const;
    virtual void getMbbc(Mbbc& out,bool tight) const;
    TimePoint getPointAtTime(double time) const;
    virtual void getMBRk(int k,MBRk &out) const;
    static std::vector<SpatialIndex::TimePoint> simplifyWithRDP(std::vector<SpatialIndex::TimePoint>& Points, double epsilon);
    virtual void getMBBCk(int k,MBBCk &out,double eps) const;

    void loadFromString(std::string s);

    std::vector<TimePoint> m_points;

    friend SIDX_DLL std::ostream operator<<(std::ostream os,const Trajectory &r);

    static const uint32_t m_dimension=2;
    virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<Trajectory> TrajectoryPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Trajectory& r);
}

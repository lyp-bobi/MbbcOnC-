//
// Created by Chuang on 2019/10/17.
//

#ifndef SPATIALINDEX_CYLINDER_H
#define SPATIALINDEX_CYLINDER_H
//
// Created by chuang on 5/29/19.
//
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL Cylinder: public Tools::IObject, public virtual IShape,public IEvolvingShape{

    public:
    Cylinder();
    ~Cylinder();
    Cylinder(const double* p,double r,double sTime,double eTime, uint32_t dimension);
    Cylinder(const Cylinder& in);
    Cylinder &operator=(const Cylinder &r);

    virtual bool operator==(const Cylinder &r) const;


    //
    // IObject interface
    //
    virtual Cylinder *clone();

    //
    // ISerializable interface
    //
    virtual uint32_t getByteArraySize() const;

    virtual void loadFromByteArray(const uint8_t *data);

    virtual void storeToByteArray(uint8_t **data, uint32_t &len);

    //
    // IEvolvingShape interface
    //
    virtual void getVMBR(Region& out) const;
    virtual void getMBRAtTime(double t, Region& out) const;


    virtual std::pair<STPoint,double> getCenterRdAtTime(double t) const;


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
    virtual double getMinimumDistance(const Point& in) const;


    virtual bool intersectsTimeRegion(const TimeRegion& in) const;
    virtual bool intersectsSTPoint(const STPoint& in) const;
    virtual bool intersectsRegion(const Region& in) const;
    virtual bool intersectsCylinder(const Cylinder& in) const;


    int checkRel(const MBC &bc) const;
    int checkRel(const Region &br) const;

//        virtual void combineCylinder(const Cylinder& in);
//        virtual bool containsCylinder(const Cylinder& in);
//        virtual void getCombinedCylinder(Cylinder& out, const Cylinder& in) const;;


    double m_startTime;
    double m_endTime;
    double* m_p;
    double m_r;
    uint32_t m_dimension=2;

    friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const Cylinder& r);

    virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<Cylinder> CylinderPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Cylinder& r);
}

#endif //SPATIALINDEX_CYLINDER_H

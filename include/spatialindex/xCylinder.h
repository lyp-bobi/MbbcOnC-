//
// Created by Chuang on 2019/10/17.
//

#ifndef SPATIALINDEX_xCylinder_H
#define SPATIALINDEX_xCylinder_H
//
// Created by chuang on 5/29/19.
//
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL xCylinder: public Tools::IObject, public virtual IxShape{

    public:
    xCylinder();
    xCylinder(const xPoint &p,double r,double sTime,double eTime);
    xCylinder(const xPoint &p,double r,double sTime,double eTime, uint32_t dimension);
    xCylinder(const xCylinder& in);
    xCylinder &operator=(const xCylinder &r);

    //
    // IObject interface
    //
    virtual xCylinder *clone();

    //
    // ISerializable interface
    //
    virtual uint32_t getByteArraySize() const;

    virtual void loadFromByteArray(const uint8_t *data);

    virtual void storeToByteArray(uint8_t **data, uint32_t &len);

    virtual void storeToByteArrayE(uint8_t** data, uint32_t& len);
    //
    // IShape interface
    //
    virtual bool intersectsShape(const IShape& in) const;
    virtual bool containsShape(const IShape& in) const;
    virtual bool touchesShape(const IShape& in) const;
    virtual void getCenter(Point& out) const;
    virtual uint32_t getDimension() const;
    virtual void getxMBR(xMBR& out) const;
    virtual double getArea() const;
    virtual double getMinimumDistance(const IShape& in) const;

    virtual double getMinimumDistance(const xMBR& in) const;
    virtual double getMinimumDistance(const xPoint& in) const;


    virtual bool intersectsxPoint(const xPoint& in) const;
    virtual bool intersectsxMBR(const xMBR& in) const;
    virtual bool intersectsxCylinder(const xCylinder& in) const;

    int checkRel(const xLine &bc) const;
    int checkRel(const xMBC &bc) const;
    int checkRel(const xMBR &br) const;
    int checkRel(const xSBB &b) const;

//        virtual void combinexCylinder(const xCylinder& in);
//        virtual bool containsxCylinder(const xCylinder& in);
//        virtual void getCombinedxCylinder(xCylinder& out, const xCylinder& in) const;;


    prex m_startTime;
    prex m_endTime;
    xPoint m_p;
    prex m_r;

    friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const xCylinder& r);

    virtual void makeInfinite(uint32_t dimension);
    private:
};
typedef Tools::PoolPointer<xCylinder> xCylinderPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const xCylinder& r);
}

#endif //SPATIALINDEX_xCylinder_H

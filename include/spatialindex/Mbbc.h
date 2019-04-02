//
// Created by chuang on 4/1/19.
//

#ifndef SPATIALINDEX_MBBC_H
#define SPATIALINDEX_MBBC_H
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL Mbbc: public Tools::IObject, public virtual IShape,public IEvolvingShape{

    public:
        Mbbc(const Region &smbr, const Region &embr, const Region &vmbr, const Region &pmbr, double tStart, double tEnd);
        Mbbc(const Mbbc& in);
        virtual Mbbc &operator=(const Mbbc &r);

        virtual bool operator==(const Mbbc &r) const;


        //
        // IObject interface
        //
        virtual Mbbc *clone();

        //
        // ISerializable interface
        //
        virtual uint32_t getByteArraySize();

        virtual void loadFromByteArray(const byte *data);

        virtual void storeToByteArray(byte **data, uint32_t &len);

        //
        // IEvolvingShape interface
        //
        virtual void getVMBR(Region& out) const;
        virtual void getMBRAtTime(double t, Region& out) const;


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


        virtual bool intersectsTimeRegion(const TimeRegion& in) const;
        virtual bool intersectsRegion(const Region& in) const;
        virtual bool intersectsLineSegment(const LineSegment& in) const;
        virtual bool containsPoint(const Point& in) const;
        virtual bool intersectsMbbc(const Mbbc& in) const;
        Region m_smbr;
        Region m_embr;
        Region m_vmbr;
        Region m_pmbr;
        double m_startTime;
        double m_endTime;

        friend SIDX_DLL std::ostream operator<<(std::ostream os,const MovingRegion &r);

    private:

    };
typedef Tools::PoolPointer<Mbbc> MbbcPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Mbbc& r);
}
#endif //SPATIALINDEX_MBBC_H

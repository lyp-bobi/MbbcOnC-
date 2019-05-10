//
// Created by chuang on 4/1/19.
//

#pragma once

namespace SpatialIndex
{
    class SIDX_DLL Mbbc: public Tools::IObject, public virtual IShape,public IEvolvingShape{

    public:
        Mbbc();
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

        virtual void loadFromByteArray(const uint8_t *data);

        virtual void storeToByteArray(uint8_t **data, uint32_t &len);

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
        virtual double getMinimumDistance(const Region& in) const;
        virtual double getMinimumDistance(const TimePoint& in) const;


        virtual bool intersectsTimeRegion(const TimeRegion& in) const;
        virtual bool intersectsTimePoint(const TimePoint& in) const;
        virtual bool intersectsRegion(const Region& in) const;
        virtual bool intersectsLineSegment(const LineSegment& in) const;
        virtual bool containsPoint(const Point& in) const;
        virtual bool intersectsMbbc(const Mbbc& in) const;

        virtual void combineMbbc(const Mbbc& in);
        virtual bool containsMbbc(const Mbbc& in);
        virtual void getCombinedMbbc(Mbbc& out, const Mbbc& in) const;

        const std::string toString() const;


        Region m_smbr;
        Region m_embr;
        Region m_xmbr;
        Region m_ymbr;
        Region m_vmbr;
        Region m_wmbr;
        double m_startTime;
        double m_endTime;

        friend SIDX_DLL std::ostream operator<<(std::ostream os,const MovingRegion &r);

        static const uint32_t m_dimension=2;
        virtual void makeInfinite();
    private:
    };
typedef Tools::PoolPointer<Mbbc> MbbcPtr;
SIDX_DLL std::ostream& operator<<(std::ostream& os, const Mbbc& r);
}

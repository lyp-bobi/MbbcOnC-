//
// Created by Chuang on 2019/5/15.
//


#pragma once

namespace SpatialIndex
{
    class SIDX_DLL MBRk: public Tools::IObject, public virtual IShape,public IEvolvingShape{

        public:
        inline int getPhase(double t) const;
        MBRk();
        MBRk(const std::vector<Region> mbrs, double tStart, double tEnd);
        MBRk(const MBRk& in);
        virtual MBRk &operator=(const MBRk &r);

        virtual bool operator==(const MBRk &r) const;


        //
        // IObject interface
        //
        virtual MBRk *clone();

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
        virtual bool intersectsMBRk(const MBRk& in) const;

        virtual void combineMBRk(const MBRk& in);
        virtual bool containsMBRk(const MBRk& in);
        virtual void getCombinedMBRk(MBRk& out, const MBRk& in) const;


        int m_k;
        std::vector<Region> m_mbrs;
        double m_startTime;
        double m_endTime;

        friend SIDX_DLL std::ostream operator<<(std::ostream os,const MBRk &r);

        static const uint32_t m_dimension=2;
        virtual void makeInfinite();
        private:
    };
    typedef Tools::PoolPointer<MBRk> MBRkPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const MBRk& r);
}


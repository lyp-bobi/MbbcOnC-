//
// Created by Chuang on 2019/5/17.
//


#pragma once

namespace SpatialIndex
{
    class SIDX_DLL MBBCk: public Tools::IObject, public virtual IShape,public IEvolvingShape{

        public:
        int getPhase(double t) const;
        MBBCk();
        MBBCk(int k);
        MBBCk(const std::vector<Region> mbrs,const std::vector<Region> vmbrs,const std::vector<Region> wmbrs, double tStart, double tEnd);
        MBBCk(const MBBCk& in);
        virtual MBBCk &operator=(const MBBCk &r);

        virtual bool operator==(const MBBCk &r) const;


        //
        // IObject interface
        //
        virtual MBBCk *clone();

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
        virtual bool intersectsMBBCk(const MBBCk& in) const;

        virtual void combineMBBCk(const MBBCk& in);
        virtual bool containsMBBCk(const MBBCk& in);
        virtual void getCombinedMBBCk(MBBCk& out, const MBBCk& in) const;




        uint32_t m_k;
        std::vector<Region> m_mbrs;
        std::vector<Region> m_vmbrs;
        std::vector<Region> m_wmbrs;
        double m_startTime;
        double m_endTime;

        friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const MBBCk& r);

        static const uint32_t m_dimension=2;
        virtual void makeInfinite(uint32_t dimension,int k);
        private:
    };
    typedef Tools::PoolPointer<MBBCk> MBBCkPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const MBBCk& r);
}


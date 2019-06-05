//
// Created by chuang on 5/29/19.
//
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL MBC: public Tools::IObject, public virtual IShape,public IEvolvingShape{

    public:
        MBC();
        ~MBC();
        MBC(const double* pLow, const double* pHigh,double sTime,double eTime, uint32_t dimension,double rd,double rv);
        MBC(const MBC& in);
        MBC &operator=(const MBC &r);

        virtual bool operator==(const MBC &r) const;


        //
        // IObject interface
        //
        virtual MBC *clone();

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
        virtual bool intersectsMBC(const MBC& in) const;

        virtual void combineMBC(const MBC& in);
        virtual bool containsMBC(const MBC& in);
        virtual void getCombinedMBC(MBC& out, const MBC& in) const;
        virtual int getOrient() const;

        virtual void getTimeMBR(TimeRegion& out) const;

        double m_startTime;
        double m_endTime;
        double m_rd;
        double m_rv;
        double* m_pLow;
        double* m_pHigh;
        uint32_t m_dimension=2;

        friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const MBC& r);

        virtual void makeInfinite(uint32_t dimension);
    private:
    };
    class SIDX_DLL MBCs: public Tools::IObject, public virtual IShape,public IEvolvingShape{
    public:
        //
        // IObject interface
        //
        virtual MBCs *clone();

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

        virtual void getTimeMBR(TimeRegion& out) const;

        MBCs()= default;
        MBCs(const MBCs &in);
        uint32_t m_dimension=2;
        std::vector<id_type> m_ids;
        std::vector<MBC> m_mbcs;
    };
    typedef Tools::PoolPointer<MBC> MBCPtr;
    typedef Tools::PoolPointer<MBCs> MBCsPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const MBC& r);
}
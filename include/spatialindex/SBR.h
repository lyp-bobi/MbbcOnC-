//
// Created by chuang on 5/29/19.
//
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL SBR: public Tools::IObject, public virtual IShape,public IEvolvingShape{

    public:
        SBR();
        ~SBR();
        SBR(const double* pLow, const double* pHigh,double sTime,double eTime, uint32_t dimension,double rd);
        SBR(const SBR& in);
        SBR &operator=(const SBR &r);

        virtual bool operator==(const SBR &r) const;


        //
        // IObject interface
        //
        virtual SBR *clone();

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
        virtual bool touchesSBR(const SBR& in) const;
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

        virtual bool containsSBR(const SBR& in);
        virtual int getOrient() const;

        virtual double getMargin() const;
        virtual void combineMBC(const MBC& in);
        virtual void combineSBR(const SBR& in);
        virtual void getFromMBC(const MBC& in,double tstart,double tend);
        virtual void getCombinedSBR(SBR& out, const SBR& in) const;
        static SBR getSBR(std::vector<SBR>& in);
        virtual double getIntersectingArea(const SBR& in) const;

        double m_startTime;
        double m_endTime;
        double m_rd;
        double* m_pLow;
        double* m_pHigh;
        uint32_t m_dimension=2;

        friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const SBR& r);

        virtual void makeInfinite(uint32_t dimension);
    private:
    };
    typedef Tools::PoolPointer<SBR> SBRPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const SBR& r);
}

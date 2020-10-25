//
// Created by chuang on 5/29/19.
//
#pragma once

extern bool bCompactMBC;

namespace SpatialIndex
{
    class SIDX_DLL xMBC: public Tools::IObject, public virtual IxShape{

    public:
        xMBC();
        xMBC(const xMBC& in);
        xMBC(xPoint ps,xPoint pe, prex rd, prex rv);
        xMBC &operator=(const xMBC &r);

        virtual bool operator==(const xMBC &r) const;


        //
        // IObject interface
        //
        virtual xMBC *clone();

        //
        // ISerializable interface
        //
        virtual uint32_t getByteArraySize() const;

        virtual void loadFromByteArray(const uint8_t *data);

        virtual void storeToByteArray(uint8_t **data, uint32_t &len);

        //
        // IEvolvingShape interface
        //
        virtual void getVMBR(xMBR& out) const;
        virtual void getMBRAtTime(double t, xMBR& out) const;


        virtual std::pair<xPoint,double> getCenterRdAtTime(double t) const;


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
        virtual bool intersectsxMBC(const xMBC& in) const;

//        virtual void combinexMBC(const xMBC& in);
//        virtual bool containsxMBC(const xMBC& in);
//        virtual void getCombinedxMBC(xMBC& out, const xMBC& in) const;;
        virtual bool prevalidate(const xMBR& in) const;

        prex m_rd;
        prex m_rv;
        xPoint m_ps,m_pe;


        friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const xMBC& r);

        virtual void makeInfinite(uint32_t dimension);
    private:
    };
    typedef Tools::PoolPointer<xMBC> xMBCPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const xMBC& r);
}

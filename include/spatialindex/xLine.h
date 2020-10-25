//
// Created by chuang on 5/29/19.
//
#pragma once

namespace SpatialIndex
{
    class SIDX_DLL xLine: public Tools::IObject, public virtual IxShape{

    public:
        xLine();
        xLine(const xLine& in);
        xLine(xPoint ps,xPoint pe, prex rd, prex rv);
        xLine &operator=(const xLine &r);

        virtual bool operator==(const xLine &r) const;


        //
        // IObject interface
        //
        virtual xLine *clone();

        //
        // ISerializable interface
        //
        virtual uint32_t getByteArraySize() const;

        virtual void loadFromByteArray(const uint8_t *data);

        virtual void storeToByteArray(uint8_t **data, uint32_t &len);

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
        virtual bool intersectsxLine(const xLine& in) const;



        xPoint m_ps,m_pe;


        friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const xLine& r);

        virtual void makeInfinite(uint32_t dimension);
    private:
    };
    typedef Tools::PoolPointer<xLine> xLinePtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const xLine& r);
}

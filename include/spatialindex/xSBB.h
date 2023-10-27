//
// Created by chuang on 5/29/19.
//
#pragma once


namespace SpatialIndex
{
    class SIDX_DLL xSBB{

    public:
        xSBB();
        xSBB(const xSBB& in);
        xSBB(xMBR &r);
        xSBB(xMBC &r);
        xSBB(xLine &r);
        xSBB(xMBR &r,xMBC &r2);
        xSBB(xPoint &r1, xPoint &r2);
        ~xSBB();
        xSBB &operator=(const xSBB &r);

        double m_startTime=1e300;
        double m_endTime=1e300;
        double tdist(const xPoint &p) const;

        std::pair<xPoint,double> crossSec(double t);

        virtual bool operator==(const xSBB &r) const;
        virtual xSBB *clone();

        bool hasbr=false;
        bool hasbc=false;
        bool hasbl=false;

        xMBR br;
        xMBC bc;
        xLine bl;

        std::string toString() const ;
        void loadFromString(std::string s);

        void loadFromByteArray(int type, const uint8_t* ptr);

        void init();

    private:
    };
    typedef Tools::PoolPointer<xSBB> xSBBPtr;
    SIDX_DLL std::ostream& operator<<(std::ostream& os, const xSBB& r);
}


#pragma once

namespace SpatialIndex
{
	class SIDX_DLL xPoint : public Tools::IObject, public virtual IxShape
	{
	public:
		xPoint(){};
		xPoint(prex x, prex y, prex t);
		xPoint(const xPoint& p);
		virtual ~xPoint();

		virtual xPoint& operator=(const xPoint& p);
		virtual bool operator==(const xPoint& p) const;

		//
		// IObject interface
		//
		virtual xPoint* clone();

		//
		// ISerializable interface
		//
		virtual uint32_t getByteArraySize() const;
		virtual void loadFromByteArray(const uint8_t* data);
		virtual void storeToByteArray(uint8_t** data, uint32_t& len);

        virtual void storeToByteArrayE(uint8_t** data, uint32_t& len);

		virtual void makeInfinite(uint32_t dimension);
		virtual void makeDimension(uint32_t dimension);

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

        virtual double getMinimumDistance(const xPoint& p) const;



	public:
        prex m_x,m_y;
        prex m_t;
        inline prex m_pCoords(int i) const{
            if(i==0) return m_x;
            if(i==1) return m_y;
            if(i==2) return m_t;
            return NAN;
        }
        static xPoint makemid(const xPoint &tp1,const xPoint &tp2,double t);
		friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const xPoint& pt);
	}; // xPoint
    typedef Tools::PoolPointer<xPoint> xPointPtr;
	SIDX_DLL std::ostream& operator<<(std::ostream& os, const xPoint& pt);
}

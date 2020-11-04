/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2004, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#pragma once

namespace SpatialIndex
{
	class SIDX_DLL xMBR : public Tools::IObject, public virtual IxShape
	{
	public:
		xMBR(){}
		xMBR(const xMBR& p);
		virtual ~xMBR(){}


		virtual xMBR& operator=(const xMBR& r);
		virtual bool operator==(const xMBR&) const;

		//
		// IObject interface
		//
		virtual xMBR* clone();

		//
		// ISerializable interface
		//
		virtual uint32_t getByteArraySize() const;
		virtual void loadFromByteArray(const uint8_t* data);
		virtual void storeToByteArray(uint8_t** data, uint32_t& length);
        virtual void storeToByteArrayE(uint8_t** data, uint32_t& len);

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


		virtual void getCenter(xPoint& out) const;
		virtual bool intersectsxMBR(const xMBR& in) const;
		virtual bool containsxMBR(const xMBR& in) const;
		virtual bool touchesxMBR(const xMBR& in) const;
		virtual double getMinimumDistance(const xMBR& in) const;


//        virtual bool intersectsMBC(const MBC& in) const;

		virtual bool containsxPoint(const xPoint& in) const;

		virtual double getMinimumDistance(const xPoint& in) const;

		virtual double getIntersectingArea(const xMBR& in) const;
		virtual double getMargin() const;

		virtual void combinexMBR(const xMBR& in);
		virtual void combinexPoint(const xPoint& in);
		virtual void getCombinedxMBR(xMBR& out, const xMBR& in) const;

		virtual void makeInfinite(uint32_t dimension);
		virtual void makeDimension(uint32_t dimension);

	private:

	public:
        prex m_xmin,m_xmax,m_ymin,m_ymax;
        prex m_tmin, m_tmax;

        prex& m_pLow(int i){
            if(i==0) return m_xmin;
            if(i==1) return m_ymin;
            if(i==2) return m_tmin;
            throw Tools::IllegalStateException("error plow index");
        }
        prex& m_pHigh(int i){
            if(i==0) return m_xmax;
            if(i==1) return m_ymax;
            if(i==2) return m_tmax;
            throw Tools::IllegalStateException("error phigh index");
        }


		friend SIDX_DLL std::ostream& operator<<(std::ostream& os, const xMBR& r);
	}; // xMBR
	
	typedef Tools::PoolPointer<xMBR> xMBRPtr;
	SIDX_DLL std::ostream& operator<<(std::ostream& os, const xMBR& r);
}

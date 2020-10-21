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

#include <spatialindex/SpatialIndex.h>

#include <cstring>
#include <string>
#include <cmath>
#include <limits>

using namespace SpatialIndex;


xMBR::xMBR(const xMBR &r)
:m_xmin(r.m_xmin),m_xmax(r.m_xmax),m_ymin(r.m_ymin),
    m_ymax(r.m_ymax),m_tmin(r.m_tmin),m_tmax(r.m_tmax)
{
}

xMBR& xMBR::operator=(const xMBR& r)
{
	m_xmin = r.m_xmin;
	m_xmax = r.m_xmax;
	m_ymin = r.m_ymin;
	m_ymax = r.m_ymax;
	m_tmin = r.m_tmin;
	m_tmax = r.m_tmax;
	return *this;
}

bool xMBR::operator==(const xMBR& r) const
{
    if (
        m_xmin < r.m_xmin - std::numeric_limits<double>::epsilon() ||
        m_xmin > r.m_xmin + std::numeric_limits<double>::epsilon() ||
        m_xmax < r.m_xmax - std::numeric_limits<double>::epsilon() ||
        m_xmax > r.m_xmax + std::numeric_limits<double>::epsilon() ||
        m_ymin < r.m_ymin - std::numeric_limits<double>::epsilon() ||
        m_ymin > r.m_ymin + std::numeric_limits<double>::epsilon() ||
        m_ymax < r.m_ymax - std::numeric_limits<double>::epsilon() ||
        m_ymax > r.m_ymax + std::numeric_limits<double>::epsilon() ||
        m_tmin < r.m_tmin - std::numeric_limits<double>::epsilon() ||
        m_tmin > r.m_tmin + std::numeric_limits<double>::epsilon() ||
        m_tmax < r.m_tmax - std::numeric_limits<double>::epsilon() ||
        m_tmax > r.m_tmax + std::numeric_limits<double>::epsilon()
                )
        return false;
	return true;
}

//
// IObject interface
//
xMBR* xMBR::clone()
{
	return new xMBR(*this);
}

//
// ISerializable interface
//
uint32_t xMBR::getByteArraySize() const
{
	return 6*sizeof(prex);
}

void xMBR::loadFromByteArray(const uint8_t* ptr)
{
    memcpy(&m_xmin, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_xmax, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_ymin, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_ymax, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_tmin, ptr, sizeof(prex));
    ptr += sizeof(prex);
    memcpy(&m_tmax, ptr, sizeof(prex));
    ptr += sizeof(prex);
}

void xMBR::storeToByteArray(uint8_t** data, uint32_t& len)
{
	len = getByteArraySize();
	*data = new uint8_t[len];
	uint8_t* ptr = *data;

    memcpy(ptr, &m_xmin, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_xmax, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_ymin, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_ymax, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_tmin, sizeof(double));
    ptr += sizeof(double);
    memcpy(ptr, &m_tmax, sizeof(double));
    ptr += sizeof(double);
    
}
//
// IShape interface
//
bool xMBR::intersectsShape(const IShape& s) const
{
    const xMBC* pbc = dynamic_cast<const xMBC*>(&s);
    if (pbc != 0) return pbc->intersectsxMBR(*this);

    const xCylinder* pcy = dynamic_cast<const xCylinder*>(&s);
    if (pcy != 0) return pcy->intersectsxMBR(*this);

	const xMBR* pr = dynamic_cast<const xMBR*>(&s);
	if (pr != 0) return intersectsxMBR(*pr);

	return s.intersectsShape(*this);

	throw Tools::IllegalStateException(
		"xMBR::intersectsShape: Not implemented yet!"
	);
}

bool xMBR::containsShape(const IShape& s) const
{
	const xMBR* pr = dynamic_cast<const xMBR*>(&s);
	if (pr != 0) return containsxMBR(*pr);

	const xPoint* ppt = dynamic_cast<const xPoint*>(&s);
	if (ppt != 0) return containsxPoint(*ppt);

	throw Tools::IllegalStateException(
		"xMBR::containsShape: Not implemented yet!"
	);
}

bool xMBR::touchesShape(const IShape& s) const
{
	const xMBR* pr = dynamic_cast<const xMBR*>(&s);
	if (pr != 0) return touchesxMBR(*pr);

	throw Tools::IllegalStateException(
		"xMBR::touchesShape: Not implemented yet!"
	);
}

void xMBR::getCenter(Point& out) const
{
}

uint32_t xMBR::getDimension() const
{
	return 2;
}

void xMBR::getMBR(Region& out) const
{
}

double xMBR::getArea() const
{

	return (m_xmax-m_xmin)*(m_ymax-m_ymin)*(m_tmax-m_tmin);
}

double xMBR::getMinimumDistance(const IShape& s) const
{
    const MBC* pbc = dynamic_cast<const MBC*>(&s);
    if (pbc != 0) return pbc->getMinimumDistance(*this);

	const xMBR* pr = dynamic_cast<const xMBR*>(&s);
	if (pr != 0) return getMinimumDistance(*pr);

	const xPoint* ppt = dynamic_cast<const xPoint*>(&s);
	if (ppt != 0) return getMinimumDistance(*ppt);

	throw Tools::IllegalStateException(
		"xMBR::getMinimumDistance: Not implemented yet!"
	);
}

bool xMBR::intersectsxMBR(const xMBR& r) const
{
    if(m_xmin>r.m_xmax || m_xmax<r.m_xmin) return false;
    if(m_ymin>r.m_ymax || m_ymax<r.m_ymin) return false;
    if(m_tmin>r.m_tmax || m_tmax<r.m_tmin) return false;
	return true;
}

bool xMBR::containsxMBR(const xMBR& r) const
{
    if(m_xmin>r.m_xmin || m_xmax<r.m_xmax) return false;
    if(m_ymin>r.m_ymin || m_ymax<r.m_ymax) return false;
    if(m_tmin>r.m_tmin || m_tmax<r.m_tmax) return false;
	return true;
}

bool xMBR::touchesxMBR(const xMBR& r) const
{
    return false;
}


double xMBR::getMinimumDistance(const xMBR& r) const
{
	double ret = 0.0;
	
	double x = 0.0;
    if (r.m_xmax < m_xmin){x = std::abs(r.m_xmax - m_xmin);}
    else if (m_xmax < r.m_xmin){x = std::abs(r.m_xmin - m_xmax);}
    ret += x * x;
    if (r.m_ymax < m_ymin){x = std::abs(r.m_ymax - m_ymin);}
    else if (m_ymax < r.m_ymin){x = std::abs(r.m_ymin - m_ymax);}
    ret += x * x;
    if (r.m_tmax < m_tmin){x = std::abs(r.m_tmax - m_tmin);}
    else if (m_tmax < r.m_tmin){x = std::abs(r.m_tmin - m_tmax);}
    ret += x * x;

	return std::sqrt(ret);
}


//bool xMBR::intersectsMBC(const SpatialIndex::MBC &in) const {
//    if (m_pLow[m_dimension] > in.m_endTime || m_pHigh[m_dimension] < in.m_startTime) return false;
//    if (in.m_pLow[m_dimension] == in.m_pHigh[m_dimension]) {
//        double x1=in.m_pLow[0],y1=in.m_pLow[1],t1=in.m_startTime;
//        double x2=in.m_pHigh[0],y2=in.m_pHigh[1],t2=in.m_endTime;
//        double xt=makemidmacro(x1,t1,x2,t2,m_pLow[m_dimension]),
//            yt=makemidmacro(y1,t1,y2,t2,m_pLow[m_dimension]);
//        double xy[2]={xt,yt};
//        double mbcr=std::min(std::min(in.m_rd,(m_pLow[m_dimension]-in.m_startTime)*in.m_rv),(in.m_endTime-m_pLow[m_dimension])*in.m_rv);
//        return xPoint(xy,2).getMinimumDistance(xMBR(in.m_pLow, in.m_pHigh, m_dimension)) <= mbcr + 1e-7;
//        xMBR br = xMBR(in.m_pLow, in.m_pHigh, 2);
//        return timed.first.getMinimumDistance(br) <= timed.second;
//    }
//}
bool xMBR::containsxPoint(const xPoint& p) const
{

    if (m_xmin > p.m_x|| m_xmax < p.m_x) return false;
    if (m_ymin > p.m_y|| m_ymax < p.m_y) return false;
    if (m_tmin > p.m_t|| m_tmax < p.m_t) return false;
	return true;
}


double xMBR::getMinimumDistance(const xPoint& p) const
{
	double ret = 0.0;
    if (p.m_x < m_xmin)
    {ret += std::pow(m_xmin - p.m_x, 2.0);}
    else if (p.m_x > m_xmax)
    {ret += std::pow(p.m_x - m_xmax, 2.0);}

    if (p.m_y < m_ymin)
    {ret += std::pow(m_ymin - p.m_y, 2.0);}
    else if (p.m_y > m_ymax)
    {ret += std::pow(p.m_y - m_ymax, 2.0);}

	return std::sqrt(ret);
}

double xMBR::getIntersectingArea(const xMBR& r) const
{
	double ret = 1.0;
	double f1, f2;


    if (m_xmin > r.m_xmax || r.m_xmax < r.m_xmin) return 0.0;
    f1 = std::max(m_xmin, r.m_xmin);
    f2 = std::min(r.m_xmax, r.m_xmax);
    ret *= f2 - f1;
    if (m_ymin > r.m_ymax || r.m_ymax < r.m_ymin) return 0.0;
    f1 = std::max(m_ymin, r.m_ymin);
    f2 = std::min(r.m_ymax, r.m_ymax);
    ret *= f2 - f1;
    if (m_tmin > r.m_tmax || r.m_tmax < r.m_tmin) return 0.0;
    f1 = std::max(m_tmin, r.m_tmin);
    f2 = std::min(r.m_tmax, r.m_tmax);
    ret *= f2 - f1;

	return ret;
}

/*
 * Returns the margin of a xMBR. It is calcuated as the sum of  2^(d-1) * width, in each dimension.
 * It is actually the sum of all edges, no matter what the dimensionality is.
*/
double xMBR::getMargin() const
{
	double mul = 4;
	double margin = 0.0;
	
	margin += (m_xmax - m_xmin) * mul;
    margin += (m_ymax - m_ymin) * mul;
    margin += (m_tmax - m_tmin) * mul;

	return margin;
}

void xMBR::combinexMBR(const xMBR& r)
{
    m_xmin = std::min(m_xmin, r.m_xmin);
    m_xmax = std::max(m_xmax, r.m_xmax);
    m_ymin = std::min(m_ymin, r.m_ymin);
    m_ymax = std::max(m_ymax, r.m_ymax);
    m_tmin = std::min(m_tmin, r.m_tmin);
    m_tmax = std::max(m_tmax, r.m_tmax);
}

void xMBR::combinexPoint(const xPoint& p)
{
    m_xmin = std::min(m_xmin, p.m_x);
    m_xmax = std::max(m_xmax, p.m_x);
    m_ymin = std::min(m_ymin, p.m_y);
    m_ymax = std::max(m_ymax, p.m_y);
    m_tmin = std::min(m_tmin, p.m_t);
    m_tmax = std::max(m_tmax, p.m_t);
}

void xMBR::getCombinedxMBR(xMBR& out, const xMBR& in) const
{
	out = *this;
	out.combinexMBR(in);
}


void xMBR::makeInfinite(uint32_t dimension)
{
	m_xmin=m_ymax=m_tmin= 1e30;
	m_xmax=m_ymax=m_tmax = -1e30;
}

void xMBR::makeDimension(uint32_t dimension)
{
}

std::ostream& SpatialIndex::operator<<(std::ostream& os, const xMBR& r)
{
	uint32_t i;

	os<<r.m_xmin<<","<<r.m_xmax<<","<<r.m_ymin<<","<<r.m_ymax<<","<<r.m_tmin<<","<<r.m_tmax;

	return os;
}

/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2003, Marios Hadjieleftheriou
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
#ifndef SPATIALINDEX_H
#define SPATIALINDEX_H

//#define TJDEBUG
#include "tools/Tools.h"
#include "tools/MutablePriorityQueue.h"
#include "tools/DiskMultiMap.h"
#include <chrono>
#include "tools/json.hpp"
#include <cmath>

using json = nlohmann::json;

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define sq(x) (x)*(x)
#define cube(x) (x)*(x)*(x)
#define makemidmacro(x1,t1,x2,t2,t) ((t)-(t1))/((t2)-(t1))*(x2)+((t2)-(t))/((t2)-(t1))*(x1)

#define prex double  //index precision

#define prexp double //storage precision

namespace SpatialIndex
{
	class Point;
	class LineSegment;
	class Region;
	class xMBR;
    class xTrajectory;


	struct DISTE{
	    double opt=0;
	    double pes=0;
	    bool infer = false;
	    DISTE(){};
	    explicit DISTE(double optimistic,double pessimistic, bool isinferred)
	    {
#ifndef NDEBUG
	        assert(optimistic<=pessimistic);
#endif
	        opt=optimistic;pes=pessimistic;infer=isinferred;
	    }
	    explicit DISTE(double exact){
	        opt=pes=exact;
	        infer = false;
	    }
	    DISTE operator+(const DISTE &d2) const{
	        DISTE res;
            res.opt=opt+d2.opt;
            res.pes=pes+d2.pes;
            res.infer = infer||d2.infer;
            return res;
	    }

        bool operator<(const DISTE &d2) const{
            return opt<d2.opt;
        }
	};

	typedef int64_t id_type;

	SIDX_DLL enum CommandType
	{
		CT_NODEREAD = 0x0,
		CT_NODEDELETE,
		CT_NODEWRITE
	};

	class SIDX_DLL InvalidPageException : public Tools::Exception
	{
	public:
		InvalidPageException(id_type id);
		virtual ~InvalidPageException() {}
		virtual std::string what();

	private:
		std::string m_error;
	}; // InvalidPageException

	//
	// Interfaces
	//

	class SIDX_DLL IShape : public Tools::ISerializable
	{
	public:
		virtual bool intersectsShape(const IShape& in) const = 0;
		virtual bool containsShape(const IShape& in) const = 0;
		virtual bool touchesShape(const IShape& in) const = 0;
		virtual void getCenter(Point& out) const = 0;
		virtual uint32_t getDimension() const = 0;
		virtual void getMBR(Region& out) const = 0;
		virtual double getArea() const = 0;
		virtual double getMinimumDistance(const IShape& in) const = 0;
		virtual ~IShape() {}
	}; // IShape


    class SIDX_DLL IxShape : public IShape
    {
    public:
        void getMBR(Region& out) const {throw Tools::NotSupportedException("Use getxMBR instead!");}
        virtual void getxMBR(xMBR& out) const = 0;
    }; // IShape

	class SIDX_DLL ITimeShape : public Tools::IInterval
	{
	public:
		virtual bool intersectsShapeInTime(const ITimeShape& in) const = 0;
		virtual bool intersectsShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const = 0;
		virtual bool containsShapeInTime(const ITimeShape& in) const = 0;
		virtual bool containsShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const = 0;
		virtual bool touchesShapeInTime(const ITimeShape& in) const = 0;
		virtual bool touchesShapeInTime(const Tools::IInterval& ivI, const ITimeShape& in) const = 0;
		virtual double getAreaInTime() const = 0;
		virtual double getAreaInTime(const Tools::IInterval& ivI) const = 0;
		virtual double getIntersectingAreaInTime(const ITimeShape& r) const = 0;
		virtual double getIntersectingAreaInTime(const Tools::IInterval& ivI, const ITimeShape& r) const = 0;
		virtual ~ITimeShape() {}
	}; // ITimeShape

	class SIDX_DLL IEvolvingShape
	{
	public:
		virtual void getVMBR(Region& out) const = 0;
		virtual void getMBRAtTime(double t, Region& out) const = 0;
		virtual ~IEvolvingShape() {}
	}; // IEvolvingShape

	class SIDX_DLL IEntry : public Tools::IObject
	{
	public:
		virtual id_type getIdentifier() const = 0;
		virtual void getShape(IShape** out) const = 0;
		virtual ~IEntry() {}
	}; // IEntry

	class SIDX_DLL INode : public IEntry, public Tools::ISerializable
	{
	public:
		virtual uint32_t getChildrenCount() const = 0;
		virtual id_type getChildIdentifier(uint32_t index) const = 0;
		virtual void getChildData(uint32_t index, uint32_t& len, uint8_t** data) const = 0;
		virtual void getChildShape(uint32_t index, IShape** out) const = 0;
		virtual uint32_t getLevel() const = 0;
		virtual bool isIndex() const = 0;
		virtual bool isLeaf() const = 0;
		virtual ~INode() {}
	}; // INode

	class SIDX_DLL IData : public IEntry
	{
	public:
		virtual void getData(uint32_t& len, uint8_t** data) const = 0;
		virtual ~IData() {}
	}; // IData

	class SIDX_DLL IDataStream : public Tools::IObjectStream
	{
	public:
		virtual IData* getNext() = 0;
		virtual ~IDataStream() {}
	}; // IDataStream

	class SIDX_DLL ICommand
	{
	public:
		virtual void execute(const INode& in) = 0;
		virtual ~ICommand() {}
	}; // ICommand

	class SIDX_DLL INearestNeighborComparator
	{
	public:
		virtual double getMinimumDistance(const IShape& query, const IShape& entry) = 0;
		virtual double getMinimumDistance(const IShape& query, const IData& data) = 0;
		virtual ~INearestNeighborComparator() {}
	}; // INearestNeighborComparator

	class SIDX_DLL IStorageManager
	{
	public:
		virtual void loadByteArray(const id_type id, uint32_t& len, uint8_t** data) = 0;
		virtual void storeByteArray(id_type& id, const uint32_t len, const uint8_t* const data) = 0;
		virtual void deleteByteArray(const id_type id) = 0;
		virtual id_type nextPage(){return 0;}
		virtual void flush() = 0;
		virtual ~IStorageManager() {}
		bool m_isro = false;
	}; // IStorageManager

	class SIDX_DLL IVisitor
	{
	public:
		virtual void visitNode(const INode& in) = 0;
		virtual void visitData(const IData& in) = 0;
		virtual void visitData(std::vector<const IData*>& v) = 0;
		virtual ~IVisitor() {}
	}; // IVisitor

	class SIDX_DLL IQueryStrategy
	{
	public:
		virtual void getNextEntry(const IEntry& previouslyFetched, id_type& nextEntryToFetch, bool& bFetchNextEntry) = 0;
		virtual ~IQueryStrategy() {}
	}; // IQueryStrategy

	class SIDX_DLL IStatistics
	{
	public:
		virtual uint64_t getReads() const = 0;
		virtual uint64_t getWrites() const = 0;
		virtual uint32_t getNumberOfNodes() const = 0;
		virtual uint64_t getNumberOfData() const = 0;
		virtual ~IStatistics() {}
	}; // IStatistics

    SIDX_DLL enum DataType
    {
        BoundingBoxType= 0x0,
        TrajectoryType=0x1
    };
	class SIDX_DLL ISpatialIndex
	{
	public:
		virtual void insertData(uint32_t len, const uint8_t* pData, const IShape& shape, id_type shapeIdentifier) = 0;
		virtual bool deleteData(const IShape& shape, id_type shapeIdentifier) = 0;
		virtual void containsWhatQuery(const IShape& query, IVisitor& v)  = 0;
		virtual void intersectsWithQuery(const IShape& query, IVisitor& v) = 0;
		virtual void pointLocationQuery(const Point& query, IVisitor& v) = 0;
		virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v, INearestNeighborComparator& nnc) = 0;
		virtual void nearestNeighborQuery(uint32_t k, const IShape& query, IVisitor& v) = 0;
		virtual void selfJoinQuery(const IShape& s, IVisitor& v) = 0;
		virtual void queryStrategy(IQueryStrategy& qs) = 0;
		virtual void getIndexProperties(Tools::PropertySet& out) const = 0;
		virtual void addCommand(ICommand* in, CommandType ct) = 0;
		virtual bool isIndexValid() = 0;
		virtual void getStatistics(IStatistics** out) const = 0;
		virtual ~ISpatialIndex() {}

		DataType m_DataType =BoundingBoxType;

	}; // ISpatialIndex

	namespace StorageManager
	{
		SIDX_DLL enum StorageManagerConstants
		{
			EmptyPage = -0x1,
			NewPage = -0x1
		};

		class SIDX_DLL IBuffer : public IStorageManager
		{
		public:
			virtual uint64_t getHits() = 0;
			virtual void clear() = 0;
			virtual ~IBuffer() {}
		}; // IBuffer

		SIDX_DLL  IStorageManager* returnMemoryStorageManager(Tools::PropertySet& in);
		SIDX_DLL  IStorageManager* createNewMemoryStorageManager();

		SIDX_DLL  IStorageManager* returnDiskStorageManager(Tools::PropertySet& in);
		SIDX_DLL  IStorageManager* createNewDiskStorageManager(std::string& baseName, uint32_t pageSize);
		SIDX_DLL  IStorageManager* loadDiskStorageManager(std::string& baseName);

		SIDX_DLL  IBuffer* returnRandomEvictionsBuffer(IStorageManager& ind, Tools::PropertySet& in);
		SIDX_DLL  IBuffer* createNewRandomEvictionsBuffer(IStorageManager& in, uint32_t capacity, bool bWriteThrough);
	}

	//
	// Global functions
	//
	SIDX_DLL  std::ostream& operator<<(std::ostream&, const ISpatialIndex&);
	SIDX_DLL  std::ostream& operator<<(std::ostream&, const IStatistics&);


	class trajStat{
	private:
		trajStat(){};
		~trajStat(){delete singleton;}
		static trajStat* singleton;
	public:
		double bt=10000;
		double M=0;//total time
		long lineCount=0;
		long trajCount=0;
		double tl=0;//time len of a segment
		double jt=0; // averaged traj len
		double v=0; //averaged speed
		double vmax = 0; //maximum speed
		double minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,mint=1e300,maxt=-1e300;
		double Dx=0,Dy=0,Dt=0;
		double Sr=0; //radius
		double P=0; //temporal length
		double dist=0;
		double Df =2; //fractral dimension
		double f = 50,fp=170; //fanout, point capacity
		bool regular = false;
		string dataset = "";
		static trajStat* instance();
		void init(){
			bt=0;
			M=0;
			long lineCount=0;
			long trajCount=0;
			tl=0;//time len of a segment
			jt=0;
			v=0;
			minx=1e300,maxx=-1e300,miny=1e300,maxy=-1e300,mint=1e300,maxt=-1e300;
			Dx=0,Dy=0,Dt=0;
			Sr=0;
			P=0;
			dist=0;
			Df =2;
			f = 50;
			vmax=0;
		}
		void usedata(string str){
			dataset = str;
            if (dataset == "glexpand"){
                dataset="gl";
                M = 1.31808e+09;
                lineCount = 200692990;
                trajCount = 206190;
                tl = 6.56764;
                jt = 6392.55;
                Sr = 0.5;
                Df = 1.25;
                P = 60000;
                std::cerr<<"use glexpand sta\n";
            }
			else if (dataset == "gl"){
				M = 1.31808e+08;
				lineCount = 20069299;
				trajCount = 20619;
				tl = 6.56764;
				jt = 6392.55;
				Sr = 0.5;
				Df = 1.25;
				P = 60000;
				std::cerr<<"use gl sta\n";
			}
            else if(dataset == "tdexpand"){
                dataset="td";
                M = 5.20574e+10;
                lineCount = 162059560;
                trajCount = 102670;
                tl = 321.224;
                jt = 507037;
                Sr = 0.5;
                Df = 1.56;
                P = 5333150;
                std::cerr<<"use tdexpand sta\n";
            }
			else if(dataset == "td"){
				M = 5.20574e+09;
				lineCount = 16205956;
				trajCount = 10267;
				tl = 321.224;
				jt = 507037;
				Sr = 0.5;
				Df = 1.56;
				P = 533315;
				std::cerr<<"use td sta\n";
			}
			else if (dataset == "od"){
				Sr = 10000;
				Df = 1.4;
				P = 5000;
				std::cerr<<"use od sta\n";
			}
		}

		double rd(double x){
			if (dataset == "gl"){
				return 2e-05 * x + 0.0023;
			}else if(dataset == "td"){
				return 1E-05 * x + 0.0003;
			}
			else if (dataset == "od"){
				return -0.0326 * x*x + 20.289*x - 78.702;
			}
			return vv()*x;
		}

		double dd(double x){
			if (dataset == "gl"){
				return 4.5e-5*x - 1e-5/4000*x*x;
			}else if(dataset == "td"){
				return 4e-5*x;
			}
			else if (dataset == "od"){
				return 55*x;
			}
			return vv()*x;
		}

		double vv(){
			if (dataset == "gl"){
				return dd(bt)/bt;
			}else if(dataset == "td"){
				return dd(bt)/bt;
			}
			else if (dataset == "od"){
				return dd(bt)/bt;
			}
			return v;
		}

		double Lx(){
			if (dataset == "gl"){

			}else if(dataset == "td"){

			}
			else if (dataset == "od"){

			}else{
				return pow(M_PI*Sr*Sr*P*M*dd(bt)/bt/bt/f,0.33) + dd(bt);
			}
		}

		double Lt(){
			if (dataset == "gl"){

			}else if(dataset == "td"){

			}
			else if (dataset == "od"){

			}else{
				return (pow(M_PI*Sr*Sr*P*M*dd(bt)/bt/bt/f,0.33) + dd(bt))*bt/dd(bt);
			}
		}

		double knncost(double _bt,int k,double qt, bool out = false){
			bt = _bt;
			double leafs = M/bt/f;
			double nt = pow(leafs * sq(vv())*sq(P) / sq((M_PI * sq(Sr))),1.0/3);
			double ltc = P/nt;
			double lt = ltc+bt;
			double lxc = ltc * vv();
			double lx = lt * vv();
			double Nq = min(1.0,(qt+jt) /P)  * trajCount;
			double rk =  pow(k/Nq, Df/2)*Sr;
			double Rk = rk*qt/2 + qt/4* sqrt(2*sq(dd(qt)+4*sq(rk)));
			double Rnq = Rk +rd(bt)*qt;

			double nmin = k;
			double nmax = 100*k;
			while (nmax - nmin > 0.1) {
				double mid = (nmax + nmin) / 2;
				double rn =  pow(mid/Nq, Df/2)*Sr;
				double Rn = rn*qt/2 + qt/4* sqrt(2*sq(dd(qt)+4*sq(rn)));
				if (Rn > Rnq) {
					nmax = mid;
				} else {
					nmin = mid;
				}
			}
			double nq = (nmax+nmin)/2;

			double pruneCost = sq(lx + 2*Rk/qt + dd(qt))*(lt+qt)   *(M/bt/f/P/pow(M_PI*sq(Sr),Df/2));
			double filterCost = nq*(1+qt/fp/lt);
			if (out){
				std::cerr<<"nt "<<nt <<" lt "<<lt<<" lx "<<lx<<"\n";
				std::cerr<<"expect cost is "<<pruneCost<<"\t"<<filterCost<<"\n";
			}
			return pruneCost+filterCost/8;
		}
		void set(double _bt, double _M, long _lineCount, long _trajCount,
				 double _tl, double _jt, double _v, double _minx, double _maxx,
				 double _miny, double _maxy, double _mint, double _maxt,
				 double _Dx, double _Dy, double _Dt, double _dist)
		{
			bt=_bt; M=_M; lineCount=_lineCount; trajCount=_trajCount;
			tl=_tl; jt=_jt; v=_v; minx=_minx; maxx=_maxx;
			miny=_miny; maxy=_maxy; mint=_mint; maxt=_maxt;
			Dx=_Dx; Dy=_Dy; Dt=_Dt; dist=_dist;
		}
		void output(){
			std::cerr<<bt<<","<< M<<","<< lineCount<<","<< trajCount<<","<<
					 tl<<","<< jt<<","<< v<<","<< minx<<","<< maxx<<","<<
					 miny<<","<< maxy<<","<< mint<<","<< maxt<<","<<
					 Dx<<","<< Dy<<","<< Dt<<","<< dist<<"\n";
		}
		string toString(){
			ostringstream ostream;
			ostream<<bt<<" "<< M<<" "<< lineCount<<" "<< trajCount<<" "<<
				   tl<<" "<< jt<<" "<< v<<" "<< minx<<" "<< maxx<<" "<<
				   miny<<" "<< maxy<<" "<< mint<<" "<< maxt<<" "<<
				   Dx<<" "<< Dy<<" "<< Dt<<" "<< dist <<" "<<Sr<<" "<<P<<" "
				<<Df<<" "<<f<<" "<<vmax;
			return ostream.str();
		}
		void fromString(string s){
			istringstream istream(s);
			istream>>bt>> M>> lineCount>> trajCount>>
				   tl>> jt>> v>> minx>> maxx>>
				   miny>> maxy>> mint>> maxt>>
				   Dx>> Dy>> Dt>> dist>>Sr>>P>>Df>>f>>vmax;
			return;
		}
		friend SIDX_DLL std::ostream& operator<<(std::ostream& os,const trajStat &r);
	};
	SIDX_DLL std::ostream& operator<<(std::ostream& os, const trajStat& r);
}

extern SpatialIndex::trajStat *tjstat;

#include "Point.h"
#include "Region.h"
#include "LineSegment.h"
#include "TimePoint.h"
#include "STPoint.h"
#include "TimeRegion.h"
#include "MovingPoint.h"
#include "MovingRegion.h"
#include "MBC.h"
#include "Cylinder.h"
#include "Trajectory.h"
#include "ShapeList.h"


#include "xPoint.h"
#include "xMBR.h"
#include "xLine.h"
#include "xMBC.h"
#include "xSBB.h"
#include "xCylinder.h"
#include "xTrajectory.h"

#include "MBCRTree.h"
#include "xRTree.h"
#include "SBBForest.h"
#include "Version.h"


// typedef SpatialIndex::Tools::PoolPointer<Region> RegionPtr;
// typedef SpatialIndex::Tools::PoolPointer<Point> PointPtr;
// typedef SpatialIndex::Tools::PoolPointer<TimeRegion> TimeRegionPtr;
// typedef SpatialIndex::Tools::PoolPointer<MovingRegion> MovingRegionPtr;
#endif
//
// Created by Chuang on 2019/6/13.
//

#pragma once

namespace SpatialIndex {
    SIDX_DLL enum LeafBoundingType
    {
        LeafBoundByMBR = 0x1,
        LeafBoundByMBC = 0x2
    };
    class SIDX_DLL ShapeList: public Tools::IObject, public virtual IShape,public IEvolvingShape {
    public:

        //
        // IObject interface
        //
        virtual ShapeList *clone();

        //
        // ISerializable interface
        //
        virtual uint32_t getByteArraySize() const;

        virtual void loadFromByteArray(const uint8_t *data);

        virtual void storeToByteArray(uint8_t **data, uint32_t &len);

        //
        // IEvolvingShape interface
        //
        virtual void getVMBR(Region &out) const;

        virtual void getMBRAtTime(double t, Region &out) const;


        //
        // IShape interface
        //
        virtual bool intersectsShape(const IShape &in) const;

        virtual bool containsShape(const IShape &in) const;

        virtual bool touchesShape(const IShape &in) const;

        virtual void getCenter(Point &out) const;

        virtual uint32_t getDimension() const;

        virtual void getMBR(Region &out) const;

        virtual double getArea() const;

        virtual double getMinimumDistance(const IShape &in) const;

        virtual void getTimeMBR(TimeRegion &out) const;

        ShapeList() = default;

        ShapeList(const ShapeList &in);

        void insert(MBC* shape);
        void insert(Region* shape);

        void insert(IShape* shape);

        uint32_t m_dimension = 3;

        uint32_t m_datatype=01;

        std::vector <MBC*> m_MBCList;
        std::vector <Region*> m_MBRList;
    };
}
//
// Created by Chuang on 2019/6/13.
//

#pragma once

namespace SpatialIndex {
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

        uint32_t m_dimension = 3;
        std::vector <IShape*> m_ShapeList;
    };
}
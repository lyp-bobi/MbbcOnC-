#include <spatialindex/SpatialIndex.h>

using namespace SpatialIndex;
using namespace std;

/* 
 * Test the Geometry
 * Nowhere near complete, but it's something
 */
int main(int argc, char** argv) {
    //define points
    double c1[2] = {1.0, 0.0};
    double c2[2] = {3.0, 2.0};
    double c3[2] = {2.0, 0.0};
    double c4[2] = {2.0, 4.0};
    double c5[2] = {1.0, 1.0};
    double c6[2] = {2.5, 3.0};
    double c7[2] = {1.0, 2.0};
    double c8[2] = {0.0, -1.0};
    double c9[2] = {4.0, 3.0};
    Point p1 = Point(&c1[0], 2); 
    Point p2 = Point(&c2[0], 2); 
    Point p3 = Point(&c3[0], 2); 
    Point p4 = Point(&c4[0], 2); 
    Point p5 = Point(&c5[0], 2); 
    Point p6 = Point(&c6[0], 2); 
    Point p7 = Point(&c7[0], 2); 
    Point p8 = Point(&c8[0], 2); 
    Point p9 = Point(&c9[0], 2); 
    
    double c3a[2] = {2.0, 3.0};
    Point p3a = Point(&c3a[0], 2); 
    
    //Now Test LineSegment intersection
    LineSegment ls1 = LineSegment(p1, p2);
    LineSegment ls2 = LineSegment(p3, p4);
    LineSegment ls3 = LineSegment(p3a, p4);

    if (!ls1.intersectsShape(ls2)) {
        cerr << "Test failed:  intersectsShape returned false, but should be true." << endl;
        cerr << ls1 << ", " << ls2 << endl;
        return -1;
    }

    if (ls1.intersectsShape(ls3)) {
        cerr << "Test failed:  intersectsShape returned true, but should be false." << endl;
        cerr << ls1 << ", " << ls3 << endl;
        return -1;
    }

    //Now LineSegment Region intersection
    Region r1 = Region(p5, p6);
    Region r2 = Region(p7, p6);
    Region r3 = Region(p8, p9);
    
    if (!r1.intersectsShape(ls1) || !ls1.intersectsShape(r1)) {
        cerr << "Test failed:  intersectsShape returned false, but should be true." << endl;
        cerr << r1 << ", " << ls1 << endl;
        return -1;
    }

    if (r2.intersectsShape(ls1) || ls1.intersectsShape(r2)) {
        cerr << "Test failed:  intersectsShape returned true, but should be false." << endl;
        cerr << r2 << ", " << ls1 << endl;
        return -1;
    }

    // This is the contains test
    if (!r3.intersectsShape(ls1) || !ls1.intersectsShape(r3)) {
        cerr << "Test failed:  intersectsShape returned false, but should be true." << endl;
        cerr << r3 << ", " << ls1 << endl;
        return -1;
    }

    //line distace test
//    double v1[2]={0,1},v2[2]={1,0},v3[2]={0,1},v4[2]={0,0};
//    double pLow[3]={0,0,0},pHigh[3]={1,1,1};
//    cout<<Trajectory::line2lineIED(TimePoint(v1,0,0,2),TimePoint(v2,1,1,2),
//            TimePoint(v3,0,0,2),TimePoint(v4,1,1,2))<<endl;
//    cout<<Trajectory::line2MBRDistance(TimePoint(v1,0,0,2),TimePoint(v2,1,1,2),
//            Region(pLow,pHigh,3))<<endl;

//    Trajectory traja;
//    traja.loadFromString("17778.000000,17543.000000,485.000000 17787.789971,17567.086437,486.000000 17797.579943,17591.172875,487.000000 17807.369914,17615.259312,488.000000 17817.159885,17639.345749,489.000000 17826.949856,17663.432187,490.000000 17836.739828,17687.518624,491.000000 17840.106177,17712.658699,492.000000 17838.523750,17738.610499,493.000000 17836.941323,17764.562299,494.000000 17835.125230,17790.497241,495.000000 17832.966047,17816.407431,496.000000 17830.806865,17842.317620,497.000000 17824.247663,17867.234211,498.000000 17815.063174,17891.557968,499.000000 17805.878685,17915.881726,500.000000 17796.694195,17940.205483,501.000000 17787.509706,17964.529240,502.000000 17778.325217,17988.852997,503.000000 17769.140728,18013.176754,504.000000 17759.956239,18037.500511,505.000000 17750.771749,18061.824268,506.000000 17741.587260,18086.148025,507.000000 17737.028723,18111.571502,508.000000 17734.284364,18137.426259,509.000000 17731.540004,18163.281017,510.000000 17728.795644,18189.135774,511.000000 17726.051284,18214.990532,512.000000 17723.306925,18240.845289,513.000000 17720.562565,18266.700047,514.000000 17723.487434,18292.373955,515.000000 17727.874085,18318.001231,516.000000 17732.260735,18343.628507,517.000000 17736.647386,18369.255784,518.000000 17743.094715,18394.339212,519.000000 17751.925472,18418.793613,520.000000 17760.756228,18443.248015,521.000000 17769.586984,18467.702416,522.000000 17771.982950,18492.994123,523.000000 17644.786266,18535.148644,524.000000 17517.589582,18577.303165,525.000000 17389.031479,18614.287901,526.000000 17257.570322,18640.248718,527.000000 17125.008294,18647.900622,528.000000 16992.558917,18630.695301,529.000000 16862.708040,18597.608299,530.000000 16734.756447,18558.039141,531.000000 16607.582107,18515.817260,532.000000 16562.000000,18505.000000,533.000000");
//    double plow[3]={16083.3,16481.8,485};
//    double pHigh[3]={17853,17749.1,485};
//    Region* query=new Region(plow,pHigh,3);
//    Region bra;
//    MBC bca;
//    traja.getMBRfull(bra);
//    traja.getMBC(bca);
//    cout<<*query<<endl<<bra<<endl<<bca<<endl;
//    query->intersectsShape(bra);
//    query->intersectsShape(bca);
//    double rlow[3]={11790.8,20103.8,79};
//    double rhigh[3]={12863.3,21411.3,99};
//    double clow[2]={16521.4,21944.1};
//    double chigh[2]={12213.8,18130.7};
//    Region br(rlow,rhigh,3);
//    MBC bc(clow,chigh,54,106,2,1330.78,119.552);
//    std::cerr<<bc.intersectsRegion(br);
}

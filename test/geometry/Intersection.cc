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
//    cout<<Trajectory::line2lineDistance(TimePoint(v1,0,0,2),TimePoint(v2,1,1,2),
//            TimePoint(v3,0,0,2),TimePoint(v4,1,1,2))<<endl;
//    cout<<Trajectory::line2MBRDistance(TimePoint(v1,0,0,2),TimePoint(v2,1,1,2),
//            Region(pLow,pHigh,3))<<endl;

//    Trajectory traja;
//    traja.loadFromString("18411.000000,10663.000000,2.000000 18322.026461,10660.829914,3.000000 18234.093116,10664.296138,4.000000 18151.755363,10\n"
//                         "698.082146,5.000000 18069.417609,10731.868154,6.000000 17987.079856,10765.654163,7.000000 17905.761503,10801.455636,8.00\n"
//                         "0000 17829.754707,10847.758627,9.000000 17753.747911,10894.061617,10.000000 17691.686337,10957.295917,11.000000 17641.97\n"
//                         "6863,11029.943447,12.000000 17603.790806,11110.335145,13.000000 17582.268105,11196.655236,14.000000 17601.970285,11280.3\n"
//                         "58090,15.000000 17638.307347,11361.602277,16.000000 17674.644410,11442.846464,17.000000 17725.915082,11511.713392,18.000\n"
//                         "000 17761.872463,11593.126332,19.000000 17662.762080,11640.464504,20.000000 17552.693620,11715.492146,21.000000 17461.44\n"
//                         "8520,11807.903475,22.000000 17358.742410,11893.970048,23.000000 17283.306144,12003.572296,24.000000 17213.644600,12118.0\n"
//                         "37417,25.000000 17144.800148,12232.999753,26.000000 17075.857717,12347.903805,27.000000 16997.812880,12454.519941,28.000\n"
//                         "000 16872.144546,12482.939030,29.000000 16858.000000,12613.957540,30.000000 16834.993813,12741.625067,31.000000 16766.42\n"
//                         "1798,12856.750559,32.000000 16699.754601,12972.982362,33.000000 16633.736186,13089.590976,34.000000 16567.717771,13206.1\n"
//                         "99590,35.000000 16513.082828,13328.518340,36.000000 16458.112852,13450.724321,37.000000 16402.899820,13572.820397,38.000\n"
//                         "000 16352.022782,13696.693249,39.000000 16305.454042,13822.340981,40.000000 16258.465783,13947.829872,41.000000 16210.07\n"
//                         "1202,14072.785186,42.000000 16159.975361,14197.068760,43.000000 16109.879521,14321.352334,44.000000 16070.000000,14420.0\n"
//                         "00000,45.000000");
//                    double a[3]={13233,9726,2},b[3]={20708,14420,246};
//                    Region rg(a,b,3);
//                    cout<<traja.getMinimumDistance(rg);
}

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

    //Mbbc test
    double sLow[2]={39.996512,116.322447};
    double sHigh[2]={39.996512,116.322447};
    double eLow[2]={40.006551,116.329964};
    double eHigh[2]={40.006551,116.329964};
    double vLow[2]={0.000001,-0.000008};
    double vHigh[2]={0.000030,0.000017};
    double pLow[2]={39.996534,116.322148};
    double pHigh[2]={39.999769,116.324886};
    double oLow[2]={39.996534,116.322148};
    double oHigh[2]={39.999769,116.324886};
    //just for test, should use delete to free pointers
    Mbbc m=Mbbc(Region(sLow,sHigh,2),Region(eLow,eHigh,2),
                     Region(vLow,vHigh,2),Region(pLow,pHigh,2),0.0,900.0);
    Region r(pLow,pHigh,2);
    TimeRegion t(oLow,oHigh,0.5,0.5,2);
//    cout<<m.intersectsTimeRegion(t);
    uint8_t* d;
    uint32_t l;
//    r.storeToByteArray(&d,l);
//    cout<<d[30]<<endl;
//    r.loadFromByteArray(d);
    m.storeToByteArray(&d,l);

    m.loadFromByteArray(d);
//    cout<<m.toString();
    double pc[2];
    vector<TimePoint> points;
    pc[0]=39.990873;pc[1]=116.317939;points.emplace_back(TimePoint(pc,680,680,2));
    pc[0]=39.991146;pc[1]=116.317736;points.emplace_back(TimePoint(pc,685,685,2));
    pc[0]=39.991179;pc[1]=116.317958;points.emplace_back(TimePoint(pc,690,690,2));
    pc[0]=39.991291;pc[1]=116.318188;points.emplace_back(TimePoint(pc,695,695,2));
    pc[0]=39.991339;pc[1]=116.318216;points.emplace_back(TimePoint(pc,700,700,2));
    pc[0]=39.99135;pc[1]=116.318202;points.emplace_back(TimePoint(pc,705,705,2));
    pc[0]=39.991341;pc[1]=116.318231;points.emplace_back(TimePoint(pc,710,710,2));
    pc[0]=39.991326;pc[1]=116.318289;points.emplace_back(TimePoint(pc,715,715,2));
    pc[0]=39.991324;pc[1]=116.318438;points.emplace_back(TimePoint(pc,720,720,2));
    pc[0]=39.991322;pc[1]=116.31864;points.emplace_back(TimePoint(pc,725,725,2));
    pc[0]=39.991327;pc[1]=116.318837;points.emplace_back(TimePoint(pc,730,730,2));
    pc[0]=39.991313;pc[1]=116.318953;points.emplace_back(TimePoint(pc,735,735,2));
    pc[0]=39.991345;pc[1]=116.319096;points.emplace_back(TimePoint(pc,740,740,2));
    pc[0]=39.991377;pc[1]=116.31919;points.emplace_back(TimePoint(pc,745,745,2));
    pc[0]=39.991434;pc[1]=116.31931;points.emplace_back(TimePoint(pc,750,750,2));
    pc[0]=39.991416;pc[1]=116.319416;points.emplace_back(TimePoint(pc,755,755,2));
    pc[0]=39.991453;pc[1]=116.31961;points.emplace_back(TimePoint(pc,760,760,2));
    pc[0]=39.991559;pc[1]=116.319749;points.emplace_back(TimePoint(pc,765,765,2));
    pc[0]=39.991624;pc[1]=116.319851;points.emplace_back(TimePoint(pc,770,770,2));
    pc[0]=39.991617;pc[1]=116.319928;points.emplace_back(TimePoint(pc,775,775,2));
    pc[0]=39.991737;pc[1]=116.320008;points.emplace_back(TimePoint(pc,780,780,2));
    pc[0]=39.99179;pc[1]=116.320071;points.emplace_back(TimePoint(pc,785,785,2));
    pc[0]=39.991922;pc[1]=116.320241;points.emplace_back(TimePoint(pc,790,790,2));
    pc[0]=39.992099;pc[1]=116.320387;points.emplace_back(TimePoint(pc,795,795,2));
    pc[0]=39.992271;pc[1]=116.320474;points.emplace_back(TimePoint(pc,800,800,2));
    pc[0]=39.992403;pc[1]=116.320479;points.emplace_back(TimePoint(pc,805,805,2));
    pc[0]=39.992503;pc[1]=116.320392;points.emplace_back(TimePoint(pc,810,810,2));
    pc[0]=39.992621;pc[1]=116.320306;points.emplace_back(TimePoint(pc,815,815,2));
    pc[0]=39.992755;pc[1]=116.320391;points.emplace_back(TimePoint(pc,820,820,2));
    pc[0]=39.992876;pc[1]=116.32036;points.emplace_back(TimePoint(pc,825,825,2));
    pc[0]=39.993004;pc[1]=116.320402;points.emplace_back(TimePoint(pc,830,830,2));
    pc[0]=39.993094;pc[1]=116.320398;points.emplace_back(TimePoint(pc,835,835,2));
    pc[0]=39.993234;pc[1]=116.320376;points.emplace_back(TimePoint(pc,840,840,2));
    pc[0]=39.993469;pc[1]=116.320318;points.emplace_back(TimePoint(pc,845,845,2));
    pc[0]=39.993673;pc[1]=116.320364;points.emplace_back(TimePoint(pc,850,850,2));
    pc[0]=39.993807;pc[1]=116.320333;points.emplace_back(TimePoint(pc,855,855,2));
    pc[0]=39.993969;pc[1]=116.320326;points.emplace_back(TimePoint(pc,860,860,2));
    pc[0]=39.994174;pc[1]=116.320295;points.emplace_back(TimePoint(pc,865,865,2));
    pc[0]=39.994401;pc[1]=116.320241;points.emplace_back(TimePoint(pc,870,870,2));
    pc[0]=39.994525;pc[1]=116.320202;points.emplace_back(TimePoint(pc,875,875,2));
    pc[0]=39.994604;pc[1]=116.320265;points.emplace_back(TimePoint(pc,880,880,2));
    pc[0]=39.994705;pc[1]=116.32041;points.emplace_back(TimePoint(pc,885,885,2));
    pc[0]=39.994748;pc[1]=116.320628;points.emplace_back(TimePoint(pc,890,890,2));
    pc[0]=39.994779;pc[1]=116.320831;points.emplace_back(TimePoint(pc,895,895,2));

    Trajectory traj(points);
//    cout<<"m and r is"<<m.toString()<<r.toString();
    double qLow[2]={39.993017,116.320135};
    double qHigh[2]={39.994017,116.321135};
    TimeRegion tm=TimeRegion(qLow,qHigh,855,855,2);
//    cout<<traj.getMinimumDistance(m);
//    cout<<"end here\n\n\n\n";
//    cout<<traj.getMinimumDistance(r);

    Trajectory traja,trajb;
    string s1="11535.221554,12974.107786,75.000000 11535.221554,12974.107786,75.000000 11505.816516,12985.490381,76.000000 11476.411478\n"
              ",12996.872976,77.000000 11447.006440,13008.255572,78.000000 11417.601402,13019.638167,79.000000 11388.007022,13030.47669\n"
              "8,80.000000 11357.888518,13039.809192,81.000000 11327.770014,13049.141686,82.000000 11297.651510,13058.474180,83.000000\n"
              "11267.533006,13067.806674,84.000000 11237.414502,13077.139168,85.000000 11207.295998,13086.471663,86.000000 11177.177494\n"
              ",13095.804157,87.000000 11147.058990,13105.136651,88.000000 11116.940486,13114.469145,89.000000 11086.821982,13123.80163\n"
              "9,90.000000 11056.703478,13133.134133,91.000000 11026.584974,13142.466628,92.000000 10996.466471,13151.799122,93.000000\n"
              "10966.347967,13161.131616,94.000000 10936.229463,13170.464110,95.000000 10906.110959,13179.796604,96.000000 10875.992455\n"
              ",13189.129099,97.000000 10845.873951,13198.461593,98.000000 10815.755447,13207.794087,99.000000";
    traja.loadFromString(s1);
    Mbbc bc;
    traja.getMbbc(bc,true);
//    cout<<"traj dist is"<<trajb.getMinimumDistance(trajb)<<endl;
//    cout<<"bc dist is"<<traja.getMinimumDistance(bc)<<endl;
    return 0;
}

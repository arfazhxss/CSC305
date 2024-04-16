////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v)
{
    // TODO
    return 0;
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans)
{
    // TODO
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    // TODO
    Point outside(0, 0);
    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    // TODO
    return true;
}

////////////////////////////////////////////////////////////////////////////////

struct Compare
{
    Point p0; // Leftmost point of the poly
    bool operator()(const Point &p1, const Point &p2)
    {
        double angle1 = std::atan2(p1.imag(), p1.real());
        double angle2 = std::atan2(p2.imag(), p2.real());
        return angle1 < angle2;
    }
};

bool inline salientAngle(Point &a, Point &b, Point &c)
{
    double aX = a.real();
    double aY = a.imag();
    double bX = b.real();
    double bY = b.imag();
    double cX = c.real();
    double cY = c.imag();

     

    // TODO
    if(((aX * (bY - cY)) + (bX * (cY - aY)) + (cX * (aY - bY))) > 0){
        return true;
    }
    return false;
}

Polygon convex_hull(std::vector<Point> &points)
{
    Compare order;
    // TODO

    int ymin = points[0].imag();
    int minIndex = 0;
    for (int i = 1; i < points.size(); ++i) {
        int y = points[i].imag();
        if (y < ymin || ymin == y && points[i].real() < points[minIndex].real()) {
            ymin = points[i].imag();
            minIndex = i;
        }
    }

    std::swap(points[0], points[minIndex]);

    order.p0 = Point(0, 0);
    std::sort(points.begin()+1, points.end(), order);

    Polygon hull;
    // TODO

    std::stack<Point> Stack;
    Stack.push(points[0]);
    Stack.push(points[1]);
    Stack.push(points[2]);


    // Process the remaining points
    for (int i = 3; i < points.size(); ++i) {
        while (Stack.size() > 1) {
            Point p2 = Stack.top();
            Stack.pop();
            Point p1 = Stack.top();

            if (salientAngle(p1, p2, points[i]) == true) {
                Stack.push(p2);
                break;
            }
        }
         Stack.push(points[i]);
    }

    // Store the points from the stack into a vector (convex hull)
    while (!Stack.empty()) {
        hull.push_back(Stack.top());
        Stack.pop();
    }

    return hull;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename)
{
    std::vector<Point> points;
    std::ifstream in(filename);
    std::string line;

    double x,y,z;

    std::getline(in, line);
    while (std::getline(in, line)) {

        in >> x >> y >> z;
        Point point(x,y);
        points.push_back(point);

}    
// TODO
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Point> &points)
{

}

Polygon load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    // TODO
    return {};
}

void save_obj(const std::string &filename, Polygon &poly)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("failed to open file " + filename);
    }
    out << std::fixed;
    for (const auto &v : poly)
    {
        out << "v " << v.real() << ' ' << v.imag() << " 0\n";
    }
    for (size_t i = 0; i < poly.size(); ++i)
    {
        out << "l " << i + 1 << ' ' << 1 + (i + 1) % poly.size() << "\n";
    }
    out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Point> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    // Polygon poly = load_obj(poly_path);
    // std::vector<Point> result;
    // for (size_t i = 0; i < points.size(); ++i)
    // {
    //     if (is_inside(poly, points[i]))
    //     {
    //         result.push_back(points[i]);
    //     }
    // }
    // save_xyz("output.xyz", result);

    // ////////////////////////////////////////////////////////////////////////////////
    // //Convex hull

    Polygon hull = convex_hull(points);
    save_obj("output.obj", hull);

    return 0;
}

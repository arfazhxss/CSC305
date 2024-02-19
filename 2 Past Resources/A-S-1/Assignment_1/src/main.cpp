#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

const std::string root_path = "../data";

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point& u, const Point& v)
{
    return u.real() * v.imag() - u.imag() * v.real();
}

// Return true if [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point& a, const Point& b, const Point& c, const Point& d, Point& ans)
{
    double det1 = det(b - a, c - a);
    double det2 = det(b - a, d - a);
    double det3 = det(d - c, a - c);
    double det4 = det(d - c, b - c);

    if (((det1 > 0 && det2 < 0) || (det1 < 0 && det2 > 0)) &&
        ((det3 > 0 && det4 < 0) || (det3 < 0 && det4 > 0)))
    {
        ans = a + (b - a) * (det3 / (det3 - det4));
        return true;
    }

    return false;
}

bool is_inside(const Polygon& poly, const Point& query)
{
    // 1. Compute bounding box and set the coordinate of a point outside the polygon
    double maxX = -std::numeric_limits<double>::infinity();
    double maxY = -std::numeric_limits<double>::infinity();
    double minX = std::numeric_limits<double>::infinity();
    double minY = std::numeric_limits<double>::infinity();

    for (const Point& p : poly)
    {
        maxX = std::max(maxX, p.real());
        maxY = std::max(maxY, p.imag());
        minX = std::min(minX, p.real());
        minY = std::min(minY, p.imag());
    }

    Point outside(maxX + 1, maxY + 1);

    // 2. Cast a ray from the query point to the 'outside' point and count the number of intersections
    int intersectionCount = 0;

    for (size_t i = 0; i < poly.size(); ++i)
    {
        const Point& p1 = poly[i];
        const Point& p2 = poly[(i + 1) % poly.size()];

        if (p1.imag() == p2.imag())
            continue; // Skip horizontal edges

        if (query.imag() >= std::min(p1.imag(), p2.imag()))
        {
            if (query.imag() <= std::max(p1.imag(), p2.imag()))
            {
                if (query.real() <= std::max(p1.real(), p2.real()))
                {
                    if (p1.imag() != p2.imag())
                    {
                        double xIntersection = (query.imag() - p1.imag()) * (p2.real() - p1.real()) /
                                                (p2.imag() - p1.imag()) + p1.real();

                        if (p1.real() == p2.real() || query.real() <= xIntersection)
                            ++intersectionCount;
                    }
                }
            }
        }
    }

    return intersectionCount % 2 == 1;
}

struct Compare
{
    Point p0; // Leftmost point of the poly

    bool operator()(const Point& p1, const Point& p2)
    {
        double angle1 = std::arg(p1 - p0);
        double angle2 = std::arg(p2 - p0);

        if (angle1 == angle2)
            return std::abs(p1 - p0) < std::abs(p2 - p0);
        else
            return angle1 < angle2;
    }
};

bool inline salientAngle(const Point& a, const Point& b, const Point& c)
{
    double angle1 = std::arg(b - a);
    double angle2 = std::arg(c - b);

    double angleDifference = std::fmod(angle2 - angle1 + 5 * M_PI, 2 * M_PI);

    return angleDifference > M_PI;
}

Polygon convex_hull(std::vector<Point>& points)
{
    if (points.size() < 3)
        return points;

    // Find the leftmost point as the reference point
    auto it = std::min_element(points.begin(), points.end(),
                               [](const Point& p1, const Point& p2) {
                                   return p1.real() < p2.real() || (p1.real() == p2.real() && p1.imag() < p2.imag());
                               });
    Point p0 = *it;

    // Sort the points based on the angle with respect to p0
    std::sort(points.begin(), points.end(),
              [&p0](const Point& p1, const Point& p2) {
                  double angle1 = std::arg(p1 - p0);
                  double angle2 = std::arg(p2 - p0);
                  if (angle1 == angle2)
                      return std::abs(p1 - p0) < std::abs(p2 - p0);
                  else
                      return angle1 < angle2;
              });

    // Build the convex hull
    Polygon hull;
    hull.push_back(points[0]);
    hull.push_back(points[1]);

    for (size_t i = 2; i < points.size(); ++i)
    {
        while (hull.size() >= 2 &&
               !salientAngle(hull[hull.size() - 2], hull[hull.size() - 1], points[i]))
        {
            hull.pop_back();
        }

        hull.push_back(points[i]);
    }

    return hull;
}


std::vector<Point> load_xyz(const std::string& filename)
{
    std::string line;
    
    std::ifstream in(filename);
    std::getline(in, line);
    if (!in)
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return {};
    }

    std::vector<Point> points;
    double x, y, z;

    while (in >> x >> y >> z)
    {
        points.emplace_back(x, y);
    }

    in.close();

    return points;
}

void save_xyz(const std::string& filename, const std::vector<Point>& points)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const Point& point : points)
    {
        out << point.real() << " " << point.imag() << " 0.0\n";
    }

    out.close();
}

Polygon load_obj(const std::string& filename)
{
    Polygon polygon;
    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return polygon;
    }

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue; // Skip empty lines and comments

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "v")
        {
            double x, y, z;
            iss >> x >> y >> z;
            Point point(x, y);
            polygon.push_back(point);
        }
    }

    in.close();

    return polygon;
}

void save_obj(const std::string& filename, const Polygon& poly)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    out << std::fixed;
    for (const auto& v : poly)
    {
        out << "v " << v.real() << ' ' << v.imag() << " 0\n";
    }
    for (size_t i = 0; i < poly.size(); ++i)
    {
        out << "l " << i + 1 << ' ' << 1 + (i + 1) % poly.size() << "\n";
    }
    out << std::endl;
}

int main(int argc, char* argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Point> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    // Point in polygon
    Polygon poly = load_obj(poly_path);
    std::vector<Point> result;
    for (const Point& point : points)
    {
        if (is_inside(poly, point))
        {
            result.push_back(point);
        }
    }
    save_xyz("output.xyz", result);

    ////////////////////////////////////////////////////////////////////////////////
    // Convex hull
    Polygon hull = convex_hull(points);
    save_obj("output.obj", hull);

    return 0;
}

#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>

using namespace Eigen;

const std::string root_path = DATA_DIR;

double inline det(const Vector2d &u, const Vector2d &v) { return u.x() * v.y() - u.y() * v.x(); }

bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    Vector2d AB = b - a;
    Vector2d CD = d - c;
    Vector2d AC = c - a;

    // Check if the line segments are not collinear
    // If not, calculate the parameter values for the intersection point
    if (det(AB, CD) != 0)
    {
        double t = det(AC, CD) / det(AB, CD);
        double u = det(AC, AB) / det(AB, CD);
        return (t >= 0 && t <= 1 && u >= 0 && u <= 1);
    }

    return false; 
}

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    double maxX = poly[0].x();
    double maxY = poly[0].y();
    double minX = poly[0].x();
    double minY = poly[0].y();

    for (const auto &point : poly)
    {
        minX = std::min(minX, point.x());
        maxX = std::max(maxX, point.x());
        minY = std::min(minY, point.y());
        maxY = std::max(maxY, point.y());
    }

    Vector2d outside(maxX + 1.0, query(0)); // Set x-coordinate just outside the bounding box

    // 2. Cast a ray from the query point to the 'outside' point, count the number of intersections
    int intersectionCount = 0;

    for (size_t i = 0; i < poly.size(); ++i)
    {
        Vector2d p1 = poly[i];
        Vector2d p2 = poly[(i + 1) % poly.size()];

        // Check if the ray intersects the edge
        if ((p1(1) > query(1)) != (p2(1) > query(1)) &&
            (query(0) < (p2(0) - p1(0)) * (query(1) - p1(1)) / (p2(1) - p1(1)) + p1(0)))
        {
            intersectionCount++;
        }
    }

    // 3. Determine if the query point is inside the polygon based on the number of intersections
    return (intersectionCount % 2) == 1;
}

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error opening file: " << filename <<std::endl;
        return points;
    }

    double x, y;
    while (in >> x >> y) {
        points.push_back(Vector2d(x, y));
    }

    in.close();
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::string full_path = root_path + filename;

    std::ofstream out(full_path);

    if (!out.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const auto &point : points) { out << point.x() << " " << point.y() << " 0" << std::endl; }
    
    out.close();
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);
    std::vector<Vector2d> poly = load_obj(poly_path);

    std::vector<Vector2d> result;

    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);
    return 0;
}

// C++ include
#include <iostream>
#include <string>
#include <vector>
#include "utils.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d cameraOrigin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, 
    // pointing in the direction -z and 
    // covering the unit square (-1,1) in x and y
    // 
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    const Vector3d lightPosition(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d rayOrigin = pixel_center;
            const Vector3d rayDirection = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(rayOrigin(0), rayOrigin(1));
            const double sphereRadius = 0.9;

            if (ray_on_xy.norm() < sphereRadius)
            {
                Vector3d rayIntersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphereRadius * sphereRadius - ray_on_xy.squaredNorm()));

                Vector3d normalRay = rayIntersection.normalized();
                C(i, j) = (lightPosition - rayIntersection).normalized().transpose() * normalRay;
                C(i, j) = std::max(C(i, j), 0.);
                A(i, j) = 1;
            }
        }
    }

    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d cameraOrigin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgramOrigin (-0.5, -0.5, 0);
    const Vector3d pgramU (0, 0.7, -10);
    const Vector3d pgramV (1, 0.4, 0);

    // Cross product between two sides  
    const Vector3d pgramNormal = pgramU.cross(pgramV);

    const Vector3d a = pgramOrigin;
    const Vector3d b = pgramOrigin + pgramU;
    const Vector3d c = pgramOrigin + pgramV;

    // Single light source
    const Vector3d lightPosition(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;
            const Vector3d rayOrigin = pixel_center;
            const Vector3d rayDirection = camera_view_direction;

            const Vector3d e = rayOrigin; // rayOrigin = e
            const Vector3d d = camera_view_direction; // rayDirection = d

            // Calculate the values of u, v, and t using the equation
            // Formula: rayOrigin + ( t × rayDirection ) = ( pgramOrigin ) + ( u × pgramU ) + ( v × pgramV )
            //
            // matrixU_X_axis = (ParalellogramOriginX-RayOriginX), (ParalellogramXWidth), (CameraViewDirectionX)
            // matrixV_X_axis = (ParalellogramYWidth), (ParalellogramOriginX-RayOriginX), (CameraViewDirectionX)
            // matrixT_X_axis = (ParalellogramYWidth), (ParalellogramWidthX), (ParalellogramOriginX-RayOriginX)
            Matrix3d matrixU, matrixV, matrixT, matrixA;
            matrixU <<  a(0) - e(0), a(0) - c(0), d(0),
                        a(1) - e(1), a(1) - c(1), d(1),
                        a(2) - e(2), a(2) - c(2), d(2);
            matrixV <<  a(0) - b(0), a(0) - e(0), d(0),
                        a(1) - b(1), a(1) - e(1), d(1),
                        a(2) - b(2), a(2) - e(2), d(2);
            matrixA <<  a(0) - b(0), a(0) - c(0), d(0),
                        a(1) - b(1), a(1) - c(1), d(1),
                        a(2) - b(2), a(2) - c(2), d(2);
            matrixT <<  a(0) - b(0), a(0) - c(0), a(0) - e(0),
                        a(1) - b(1), a(1) - c(1), a(1) - e(1),
                        a(2) - b(2), a(2) - c(2), a(2) - e(2);

            double detU = matrixU.determinant();
            double detV = matrixV.determinant();
            double detT = matrixT.determinant();
            double detA = matrixA.determinant();

            double u, v, t;
            u = detU / detA;
            v = detV / detA;
            t = detT / detA;

            // if the ray hits the parallelogram:
            if (t > 0 && u > 0 && v > 0 && u < 1 && v < 1)
            {
                Vector3d rayIntersection = rayOrigin + t * rayDirection;
                Vector3d normalRay = pgramNormal.normalized();

                C(i, j) = -(lightPosition - rayIntersection).normalized().dot(normalRay);
                C(i, j) = std::max(C(i, j), 0.0);
                A(i, j) = 1;
            }
        }
    }

    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);
    const Vector3d center(0, 0, 1);
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);
    const Vector3d pgramOrigin(-0.5, -0.5, 0);
    const Vector3d pgramU(0, 0.7, -10);
    const Vector3d pgramV(1, 0.4, 0);
    const Vector3d pgramNormal = pgramU.cross(pgramV);

    const Vector3d a = pgramOrigin;
    const Vector3d b = pgramOrigin + pgramU;
    const Vector3d c = pgramOrigin + pgramV;

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;
            const Vector3d rayOrigin = camera_origin;
            const Vector3d rayDirection = pixel_center - camera_origin;

            const Vector3d e = camera_origin;
            const Vector3d d = rayDirection;

            Matrix3d matrixU, matrixV, matrixT, matrixA;
            matrixU <<  a(0) - e(0), a(0) - c(0), d(0),
                        a(1) - e(1), a(1) - c(1), d(1),
                        a(2) - e(2), a(2) - c(2), d(2);
            matrixV <<  a(0) - b(0), a(0) - e(0), d(0),
                        a(1) - b(1), a(1) - e(1), d(1),
                        a(2) - b(2), a(2) - e(2), d(2);
            matrixA <<  a(0) - b(0), a(0) - c(0), d(0),
                        a(1) - b(1), a(1) - c(1), d(1),
                        a(2) - b(2), a(2) - c(2), d(2);
            matrixT <<  a(0) - b(0), a(0) - c(0), a(0) - e(0),
                        a(1) - b(1), a(1) - c(1), a(1) - e(1),
                        a(2) - b(2), a(2) - c(2), a(2) - e(2);

            double detU = matrixU.determinant();
            double detV = matrixV.determinant();
            double detT = matrixT.determinant();
            double detA = matrixA.determinant();

            double u, v, t;
            u = detU / detA;
            v = detV / detA;
            t = detT / detA;

            if (t > 0 && u > 0 && v > 0 && u < 1 && v < 1)
            {
                Vector3d ray_intersection = rayOrigin + t * rayDirection;
                Vector3d normalRay = pgramNormal.normalized();
                C(i, j) = -(light_position - ray_intersection).normalized().transpose() * normalRay;
                C(i, j) = std::max(C(i, j), 0.0);
                A(i, j) = 1;
            }
        }
    }

    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the red
    MatrixXd G = MatrixXd::Zero(800, 800); // Store the green
    MatrixXd B = MatrixXd::Zero(800, 800); // Store the blue
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d cameraOrigin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    const Vector3d sphereCentre(0, 0, 0);
    const double sphereRadius = 0.9;

    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specularColor(0, 0, 1);

    const Vector3d lightPosition(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            const Vector3d rayOrigin = cameraOrigin;
            const Vector3d rayDirection = (pixel_center - cameraOrigin).normalized();
            const double sphereRadius = 0.9;
            Vector3d e = rayOrigin;
            Vector3d d = rayDirection;
            double r = sphereRadius;

            double A1 = rayDirection.dot(rayDirection);
            double B1 = 2 * (rayDirection.dot(rayOrigin - sphereCentre));
            double C1 = ((rayOrigin - sphereCentre).dot(rayOrigin - sphereCentre)) - (sphereRadius * sphereRadius);

            double delta = sqrt(B1 * B1 - 4 * A1 * C1);

            if (delta >= 0)
            {
                double t = (-B1 - delta) / (2 * A1);
                Vector3d rayIntersection = rayOrigin + t * rayDirection;
                Vector3d normalRay = (rayIntersection - sphereCentre).normalized();
                Vector3d n = normalRay;
                Vector3d l = (lightPosition - rayIntersection).normalized();
                Vector3d v = (cameraOrigin - rayIntersection).normalized();
                Vector3d h = (v + l) / ((v + l).norm());

                Vector3d diffuse = diffuse_color * l.dot(n);
                Vector3d specular = specularColor * pow(n.dot(h), specular_exponent);

                R(i, j) = ambient + diffuse[0] + specular[0];
                R(i, j) = std::max(R(i, j), 0.0);

                G(i, j) = ambient + diffuse[1] + specular[1];
                G(i, j) = std::max(G(i, j), 0.0);

                B(i, j) = ambient + diffuse[2] + specular[2];
                B(i, j) = std::max(B(i, j), 0.0);

                A(i, j) = 1;
            }
        }
    }

    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}

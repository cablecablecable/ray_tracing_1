#include <complex>
#include <iostream>
#include "vec3.h"
#include "color.h"
#include "ray.h"

double hit_sphere(const point3& center, double radius, const ray& r)
{
    vec3 oc = center - r.origin();
    auto a = dot(r.direction(), r.direction());
    auto b = -2.0 * dot(r.direction(), oc);
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;

    if (discriminant >= 0)
        return (-b - std::sqrt(discriminant)) / (2.0 * a);

    return -1.0;
}

color ray_color(const ray& r)
{
    double t = hit_sphere(point3(0,0,-1), 0.5, r);

    if (t >= 0)
    {
        point3 intersection = r.at(t);
        vec3 normal = vec3(intersection - point3(0,0,-1));
        normal = unit_vector(normal);

        return 0.5 * color(normal.x() + 1.0, normal.y() + 1.0, normal.z() + 1.0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    double a = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
}

int main()
{
    //image
    double aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;

    // calculate the image height, and ensure that it's at least 1
    int image_height = (int) (image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;


    // camera
    double focal_length = 1.0;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * (double) image_width / image_height;
    point3 camera_center = point3(0, 0, 0);


    // calculate vectors across horizontal and down the vertical viewport edges
    vec3 viewport_u = vec3(viewport_width, 0, 0);
    vec3 viewport_v = vec3(0, -viewport_height, 0);


    // calculate the horizontal and vertical delta vectors from pixel to pixel
    vec3 pixel_delta_u = viewport_u / image_width;
    vec3 pixel_delta_v = viewport_v / image_height;


    // calculate location of upper left pixel
    point3 viewport_upper_left =
        camera_center - vec3(0, 0, focal_length) - (viewport_u / 2) - (viewport_v / 2);

    point3 pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);


    //render
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";


    for (int j = 0; j < image_height; j++)
    {
        std::clog << "\rScanlines remaining: " << (image_height - 1) << " " << std::flush;

            for (int i = 0; i < image_width; i++)
        {
            vec3 pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
            vec3 ray_direction = pixel_center - camera_center;
            ray r(camera_center, ray_direction);

            color pixel_color = ray_color(r);
            write_color(std::cout, pixel_color);
        }
    }

    std::clog << "\rDone.                   \n";
}
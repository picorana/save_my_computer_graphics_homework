#include "scene.h"
#include "intersect.h"

void print_vector(vec3f vec) {
	message("%f %f %f\n", vec.x, vec.y, vec.z);
}

intersection3f intersect_quad(const Scene &scene, ray3f &ray, const Surface &surface, float mindistance) {
	
	intersection3f interx = intersection3f();

	//	ray r = E + tD
	//	E -> ray origin
	//	D -> ray direction
	//	
	//	point on a ray: P(t) = E + tD
	//	point on a plane: (P(t) - C) * n = 0
	//	so (E + tD - C) * n = 0
	//	t = [(C - E) * n] / d * n

	// compute ray intersection (and ray parameter), continue if not hit
	if (dot(ray.d, surface.frame.z) == 0) {
		interx.hit = false;
		return interx;
	}

	//auto t = -ray.e.z / ray.d.z;
	auto t = dot(surface.frame.o - ray.e, surface.frame.z) / dot(ray.d, surface.frame.z);
	

	// check if computed param is within ray.tmin and ray.tmax
	if (t < ray.tmin or t > ray.tmax) {
		interx.hit = false;
		return interx;
	}

	auto p = transform_point(surface.frame, ray.eval(t));
	//auto p = ray.eval(t);
	/*if (surface.frame.o.x + surface.radius < p.x or surface.frame.o.x - surface.radius > p.x
		or surface.frame.o.y + surface.radius < p.y or surface.frame.o.y - surface.radius> p.y) {
		interx.hit = false;
		return interx;
	}*/

	if (abs(p.x) > surface.radius or abs(p.y) > surface.radius) {
		interx.hit = false;
		return interx;
	}
	
	p = transform_point_inverse(surface.frame, p);
	// if hit, set intersection record values
	if (t < mindistance) {

		auto pl = normalize(surface.frame.z);

		interx.ray_t = t;
		interx.hit = true;
		interx.mat = surface.mat;
		interx.norm = pl;
		//interx.norm = transform_normal(surface.frame, pl);
		interx.pos = p;
		//interx.pos = transform_point(scene.camera->frame, p);
		mindistance = t;
	}
	
	return interx;
}

intersection3f intersect_sphere(const Scene &scene, ray3f &ray, const Surface &surface, float mindistance) {

	intersection3f interx = intersection3f();

	auto a = lengthSqr(ray.d);
	auto b = 2 * dot(ray.d, ray.e - surface.frame.o);
	auto c = lengthSqr(ray.e - surface.frame.o) - surface.radius*surface.radius;
	auto d = b*b - 4 * a*c;

	if (d < 0) {
		interx.hit = false;
		return interx;
	}

	auto tmin = (-b - sqrt(d)) / (2 * a);
	auto tmax = (-b + sqrt(d)) / (2 * a);

	float t;

	// check if computed param is within ray.tmin and ray.tmax
	if (tmin >= ray.tmin && tmin <= ray.tmax) t = tmin;
	else if (tmax >= ray.tmin && tmax <= ray.tmax) t = tmax;
	else {
		interx.hit = false;
		return interx;
	}

	// compute ray intersection (and ray parameter), continue if not hit
	//auto p = transform_point(surface.frame, ray.eval(t));
	auto p = ray.eval(t);
	auto pl = normalize(p - surface.frame.o);

	// if hit, set intersection record values
	if (t < mindistance) {
		interx.ray_t = t;
		interx.hit = true;
		interx.mat = surface.mat;
		interx.pos = p;
		interx.norm = pl;
		mindistance = t;
	}

	return interx;
}

// intersects the scene's surfaces and return the first intrerseciton (used for raytracing homework)
intersection3f intersect_surfaces(Scene* scene, ray3f ray) {
    
	// create a default intersection record to be returned
    auto intersection = intersection3f();

	//ray = transform_ray(scene->camera->frame, ray);

	auto mindistance = INFINITY;

    // foreach surface
	for (Surface* surface : scene->surfaces) {
		
		if (surface->isquad) {		// if it is a quad
			intersection3f intersectiontmp = intersect_quad(*scene, ray, *surface, mindistance);
			if (intersectiontmp.hit) intersection = intersectiontmp;
		} else {					// if it is a sphere
			intersection3f intersectiontmp = intersect_sphere(*scene, ray, *surface, mindistance);
			if (intersectiontmp.hit) intersection = intersectiontmp;
		}
	}
                
    // record closest intersection
    return intersection;
}


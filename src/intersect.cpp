#include "scene.h"
#include "intersect.h"

intersection3f intersect_quad(const Scene &scene, const ray3f &ray, const Surface &surface, float mindistance) {
	
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
	if (ray.d.z == 0) {
		interx.hit = false;
		return interx;
	}

	auto t = -ray.e.z / ray.d.z;

	// check if computed param is within ray.tmin and ray.tmax
	if (t < ray.tmin or t > ray.tmax) {
		interx.hit = false;
		return interx;
	}

	auto p = ray.eval(t);
	if (surface.frame.o.x + surface.radius < p.x or surface.frame.o.x - surface.radius > p.x
		or surface.frame.o.y + surface.radius < p.y or surface.frame.o.y - surface.radius> p.y) {
		interx.hit = false;
		return interx;
	}

	// if hit, set intersection record values
	if (t < mindistance) {
		interx.ray_t = t;
		interx.hit = true;
		interx.mat = surface.mat;
		interx.norm = surface.frame.z;
		interx.pos = p;
		mindistance = t;
	}

	return interx;
}

// intersects the scene's surfaces and return the first intrerseciton (used for raytracing homework)
intersection3f intersect_surfaces(Scene* scene, ray3f ray) {
    
	// create a default intersection record to be returned
    auto intersection = intersection3f();

	ray = transform_ray(scene->camera->frame, ray);

	auto mindistance = INFINITY;

    // foreach surface
	for (Surface* surface : scene->surfaces) {
		// if it is a quad
		if (surface->isquad) {
			intersection3f intersectiontmp = intersect_quad(*scene, ray, *surface, mindistance);
			if (intersectiontmp.hit) intersection = intersectiontmp;
		}
		else { // if it is a sphere

			// NOPE CHECK EVERYTHING
			// CHECK. EVERYTHING.
			auto a = lengthSqr(ray.d);
			auto b = 2 * dot(ray.d, ray.e - surface->frame.o);
			auto c = lengthSqr(ray.e - surface->frame.o) - surface->radius*surface->radius;
			auto d = b*b - 4 * a*c;
			
			if (d < 0) {
				intersection.hit = false;
				//return intersection;
				break;
			}

			auto tmin = (-b - sqrt(d)) / (2 * a);
			auto tmax = (-b + sqrt(d)) / (2 * a);

			float t;

			// check if computed param is within ray.tmin and ray.tmax
			if (tmin >= ray.tmin && tmin <= ray.tmax) t = tmin;
			else if (tmax >= ray.tmin && tmax <= ray.tmax) t = tmax;
			else {
				intersection.hit = false;
				//return intersection;
				break;
			}

			// compute ray intersection (and ray parameter), continue if not hit
			auto p = ray.eval(t);
			auto pl = (p - surface->frame.o) / surface->radius;

			auto phi = atan2(pl.y, pl.x);
			auto theta = acos(pl.z);
			auto ct = cos(theta);
			auto st = sin(theta);
			auto cp = cos(phi);
			auto sp = sin(phi);

			frame3f f = frame3f();
			f.o = p;
			f.x = vec3f(sp, cp, 0);
			f.y = vec3f(ct*cp, ct*sp, st);
			f.z = pl;

			// if hit, set intersection record values
			vec3f aa = transform_point(scene->camera->frame, intersection.pos);
			if (t < mindistance) {
				intersection.ray_t = t;
				intersection.hit = true;
				intersection.mat = surface->mat;
				intersection.pos = p;
				intersection.norm = f.z;
				mindistance = t;
			}

		}
	}

            // check if this is the closest intersection, continue if not
                
    // record closest intersection
    return intersection;
}


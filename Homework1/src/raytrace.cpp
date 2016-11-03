#include "scene.h"
#include "intersect.h"
#include "vmath.h"


// compute the color corresponing to a ray by raytracing
vec3f raytrace_ray(Scene* scene, ray3f ray) {
    
	vec3f res = zero3f;

	// get scene intersection
	intersection3f intersection = intersect_surfaces(scene, ray);
    
	// if not hit, return background
	if (not intersection.hit) return scene->background;
    
	// accumulate color starting with ambient
	res += (scene->ambient * intersection.mat->kd);

	// foreach light
	for (auto light : scene->lights) {
		/*
		*	LAMBERT lighting model
		*	----------------------
		*	Cd = kd * | n * l |
		*/

		/*	
		*	BLINN ligthing model 
		*	--------------------
		*	l -> direction from the point to the light
		*	v -> direction from the point to the viewer
		*	f -> local shading frame that describes surface orientation (normal and tangent)
		* 
		*	cosine of bisector h and normal n
		*	h = (l + v)/|l + v|
		*	brdf(l, v, f) = ks * max(0, n*h)^n
		*	Cs = brdf(l, v, f) * |n * l|
		*/

		/*
		*	POINT LIGHT
		*	-----------
		*	S -> source, P -> intersection point
		*	direction: l = (S - P) / | S - P |
		*	color: L = kl / | S - P |
		*/
		
		vec3f S = light->frame.o;
		vec3f P = intersection.pos;
		vec3f l = normalize(S - P);
		vec3f v = normalize(ray.e - P);
		vec3f h = (l + v) / length(l + v);
		vec3f light_color = light->intensity / distSqr(S , P);
		vec3f light_direction = normalize((S - P)/abs(dist(S, P)));
		float dotres = abs(dot(intersection.norm, h));

		float shadowfloat = 0;
		ray3f shadowray = ray3f();
		shadowray.e = intersection.pos;
		shadowray.d = l;
		intersection3f inter = intersect_surfaces(scene, shadowray);
		if (inter.hit) shadowfloat = 0;
		else shadowfloat = 1;

		vec3f cl = light_color * shadowfloat * (intersection.mat->kd + intersection.mat->ks * pow(max(0.0, dotres), intersection.mat->n)) * abs(dot(l, intersection.norm));

		res += cl;

		
	}

	
	vec3f reflectionvector = zero3f;
	if (not (intersection.mat->kr == zero3f)) {
		ray3f reflectionray = ray3f(intersection.pos, ray.d - intersection.norm * 2 * dot(ray.d, intersection.norm));
		reflectionvector = intersection.mat->kr * raytrace_ray(scene, reflectionray);
	}

	res += reflectionvector;
	
   
	return res;
}

// raytrace an image
image3f raytrace(Scene* scene) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);
    
    // if no anti-aliasing
	if (scene->image_samples == 1) {  
		// foreach pixel
		for (int y = 0; y < scene->image_height ; y++) {
			for (int x = 0; x < scene->image_width; x++) {
				
				// compute ray-camera parameters (u,v) for the pixel
				float u = (x + 0.5f) / scene->image_width;
				float v = (y + 0.5f) / scene->image_height;

				// compute camera ray
				Camera* camera = scene->camera;

				frame3f f = camera->frame;

				auto q = (u - 0.5f) * camera->width * f.x + (v - 0.5f) * camera->height * f.y + -camera->dist * f.z;
				ray3f ray = ray3f(scene->camera->frame.o, normalize(q));

				// set pixel to the color raytraced with the ray
				image.at(x, y) = raytrace_ray(scene, ray);		
			}
		}	
	}
	else {
		int number_of_samples = scene->image_samples;
		for (int y = 0; y < scene->image_height; y++) {
			for (int x = 0; x < scene->image_width; x++) {
				vec3f c = zero3f;
				for (int yy = 0; yy < number_of_samples; yy++) {
					for (int xx = 0; xx < number_of_samples; xx++) {
						float u = (x + (xx+0.5f)/number_of_samples) / scene->image_width;
						float v = (y + (yy + 0.5f) / number_of_samples) / scene->image_height;
						Camera* camera = scene->camera;
						frame3f f = camera->frame;

						auto q = (u - 0.5f) * camera->width * f.x + (v - 0.5f) * camera->height * f.y + -camera->dist * f.z;
						ray3f ray = ray3f(scene->camera->frame.o, normalize(q));
						c += raytrace_ray(scene, ray);
					}
				}
				image.at(x, y) = c / pow(number_of_samples, 2);
			}
		}
	}
        
    return image;
}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "01_raytrace", "raytrace a scene",
            {  {"resolution", "r", "image resolution", "int", true, jsonvalue() }  },
            {  {"scene_filename", "", "scene filename", "string", false, jsonvalue("scene.json")},
               {"image_filename", "", "image filename", "string", true, jsonvalue("")}  }
        });
    auto scene_filename = args.object_element("scene_filename").as_string();
    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    auto scene = load_json_scene(scene_filename);
    if(not args.object_element("resolution").is_null()) {
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    message("rendering %s ... ", scene_filename.c_str());
    auto image = raytrace(scene);
    write_png(image_filename, image, true);
    delete scene;
    message("done\n");
}

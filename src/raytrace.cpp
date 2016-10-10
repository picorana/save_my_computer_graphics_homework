#include "scene.h"
#include "intersect.h"
#include "vmath.h"

vec3f material_brdf() {
	return one3f;
}

vec3f compute_color(const Scene &scene, const Light &light, intersection3f intersection) {
	// distSqr computes distance squared between two points, in this case the origin of the light and the 
	// point of the intersection transformed into light->frame coords
	//return light.intensity / dist(light.frame.o, transform_point_inverse(light.frame, intersection.pos));
	return light.intensity / dist(transform_point_inverse(scene.camera->frame, light.frame.o), intersection.pos);
}

vec3f compute_light_direction(const Scene &scene, const Light &light, intersection3f intersection) {
	// direction: l = (S - P) / | S - P |
	//vec3f point_in_light_frame = transform_point_inverse(light.frame, intersection.pos);
	//return (light.frame.o - point_in_light_frame) / dist(light.frame.o, point_in_light_frame);
	vec3f light_src_in_camera_frame = transform_point_inverse(scene.camera->frame, light.frame.o);
	return (light_src_in_camera_frame - intersection.pos) / dist(light_src_in_camera_frame, intersection.pos);
}

// compute the color corresponing to a ray by raytracing
vec3f raytrace_ray(Scene* scene, ray3f ray) {
    
	vec3f res = zero3f;

	// get scene intersection
	intersection3f intersection = intersect_surfaces(scene, ray);
    
	// if not hit, return background
	if (not intersection.hit) return scene->background;
    
	// accumulate color starting with ambient
	res += scene->ambient * intersection.mat->ke;

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
		
		vec3f light_color = compute_color(*scene, *light, intersection);
		if (light_color == zero3f) continue;

		vec3f light_direction = compute_light_direction(*scene, *light, intersection);
		
		//ray3f sr = transform_ray(scene->lights->f, light_shadow_ray(l, transform_point_inverse(scene->lights->f, intersection.f.o)));
		//vec3f cl = sc * material_brdf(intersection.mat, intersection.f, sr.d, -ray.d) * max(0.0f, dot(intersection.f.z, sr.d));
		
		//vec3f cl = light_color * material_brdf() * max(0.0f, dot(intersection.norm, a.d));

		vec3f cl = intersection.mat->kd * dot(intersection.norm, light_direction);

		/*
		message("light direction: x: %f y: %f z: %f\n", light_direction.x, light_direction.y, light_direction.z);
		message("light color: x: %f y: %f z: %f\n", light_color.x, light_color.y, light_color.z);
		message("intersection in light frame: x: %f y: %f z: %f\n", intersection.norm.x, intersection.norm.y, intersection.norm.z);
		message("cl: x: %f y: %f z: %f\n", cl.x, cl.y, cl.z);
		*/

		//if (cl == zero3f) continue;
		//if (opts.shadows) {
		//	if (not intersect_scene_any(scene, sr)) c += cl;
		//}*/
		//else c += cl;
		res += cl;
	}
    
        // compute light response
        // compute light direction
        // compute the material response (brdf*cos)
        // check for shadows and accumulate if needed
    
    // if the material has reflections
        // create the reflection ray
        // accumulate the reflected light (recursive call) scaled by the material reflection

    // return the accumulated color∫
    //return zero3f;
	return res;
}

// raytrace an image
image3f raytrace(Scene* scene) {
    // allocate an image of the proper size
    auto image = image3f(scene->image_width, scene->image_height);
    
    // if no anti-aliasing
	if (true) {  
		// foreach pixel
		for (int y = 0; y < scene->image_height ; y++) {
			for (int x = 0; x < scene->image_width; x++) {
				// compute ray-camera parameters (u,v) for the pixel
				
				//what is this
				float u = (x + 0.5f) / scene->image_width;
				float v = (y + 0.5f) / scene->image_height;

				// compute camera ray
				vec2f uv = vec2f(u, v);
				Camera* camera = scene->camera;
				//WHAT IS THIS
				auto q = vec3f((uv.x - 0.5f)*scene->camera->width, (uv.y - 0.5f)*scene->camera->height, -scene->camera->dist);
				ray3f ray = ray3f(zero3f, normalize(q));
				// set pixel to the color raytraced with the ray
				//ok...
				image.at(x, scene->image_height - 1 - y) = raytrace_ray(scene, ray);		
			}
		}	
	}
        
    // else
        // foreach pixel
                // init accumulated color
                // foreach sample
                        // compute ray-camera parameters (u,v) for the pixel and the sample
                        // compute camera ray
                        // set pixel to the color raytraced with the ray
                // scale by the number of samples
    

    // done
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

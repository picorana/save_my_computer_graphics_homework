#include "animation.h"
#include "tesselation.h"

// compute the frame from an animation
frame3f animate_compute_frame(FrameAnimation* animation, int time) {
    // grab keyframe interval
    auto interval = 0;
    for(auto t : animation->keytimes) if(time < t) break; else interval++;
    interval--;
    // get translation and rotation matrices
    auto t = float(time-animation->keytimes[interval])/float(animation->keytimes[interval+1]-animation->keytimes[interval]);
    auto m_t = translation_matrix(animation->translation[interval]*(1-t)+animation->translation[interval+1]*t);
    auto m_rz = rotation_matrix(animation->rotation[interval].z*(1-t)+animation->rotation[interval+1].z*t,z3f);
    auto m_ry = rotation_matrix(animation->rotation[interval].y*(1-t)+animation->rotation[interval+1].y*t,y3f);
    auto m_rx = rotation_matrix(animation->rotation[interval].x*(1-t)+animation->rotation[interval+1].x*t,x3f);
    // compute combined xform matrix
    auto m = m_t * m_rz * m_ry * m_rx;
    // return the transformed frame
    return transform_frame(m, animation->rest_frame);
}

// update mesh frames for animation
void animate_frame(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
	for (auto mesh : scene->meshes){
		if (mesh->animation==nullptr) continue;
		mesh->frame = animate_compute_frame(mesh->animation, scene->animation->time);
	}
        // if not animation, continue
        // update frame
    // foreach surface
	for (auto surface : scene->surfaces){
		if (surface->animation==nullptr) continue;
		surface->frame = animate_compute_frame(surface->animation, scene->animation->time);
		surface->_display_mesh = make_surface_mesh(surface->frame, surface->radius, surface->isquad, surface->mat);
	}
        // if not animation, continue
        // update frame
        // update the _display_mesh
}

// skinning scene
void animate_skin(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
	for (auto mesh : scene->meshes){
		if (mesh->skinning == nullptr) continue;
		for (int i = 0; i < mesh->pos.size(); i++){
			mesh->pos[i] = zero3f;
			mesh->norm[i] = zero3f;
			for (int j = 0; j < 4; j++){
				auto w = mesh->skinning->bone_weights[i][j];
				auto index = mesh->skinning->bone_ids[i][j];
				if (index < 0) continue;
				auto xform = mat4f();
				xform = mesh->skinning->bone_xforms[scene->animation->time][index];
				auto restpo = vec4f(mesh->skinning->rest_pos[i].x, mesh->skinning->rest_pos[i].y, mesh->skinning->rest_pos[i].z, 1);
				auto res = w * xform * restpo;
				auto restno = vec4f(mesh->skinning->rest_norm[i].x, mesh->skinning->rest_norm[i].y, mesh->skinning->rest_norm[i].z, 1);
				auto resn = w * xform * restno;
				mesh->pos[i] += vec3f(res.x, res.y, res.z);
				mesh->norm[i] += vec3f(resn.x, resn.y, resn.z);
			}
			mesh->norm[i] = normalize(mesh->norm[i]);
		}
	}
        // if no skinning, continue
        // foreach vertex index
            // set pos/norm to zero
            // for each bone slot (0..3)
                // get bone weight and index
                // if index < 0, continue
                // grab bone xform
                // update position and normal
            // normalize normal
}

// particle simulation
void simulate(Scene* scene) {
    // YOUR CODE GOES HERE ---------------------
    // for each mesh
	for (auto mesh : scene->meshes){
        // skip if no simulation
		if (mesh->simulation == nullptr) continue;
        // compute time per step
		auto timenow = scene->animation->dt/scene->animation->simsteps;
		//message("%d\n", timenow);
        // foreach simulation steps
		for (int i = 0; i < scene->animation->simsteps; i++){
            // compute extenal forces (gravity)

			for (int j = 0; j < mesh->pos.size(); j++){
				mesh->simulation->force[j] = scene->animation->gravity * mesh->simulation->mass[j];
			}
			
            // for each spring, compute spring force on points
			for (auto spring : mesh->simulation->springs){
				// compute spring distance and length
				
				auto spring_length = length(mesh->pos[spring.ids.y] - mesh->pos[spring.ids.x]);
				auto spring_direction = normalize(mesh->pos[spring.ids.y] - mesh->pos[spring.ids.x]);
				auto rest_length = spring.restlength;
				auto spring_relative_vel = mesh->simulation->vel[spring.ids.y] - mesh->simulation->vel[spring.ids.x];
				// compute static force
				auto static_force_on_pi = spring_direction * (spring_length - rest_length) * spring.ks;
				// accumulate static force on points
				mesh->simulation->force[spring.ids.x] += static_force_on_pi;
				mesh->simulation->force[spring.ids.y] += - static_force_on_pi;
				// compute dynamic force
				auto dynamic_force = spring.kd * spring_direction * dot(spring_relative_vel, spring_direction);
				// accumulate dynamic force on points
				mesh->simulation->force[spring.ids.x] += dynamic_force;
				mesh->simulation->force[spring.ids.y] += - dynamic_force;
			}

            // newton laws
			for (int j = 0; j < mesh->pos.size(); j++){
                // if pinned, skip
				if (mesh->simulation->pinned[j]) continue;
                // acceleration
				auto acceleration = mesh->simulation->force[j] / mesh->simulation->mass[j];
                // update velocity and positions using Euler's method
				mesh->pos[j] += mesh->simulation->vel[j] * timenow + acceleration * timenow * timenow / 2;
				mesh->simulation->vel[j] += acceleration * timenow;
                // for each mesh, check for collision
				for (auto surf : scene->surfaces){
					// compute inside tests
                    // if quad
					if (surf->isquad){
                        // compute local poisition
						auto pl = transform_point_inverse(surf->frame, mesh->pos[j]);
						auto r = surf->radius;
                        // perform inside test
						if (pl.z < 0 && -r < pl.x && pl.x < r && -r < pl.y && pl.y < r){
							mesh->pos[j] = transform_point(surf->frame, vec3f(pl.x, pl.y, 0));
							mesh->simulation->vel[j] = (mesh->simulation->vel[j] - dot(surf->frame.z, mesh->simulation->vel[j]) * surf->frame.z)*(1 - scene->animation->bounce_dump[0]) +
								(- dot(surf->frame.z, mesh->simulation->vel[j])*surf->frame.z)*(1 - scene->animation->bounce_dump[1]);
						}
                            // if inside, set position and normal
					}
					else {
                        // inside test
						auto sphere_center = surf->frame.o;
						auto sphere_radius = surf->radius;
						if (length(mesh->pos[j] - sphere_center) < sphere_radius){
							mesh->pos[j] = sphere_radius * normalize(mesh->pos[j] - sphere_center) + sphere_center;
							auto sphere_norm = normalize(mesh->pos[j] - sphere_center);
							mesh->simulation->vel[j] = (mesh->simulation->vel[j] - dot(sphere_norm, mesh->simulation->vel[j]) * sphere_norm)*(1 - scene->animation->bounce_dump[0]) +
								(-dot(sphere_norm, mesh->simulation->vel[j])*sphere_norm)*(1 - scene->animation->bounce_dump[1]);
						}

					}
                        
                    // if inside
                        // set particle position
                        // update velocity
				}
                    
			}
        // smooth normals if it has triangles or quads
			if (mesh->quad.size() != 0 || mesh->triangle.size() != 0){
				smooth_normals(mesh);
			}

		}
	}
}

// scene reset
void animate_reset(Scene* scene) {
    scene->animation->time = 0;
    for(auto mesh : scene->meshes) {
        if(mesh->animation) {
            mesh->frame = mesh->animation->rest_frame;
        }
        if(mesh->skinning) {
            mesh->pos = mesh->skinning->rest_pos;
            mesh->norm = mesh->skinning->rest_norm;
        }
        if(mesh->simulation) {
            mesh->pos = mesh->simulation->init_pos;
            mesh->simulation->vel = mesh->simulation->init_vel;
            mesh->simulation->force.resize(mesh->simulation->init_pos.size());
        }
    }
}

// scene update
void animate_update(Scene* scene) {
    if(scene->animation->time >= scene->animation->length-1) {
        if(scene->animation->loop) animate_reset(scene);
        else return;
    } else scene->animation->time ++;
    animate_frame(scene);
    if(!scene->animation->gpu_skinning) animate_skin(scene);
    simulate(scene);
}

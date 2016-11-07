#include "tesselation.h"

// make normals for each face - duplicates all vertex data
void facet_normals(Mesh* mesh) {
    // allocates new arrays
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();
    auto triangle = vector<vec3i>();
    auto quad = vector<vec4i>();
    // froeach triangle
    for(auto f : mesh->triangle) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));
        // add triangle
        triangle.push_back({nv,nv+1,nv+2});
        // add vertex data
        for(auto i : range(3)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // froeach quad
    for(auto f : mesh->quad) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));
        // add quad
        quad.push_back({nv,nv+1,nv+2,nv+3});
        // add vertex data
        for(auto i : range(4)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // set back mesh data
    mesh->pos = pos;
    mesh->norm = norm;
    mesh->texcoord = texcoord;
    mesh->triangle = triangle;
    mesh->quad = quad;
}

// smooth out normal - does not duplicate data
void smooth_normals(Mesh* mesh) {
    // PLACEHOLDER CODE - REMOVE AFTER FUNCTION IS IMPLEMENTED
    //mesh->norm.resize(mesh->pos.size());
    // YOUR CODE GOES HERE ---------------------
    // set normals array to the same length as pos and init all elements to zero
	auto norm = vector<vec3f>();
	for (auto v : mesh->norm){
		norm.push_back(vec3f(0, 0, 0));
	}
    // foreach triangle
        // compute face normal
        // accumulate face normal to the vertex normals of each face index
    // foreach quad
	for (auto q : mesh->quad){
		mesh->norm.resize(mesh->pos.size());
		// compute face normal
		for (auto vid : q){
			mesh->norm[vid] += quad_normal(mesh->pos[f.x], mesh->pos[f.y], mesh->pos[f.z], mesh->pos[f.w]);
		}
		// accumulate face normal to the vertex normals of each face index
	}
    // normalize all vertex normals
	for (int i = 0; i < norm.size(); i++){
		norm[i] = normalize(norm[i]);
	}
	mesh->norm = norm;
}

// smooth out tangents
void smooth_tangents(Mesh* polyline) {
    // set tangent array
    polyline->norm = vector<vec3f>(polyline->pos.size(),zero3f);
    // foreach line
    for(auto l : polyline->line) {
        // compute line tangent
        auto lt = normalize(polyline->pos[l.y]-polyline->pos[l.x]);
        // accumulate segment tangent to vertex tangent on each vertex
        for (auto i : range(2)) polyline->norm[l[i]] += lt;
    }
    // normalize all vertex tangents
    for (auto& t : polyline->norm) t = normalize(t);
}

// apply Catmull-Clark mesh subdivision
// does not subdivide texcoord
void subdivide_catmullclark(Mesh* subdiv) {
	
	// YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working Mesh copied from the subdiv
	auto tesselation = subdiv;
	message("\n\ninit number of verts: %d", sizeof(subdiv->pos));

    // foreach level
	for (int i = 0; i < subdiv->subdivision_catmullclark_level; i++){
		// make empty pos and quad arrays
		auto pos = vector<vec3f>();
		auto quad = vector<vec4i>();
		// create edge_map from current mesh
		auto edge_map = EdgeMap(vector<vec3i>(), tesselation->quad);
		// linear subdivision - create vertices
		// copy all vertices from the current mesh
		for (auto vert : tesselation->pos){
			pos.push_back(vert);
		}

		int evo = pos.size();
		// add vertices in the middle of each edge (use EdgeMap)
		for (auto edge : edge_map._edge_list){
			pos.push_back(tesselation->pos[edge.x] * 0.5 + tesselation->pos[edge.y] * 0.5);
			//WHAT IS HAPPENING HERE
		}
		// add vertices in the middle of each triangle
		// add vertices in the middle of each quad
		int fvo = pos.size();
		for (auto vert : tesselation->quad){
			pos.push_back(tesselation->pos[vert.x] * 0.25 + tesselation->pos[vert.y] * 0.25 + 
				tesselation->pos[vert.z] * 0.25 + tesselation->pos[vert.w] * 0.25);
		}
		// subdivision pass --------------------------------
		// compute an offset for the edge vertices
		// compute an offset for the triangle vertices
		// compute an offset for the quad vertices
		// foreach triangle
		// add three quads to the new quad array
		// foreach quad
		/*for (auto q : tesselation->quad){
			quad.push_back(q);
		}*/
		
		for (int fid = 0; fid < tesselation->quad.size(); fid++) {
			auto f = tesselation->quad[fid];
			auto ve = vec4i(edge_map.edge_index(vec2i(f.x, f.y)),
				edge_map.edge_index(vec2i(f.y, f.z)),
				edge_map.edge_index(vec2i(f.z, f.w)),
				edge_map.edge_index(vec2i(f.w, f.x))) + vec4i(evo, evo, evo, evo);
			auto vf = fid + fvo;
			quad.push_back(vec4i(f.x, ve.x, vf, ve.w));
			quad.push_back(vec4i(f.y, ve.y, vf, ve.x));
			quad.push_back(vec4i(f.z, ve.z, vf, ve.y));
			quad.push_back(vec4i(f.w, ve.w, vf, ve.z));
		}

		// add four quads to the new quad array
		// averaging pass ----------------------------------
		// create arrays to compute pos averages (avg_pos, avg_count)
		auto avg_pos = vector<vec3f>();
		auto avg_count = vector<int>();
		// arrays have the same length as the new pos array, and are init to zero
		for (auto v : pos){
			avg_pos.push_back(vec3f(0, 0, 0));
			avg_count.push_back(0);
		}

		auto npos = vector<vec3f>(pos.size(), zero3f);
		auto count = vector<int>(pos.size(), 0);
		for (int i = 0; i < quad.size(); i++) {
			auto f = quad[i];
			for (int j = 0; j < 4; j++){
				auto vid = f[j];
				npos[vid] += (pos[f.x] + pos[f.y] + pos[f.z] + pos[f.w]) / 4;
				count[vid] ++;
			}
		}
		// for each new quad
		// compute quad center using the new pos array
		// foreach vertex index in the quad
		// normalize avg_pos with its count avg_count
		for (auto i : range(pos.size())) {
			npos[i] /= count[i];
		}
		// correction pass ----------------------------------
		// foreach pos, compute correction p = p + (avg_p - p) * (4/avg_count)
		for (auto i : range(pos.size())) {
			npos[i] = pos[i] + (npos[i] - pos[i])*(4.0 / count[i]);
		}
		// set new arrays pos, quad back into the working mesh; clear triangle array

		tesselation->pos = npos;
		tesselation->quad = quad;
	}
        
    // clear subdivision
    // according to smooth, either smooth_normals or facet_normals
	if (subdiv->subdivision_catmullclark_smooth) smooth_normals(subdiv);
	//facet_normals(subdiv);
    // copy back
	subdiv = tesselation;
	

    // clear
}

// subdivide bezier spline into line segments (assume bezier has only bezier segments and no lines)
void subdivide_bezier(Mesh* bezier) {
    // YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working polyline from bezier
    // foreach level
        // make new arrays of positions and bezier segments
        // copy all the vertices into the new array (this waste space but it is easier for now)
        // foreach bezier segment
            // apply subdivision algorithm
            // prepare indices for two new segments
            // add mid point
            // add points for first segment and fix segment indices
            // add points for second segment and fix segment indices
            // add indices for both segments into new segments array
        // set new arrays pos, segments into the working lineset
    // copy bezier segments into line segments
    // clear bezier array from lines
    // run smoothing to get proper tangents
    // copy back
    // clear
}

Mesh* make_surface_mesh(frame3f frame, float radius, bool isquad, Material* mat, float offset) {
    auto mesh = new Mesh{};
    mesh->frame = frame;
    mesh->mat = mat;
    if(isquad) {
        mesh->pos = { {-radius,-radius,-offset}, {radius,-radius,-offset},
            {radius,radius,-offset}, {-radius,radius,-offset} };
        mesh->norm = {z3f,z3f,z3f,z3f};
        mesh->quad = { {0,1,2,3} };
    } else {
        map<pair<int,int>,int> vid;
        for(auto j : range(64+1)) {
            for(auto i : range(128+1)) {
                auto u = 2 * pif * i / 64.0f, v = pif * j / 32.0f;
                auto d = vec3f{cos(u)*sin(v),sin(u)*sin(v),cos(v)};
                vid[{i,j}] = mesh->pos.size();
                mesh->pos.push_back(d*radius*(1-offset));
                mesh->norm.push_back(d);
            }
        }
        for(auto j : range(64)) {
            for(auto i : range(128)) {
                mesh->quad.push_back({vid[{i,j}],vid[{i+1,j}],vid[{i+1,j+1}],vid[{i,j+1}]});
            }
        }
    }
    return mesh;
}

void subdivide_surface(Surface* surface) {
    surface->_display_mesh = make_surface_mesh(
        surface->frame, surface->radius, surface->isquad, surface->mat);
}

void subdivide(Scene* scene) {
    for(auto mesh : scene->meshes) {
        if(mesh->subdivision_catmullclark_level) subdivide_catmullclark(mesh);
        if(mesh->subdivision_bezier_level) subdivide_bezier(mesh);
    }
    for(auto surface : scene->surfaces) {
        subdivide_surface(surface);
    }
}

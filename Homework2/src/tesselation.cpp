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
            if(!mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
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
            if(!mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
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

    // set normals array to the same length as pos and init all elements to zero
	auto norm = vector<vec3f>(mesh->pos.size());
	for (auto v : norm){
		norm.push_back(vec3f(0, 0, 0));
	}

    // foreach triangle
	for (auto f : mesh->triangle){
		// compute face normal
		auto a = normalize(cross(mesh->pos[f.y] - mesh->pos[f.x], mesh->pos[f.z] - mesh->pos[f.x]));
		// accumulate face normal to the vertex normals of each face index
		for (auto i : range(3)) {
			norm[f[i]] += a;
		}
	}
        
	// foreach quad
	for (auto f : mesh->quad){
		
		// compute face normal
		auto a = normalize(cross(mesh->pos[f.y] - mesh->pos[f.x], mesh->pos[f.z] - mesh->pos[f.x]));
		auto b = normalize(cross(mesh->pos[f.z] - mesh->pos[f.x], mesh->pos[f.w] - mesh->pos[f.x]));
		auto c = normalize(a + b);
		
		// accumulate face normal to the vertex normals of each face index
		for (auto i : range(4)) {
			norm[f[i]] += c;
		}
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
		auto edge_map = EdgeMap(tesselation->triangle, tesselation->quad);
		auto hash = EdgeMap(vector<vec3i>(), quad);

		// linear subdivision - create vertices
		// copy all vertices from the current mesh
		for (auto vert : tesselation->pos){
			pos.push_back(vert);
		}
		
		// add vertices in the middle of each edge (use EdgeMap)
		int prevpos = pos.size();
		for (int j = 0; j < edge_map._edge_map.size(); j++){
			auto edge = edge_map._edge_list[j];
			if (hash.edge_index(vec2i(edge.x, edge.y)) != -1) {
				hash._add_edge(edge.x, edge.y);
				pos.push_back(tesselation->pos[edge.x] * 0.5 + tesselation->pos[edge.y] * 0.5);
				//hash._add_edge(edge.x, prevpos + j);
				//hash._add_edge(prevpos + j, edge.y);
			}
		}

		// add vertices in the middle of each triangle
		int triangleindex = 0;
		for (auto triangle : tesselation->triangle){
			auto centroid = (pos[triangle.x] + pos[triangle.y] + pos[triangle.z]) / 3;
			pos.push_back(centroid);
		}
		// add vertices in the middle of each quad
		int fvo = pos.size();
		for (auto vert : tesselation->quad){
			pos.push_back(tesselation->pos[vert.x] * 0.25 + tesselation->pos[vert.y] * 0.25 + 
				tesselation->pos[vert.z] * 0.25 + tesselation->pos[vert.w] * 0.25);
		}

		// subdivision pass --------------------------------
		// compute an offset for the edge vertices
		int edge_offset = tesselation->pos.size();
		// compute an offset for the triangle vertices
		int triangle_offset = tesselation->pos.size();
		// compute an offset for the quad vertices
		int quad_offset = tesselation->pos.size() + tesselation->triangle.size() + tesselation->quad.size();
		
		// foreach triangle
		for (int j = 0; j < tesselation->triangle.size(); j++){
			auto t = tesselation->triangle[j];
			message("\ntriangle %d:\ncentroid: %d coords: %f %f %f\nvert0: %d coords: %f %f %f\nvert1: %d coords: %f %f %f\nvert2: %d coords: %f %f %f",
				j, triangle_offset + j, pos[triangle_offset+j].x, pos[triangle_offset+j].y, pos[triangle_offset+j].z, 
				t[0], pos[t[0]].x, pos[t[0]].y, pos[t[0]].z,
				t[1], pos[t[1]].x, pos[t[1]].y, pos[t[1]].z,
				t[2], pos[t[2]].x, pos[t[2]].y, pos[t[2]].z);
			quad.push_back(
				vec4i(t[0],
				edge_offset + edge_map.edge_index(vec2i(t[0], t[1])),
				triangle_offset + j,
				edge_offset + edge_map.edge_index(vec2i(t[0], t[2])))
				);

			quad.push_back(
				vec4i(t[1],
				edge_offset + edge_map.edge_index(vec2i(t[1], t[2])),
				triangle_offset + j,
				edge_offset + edge_map.edge_index(vec2i(t[1], t[0])))
				);

			quad.push_back(
				vec4i(t[2],
				edge_offset + edge_map.edge_index(vec2i(t[2], t[0])),
				triangle_offset + j,
				edge_offset + edge_map.edge_index(vec2i(t[2], t[1])))
				);
		}
		
		// foreach quad
		for (int fid = 0; fid < tesselation->quad.size(); fid++) {
			auto f = tesselation->quad[fid];
			auto ve = vec4i(edge_map.edge_index(vec2i(f.x, f.y)),
				edge_map.edge_index(vec2i(f.y, f.z)),
				edge_map.edge_index(vec2i(f.z, f.w)),
				edge_map.edge_index(vec2i(f.w, f.x))) + vec4i(edge_offset, edge_offset, edge_offset, edge_offset);
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
	else facet_normals(subdiv);
    // copy back
	subdiv = tesselation;
	

    // clear
}

// subdivide bezier spline into line segments (assume bezier has only bezier segments and no lines)
void subdivide_bezier(Mesh* bezier) {
    // YOUR CODE GOES HERE ---------------------
    // skip is needed
    // allocate a working polyline from bezier
	auto polyline = bezier;
    // foreach level
	for (auto i : range(polyline->subdivision_bezier_level)){
        // make new arrays of positions and bezier segments
		auto pos = vector<vec3f>();
		auto segments = vector<vec4i>();

        // copy all the vertices into the new array (this waste space but it is easier for now)
		for (auto p : polyline->pos) { 
			pos.push_back(p); 
		}

        // foreach bezier segment
		for (auto line : polyline->spline){
			// apply subdivision algorithm
			// prepare indices for two new segments
			// 1 : {P0, Q0, R0, S}
			// 2 : {S, R1, Q2, P3}
			//Qi = (Pi + Pi+1)/2
			//Ri = (Qi + Qi+1)/2
			//S  = (Ri + Ri+1)/2
			auto offset = pos.size();
			auto P0 = pos[line.x];
			auto P1 = pos[line.y];
			auto P2 = pos[line.z];
			auto P3 = pos[line.w];

			// add mid point
			auto Q0 = (P0 + P1) / 2;
			auto Q1 = (P1 + P2) / 2;
			auto Q2 = (P2 + P3) / 2;

			auto R0 = (Q0 + Q1) / 2;
			auto R1 = (Q1 + Q2) / 2;

			auto S = (R0 + R1) / 2;
			// add points for first segment and fix segment indices
			pos.push_back(Q0);
			pos.push_back(R0);
			pos.push_back(S);
			// add points for second segment and fix segment indices
			pos.push_back(R1);
			pos.push_back(Q2);
			// add indices for both segments into new segments array
			segments.push_back(vec4i(line.x, offset, offset+1, offset+2));
			segments.push_back(vec4i(offset + 2, offset+3, offset+4, line.w));
		}

        // set new arrays pos, segments into the working lineset
		polyline->pos = pos;
		polyline->spline = segments;
	}
    // copy bezier segments into line segments
	for (auto line : polyline->spline){
		polyline->line.push_back(vec2i(line.x, line.y));
		polyline->line.push_back(vec2i(line.z, line.w));
	}
    // clear bezier array from lines
    // run smoothing to get proper tangents
    // copy back
	bezier = polyline;
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

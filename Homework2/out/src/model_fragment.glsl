#version 120

varying vec3 pos;                   // [from vertex shader] position in world space
varying vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
varying vec2 texcoord;              // [from vertex shader] texture coordinate

uniform vec3 camera_pos;            // camera position (center of the camera frame)

uniform vec3 ambient;               // scene ambient

uniform int lights_num;             // number of lights
uniform vec3 light_pos[16];         // light positions
uniform vec3 light_intensity[16];   // light intensities

uniform vec3 material_kd;           // material kd
uniform vec3 material_ks;           // material ks
uniform float material_n;           // material n
uniform bool material_is_lines;     // whether the material is lines or meshes

uniform bool material_kd_txt_on;    // material kd texture enabled
uniform sampler2D material_kd_txt;  // material kd texture
uniform bool material_ks_txt_on;    // material ks texture enabled
uniform sampler2D material_ks_txt;  // material ks texture
uniform bool material_norm_txt_on;    // material norm texture enabled
uniform sampler2D material_norm_txt;  // material norm texture

float dot(vec3 a, vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
float lengthSqr(vec3 a) { return dot(a,a); }
float length(vec3 a) { return sqrt(dot(a,a)); }
float dist(vec3 a, vec3 b) { return length(a-b); }
float distSqr(vec3 a, vec3 b) { return lengthSqr(a-b); }

// main
void main() {
    // re-normalize normals
    vec3 n = normalize(norm);
    vec3 c = vec3(0,0,0);   // initialize to red to see it well
    // YOUR CODE GOES HERE ---------------------
    // lookup normal map if needed
    // compute material values by looking up textures is necessary
    // accumulate ambient
    vec3 kd = material_kd;
    vec3 ks = material_ks;
    vec3 mn = norm;
    if (material_kd_txt_on) kd = material_kd * texture2D(material_kd_txt, texcoord).rgb;
    if (material_ks_txt_on) ks = material_ks * texture2D(material_ks_txt, texcoord).rgb;
    if (material_norm_txt_on) mn = normalize(texture2D(material_norm_txt, texcoord).rgb*2 - 1);
	c = c + ambient*kd;
    // foreach light
	for (int i=0; i<lights_num; i++){
        // compute point light color at pos
        // compute light direction at pos
        // compute view direction using camera_pos and pos
        // compute h
        // accumulate blinn-phong model
        vec3 S = light_pos[i];
        vec3 P = pos;
        vec3 l = normalize(S - P);
        vec3 light_color = light_intensity[i] / distSqr(S, P);
        vec3 light_direction = normalize((S - P)/abs(dist(S, P)));
        vec3 v = normalize(camera_pos - P);
        vec3 h = (l + v) / length(l + v);
        float dotres = abs(dot(norm, h));

        vec3 cl = light_color * (kd + ks * pow(max(0.0, dotres), material_n)) * abs(dot(l, mn));
        c = c + cl;
	}
    // output final color by setting gl_FragColor
    gl_FragColor = vec4(c,1);
}

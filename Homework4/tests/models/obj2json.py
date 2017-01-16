#! /usr/bin/env python -B

import argparse, json, textwrap
from math import sqrt

def _flatten(v):
    if not v: return []
    return [ vvv for vv in v for vvv in vv ]

def obj2json(filename_obj,filename_json,
             material_lowercase=False,material_override='',
             smooth_normals=False,
             flipyz=False,flipy=False,mesh_only=False):
    # parse obj
    m = 'default'
    n = 'default'
    v, vt, vn = [], [], []
    mesh = { 'name': n, 'pos': [], 'norm': [], 'texcoord': [], 'triangle': [], 'quad': [], 'material': m }
    meshes = [ mesh ]
    if not material_override:
        materials = { 'default': { 'ke': [0,0,0], 'kd': [1,1,1], 'ks': [0,0,0], 'n': 100, 'ke_txt': '', 'kd_txt': '', 'ks_txt': '', 'norm_txt': '' } }
    else:
        with open(material_override) as f: material = json.load(f)
        materials = { 'default': material }
    vids_map = {}
    for line in open(filename_obj):
        if not line or line.startswith('#') or line.isspace(): continue
        tokens = line.split()
        token, values = tokens[0], tokens[1:]
        if token == 'v': v += [ [float(t) for t in values] ]
        elif token == 'vt': vt += [ [float(t) for t in values] ]
        elif token == 'vn': vn += [ [float(t) for t in values] ]
        elif token == 'f':
            f = []
            for vids_t in values:
                if '/' not in vids_t: vids_t += '//'
                elif vids_t.count('/') == 1: vids_t += '/'
                vids = [int(vid) if vid else 0 for vid in vids_t.split('/')]
                if vids[0] < 0: vids[0] = len(v)+vids[0]+1
                if vids[1] < 0: vids[1] = len(vt)+vids[1]+1
                if vids[2] < 0: vids[2] = len(vn)+vids[2]+1
                if tuple(vids) not in vids_map:
                    mesh['pos'] += [ v[vids[0]-1] ]
                    if vids[1]: mesh['texcoord'] += [ vt[vids[1]-1] ]
                    if vids[2]: mesh['norm'] += [ vn[vids[2]-1] ]
                    vids_map[tuple(vids)] = len(mesh['pos'])-1
                f += [ vids_map[tuple(vids)] ]
            if len(f) == 3: mesh['triangle'] += [ f ]
            elif len(f) == 4: mesh['quad'] += [ f ]
            else: print 'mesh {}: face length not supported {}'.format(n,len(f))
        elif token == 'o' or token == 'g':
            if token == 'o': n = values[0]
            mesh = { 'name': n, 'pos': [], 'norm': [], 'texcoord': [], 'triangle': [], 'quad': [], 'material': m }
            meshes += [ mesh ]
            vids_map = {}
        elif token == 'usemtl':
            if material_override: continue
            m = values[0]
            if material_lowercase: m = m.lower()
            mesh = { 'name': n, 'pos': [], 'norm': [], 'texcoord': [], 'triangle': [], 'quad': [], 'material': m }
            meshes += [ mesh ]
            vids_map = {}
        elif token == 'mtllib':
            if material_override: continue
            material = None
            mname = None
            for mline in open(values[0]):
                if not mline or mline.startswith('#') or mline.isspace(): continue
                mtokens = mline.split()
                mtoken, mvalues = mtokens[0], mtokens[1:]
                if mtoken == 'newmtl':
                    mname = mvalues[0]
                    material = { 'ke': [0,0,0], 'kd': [1,1,1], 'ks': [0,0,0], 'n': 100, 'ke_txt': '', 'kd_txt': '', 'ks_txt': '', 'norm_txt': '' }
                    materials[mname] = material
                elif mtoken == 'Ke': material['ke'] = [ float(t) for t in mvalues ]
                elif mtoken == 'Kd': material['kd'] = [ float(t) for t in mvalues ]
                elif mtoken == 'Ks': material['ks'] = [ float(t) for t in mvalues ]
                elif mtoken == 'Ns': material['n'] = float(mvalues[0])
                elif mtoken == 'map_Ke': material['ke_txt'] = mvalues[0]
                elif mtoken == 'map_Kd': material['kd_txt'] = mvalues[0]
                elif mtoken == 'map_Ks': material['ks_txt'] = mvalues[0]
                else: print 'material {}: ignored statement: {}'.format(mname,mline.strip())
        else: print 'mesh {}: unknown command {}'.format(n,line.strip())

    # remove empty meshes
    meshes = [ mesh for mesh in meshes if mesh['triangle'] or mesh['quad'] ]

    # flip if needed
    if flipyz:
        for mesh in meshes:
            mesh['pos']  = [ [ vv[0], vv[2], vv[1] ] for vv in mesh['pos']  ]
            mesh['norm'] = [ [ vv[0], vv[2], vv[1] ] for vv in mesh['norm'] ]
    if flipy:
        for mesh in meshes:
            mesh['pos']  = [ [ vv[0], -vv[1], vv[2] ] for vv in mesh['pos']  ]
            mesh['norm'] = [ [ vv[0], -vv[1], vv[2] ] for vv in mesh['norm'] ]

    ## smooth normals if needed
    if smooth_normals:
        for mesh in meshes:
            if mesh['norm']: continue
            mesh['norm'] = []
            for _ in mesh['pos']: mesh['norm'].append([0.0,0.0,0.0])
            for f in mesh['triangle']+mesh['quad']:
                p0, p1, p2 = mesh['pos'][f[0]], mesh['pos'][f[1]], mesh['pos'][f[2]]
                a = [ p1[i]-p0[i] for i in range(3) ]
                b = [ p2[i]-p0[i] for i in range(3) ]
                nn = [ a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0] ]
                ll = sqrt( nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2] )
                for i in range(3): nn[i] /= ll
                for vid in f:
                    for i in range(3): mesh['norm'][vid][i] += nn[i]
            for nn in mesh['norm']:
                ll = sqrt( nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2] )
                for i in range(3): nn[i] /= ll

    # clear unused materuial elements
    for material in materials.values():
        if not material['ke_txt']: del material['ke_txt']
        if not material['kd_txt']: del material['kd_txt']
        if not material['ks_txt']: del material['ks_txt']
        if not material['norm_txt']: del material['norm_txt']
        if material['ke'] == [0,0,0]: del material['ke']

    # assign materials
    for mesh in meshes:
        mesh['material'] = materials[ mesh['material'] ]

    # flatten arrays or clear them
    for mesh in meshes:
        if mesh['pos']: mesh['pos'] = [ vvv for vv in mesh['pos'] for vvv in vv ]
        if mesh['norm']: mesh['norm'] = [ vvv for vv in mesh['norm'] for vvv in vv ]
        else: del mesh['norm']
        if mesh['texcoord']: mesh['texcoord'] = [ vvv for vv in mesh['texcoord'] for vvv in vv ]
        else: del mesh['texcoord']
        if mesh['triangle']: mesh['triangle'] = [ vvv for vv in mesh['triangle'] for vvv in vv ]
        else: del mesh['triangle']
        if mesh['quad']: mesh['quad'] = [ vvv for vv in mesh['quad'] for vvv in vv ]
        else: del mesh['quad']

    # pack json objects
    if mesh_only: json_value = meshes[0]
    else: json_value = meshes 

    # save json format
    with open(filename_json,'w') as f: json.dump(json_value,f,indent=2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('-ml','--material_lowercase',default=False,action='store_true')
    parser.add_argument('-mo','--material_override',default='')
    parser.add_argument('-s','--smooth_normals',default=False,action='store_true')
    parser.add_argument('-f','--flipyz',default=False,action='store_true')
    parser.add_argument('-y','--flipy',default=False,action='store_true')
    parser.add_argument('-m','--mesh_only',default=False,action='store_true')
    args = parser.parse_args()
    obj2json(args.filename,args.filename.replace('.obj','.json'),
        material_lowercase = args.material_lowercase,
        material_override = args.material_override,
        smooth_normals = args.smooth_normals,
        flipyz = args.flipyz,
        flipy = args.flipy,
        mesh_only = args.mesh_only)


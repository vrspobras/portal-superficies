from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse
from pathlib import Path
import json, uuid, re, csv, math
from xml.etree import ElementTree as ET
import numpy as np
from scipy.spatial import Delaunay
from pyproj import Transformer
ROOT=Path(__file__).resolve().parent; DATA=ROOT/'data'; STATIC=ROOT/'static'; HOST,PORT='0.0.0.0',8000; ALLOWED={'SAU3','SAU4'}
TRANSFORM=Transformer.from_crs('EPSG:31982','EPSG:4326',always_xy=True)
def clean_tag(tag): return tag.rsplit('}',1)[-1].lower()
def fnum(v): return float(str(v).strip().replace(',','.'))
def parse_delimited(text):
    text=text.lstrip('\ufeff').replace('\r\n','\n').replace('\r','\n'); lines=[x.strip() for x in text.split('\n') if x.strip()]
    if not lines: raise ValueError('Arquivo vazio.')
    first=lines[0]; delimiter=';' if first.count(';')>max(first.count(','),first.count('\t')) else ('\t' if first.count('\t')>first.count(',') else ',')
    rows=[next(csv.reader([ln],delimiter=delimiter)) for ln in lines]; header=[re.sub(r'[^a-z0-9]+','',c.strip().lower()) for c in rows[0]]
    if all(k in header for k in ('n','e','z')):
        i_n,i_e,i_z=header.index('n'),header.index('e'),header.index('z'); i_name=header.index('name') if 'name' in header else None; pts=[]
        for row in rows[1:]:
            if len(row)<=max(i_n,i_e,i_z): continue
            name=row[i_name].strip() if i_name is not None and i_name<len(row) else str(len(pts)+1)
            if name.lower().startswith('base:'): continue
            try:n,e,z=fnum(row[i_n]),fnum(row[i_e]),fnum(row[i_z])
            except:continue
            if all(np.isfinite(v) for v in (n,e,z)): pts.append({'id':name,'description':name,'north':n,'east':e,'z':z})
        if len(pts)>=3:return pts
    cand={'id':['numerodoponto','numero','ponto','name','nome'],'desc':['descricao','description','desc'],'north':['norte','north','n'],'east':['este','east','e'],'z':['cota','elevacao','elevation','z']};idx={}
    for key,names in cand.items():
        for name in names:
            if name in header:idx[key]=header.index(name);break
    if all(k in idx for k in ('id','north','east','z')):
        pts=[]
        for row in rows[1:]:
            if len(row)<=max(idx.values()):continue
            try:n,e,z=fnum(row[idx['north']]),fnum(row[idx['east']]),fnum(row[idx['z']])
            except:continue
            if not all(np.isfinite(v) for v in (n,e,z)):continue
            ident=row[idx['id']].strip() or str(len(pts)+1);desc=row[idx['desc']].strip() if 'desc' in idx and idx['desc']<len(row) else ident
            pts.append({'id':ident,'description':desc,'north':n,'east':e,'z':z})
        if len(pts)>=3:return pts
    pts=[]
    for row in rows:
        if len(row)<5:continue
        try:float(row[0]);ident=row[1].strip();n,e,z=fnum(row[2]),fnum(row[3]),fnum(row[4])
        except:continue
        if all(np.isfinite(v) for v in (n,e,z)):pts.append({'id':ident or str(len(pts)+1),'description':row[-1].strip() if row else ident,'north':n,'east':e,'z':z})
    if len(pts)>=3:return pts
    raise ValueError('Não foi possível identificar pelo menos 3 pontos com Norte, Este e Cota.')
def parse_landxml(text):
    try:root=ET.fromstring(text)
    except ET.ParseError as e:raise ValueError(f'XML inválido: {e}')
    points=[];id_to_idx={};faces_raw=[]
    for elem in root.iter():
        tag=clean_tag(elem.tag)
        if tag in ('cgpoint','p'):
            vals=re.split(r'[\s,;]+',(elem.text or '').strip())
            if len(vals)<3:continue
            try:a,b,z=fnum(vals[0]),fnum(vals[1]),fnum(vals[2])
            except:continue
            ident=str(elem.attrib.get('name') or elem.attrib.get('id') or len(points)+1);desc=str(elem.attrib.get('desc') or elem.attrib.get('description') or ident)
            if ident.lower().startswith('base') or ident in id_to_idx:continue
            id_to_idx[ident]=len(points);points.append({'id':ident,'description':desc,'north':b,'east':a,'z':z})
        elif tag=='face':
            vals=re.findall(r'-?\d+',(elem.text or ''))
            if len(vals)>=3:faces_raw.append(vals[:3])
    faces=[]
    for f in faces_raw:
        try:idxs=[id_to_idx[x] for x in f]
        except KeyError:
            try:idxs=[int(x) for x in f]
            except:continue
        if len(set(idxs))==3 and min(idxs)>=0 and max(idxs)<len(points):faces.append(idxs)
    if len(points)<3:raise ValueError('O XML não contém pelo menos 3 pontos reconhecíveis.')
    return points,faces
def build_tin(points,faces_override=None):
    unique=[];seen=set()
    for p in points:
        key=(round(float(p['east']),8),round(float(p['north']),8))
        if key in seen:continue
        seen.add(key);unique.append(p)
    if len(unique)<3:raise ValueError('São necessários pelo menos 3 pontos distintos.')
    xy=np.array([[p['east'],p['north']] for p in unique],float);origin=xy.mean(axis=0);local=xy-origin
    if faces_override:
        faces=[list(map(int,f)) for f in faces_override if len(f)==3 and min(f)>=0 and max(f)<len(unique)]
        if not faces:faces=Delaunay(local,qhull_options='QJ').simplices.astype(int).tolist()
    else:faces=Delaunay(local,qhull_options='QJ').simplices.astype(int).tolist()
    verts=[{'id':p['id'],'description':p.get('description',p['id']),'east':float(p['east']),'north':float(p['north']),'z':float(p['z']),'x':float(p['east']-origin[0]),'y':float(p['north']-origin[1])} for p in unique]
    return {'vertices':verts,'faces':faces,'origin':{'east':float(origin[0]),'north':float(origin[1])},'zmin':float(min(v['z'] for v in verts)),'zmax':float(max(v['z'] for v in verts)),'points':len(verts),'triangles':len(faces)}
def parse_file(filename,data):
    text=data.decode('utf-8-sig',errors='replace')
    if filename.lower().endswith(('.xml','.landxml')):pts,faces=parse_landxml(text);return pts,faces,'LandXML'
    return parse_delimited(text),None,'TXT/CSV'
def list_surfaces(sau):
    folder=DATA/sau;folder.mkdir(parents=True,exist_ok=True);out=[]
    for f in sorted(folder.glob('*.json'),key=lambda p:p.stat().st_mtime,reverse=True):
        try:d=json.loads(f.read_text());m=d['mesh'];out.append({'id':d['id'],'name':d['name'],'filename':d.get('filename',''),'created':d.get('created',''),'source':d.get('source',''),'points':m['points'],'triangles':m['triangles']})
        except:continue
    return out
def save_surface(sau,name,filename,source,pts,mesh):
    sid=uuid.uuid4().hex[:10];payload={'id':sid,'name':name or 'Nova superfície','filename':filename,'source':source,'created':__import__('datetime').datetime.now().isoformat(timespec='seconds'),'points_raw':pts,'mesh':mesh};(DATA/sau/f'{sid}.json').write_text(json.dumps(payload,ensure_ascii=False));return payload
def load_surface(sau,sid):
    if sau not in ALLOWED:raise ValueError('SAU inválida')
    p=DATA/sau/f'{sid}.json'
    if not p.exists():raise FileNotFoundError(sid)
    return json.loads(p.read_text())
def save_payload(sau,payload):(DATA/sau/f"{payload['id']}.json").write_text(json.dumps(payload,ensure_ascii=False))
def point_in_triangle(pt,a,b,c,eps=1e-9):
    px,py=pt;v0x,v0y=c[0]-a[0],c[1]-a[1];v1x,v1y=b[0]-a[0],b[1]-a[1];v2x,v2y=px-a[0],py-a[1];den=v0x*v1y-v1x*v0y
    if abs(den)<eps:return None
    u=(v2x*v1y-v1x*v2y)/den;v=(v0x*v2y-v2x*v0y)/den
    if u>=-eps and v>=-eps and u+v<=1+eps:return (1-u-v)*a[2]+v*b[2]+u*c[2]
    return None
def z_at(mesh,e,n):
    vs=mesh['vertices']
    for f in mesh['faces']:
        a,b,c=(vs[f[i]] for i in range(3));val=point_in_triangle((e,n),(a['east'],a['north'],a['z']),(b['east'],b['north'],b['z']),(c['east'],c['north'],c['z']))
        if val is not None:return float(val)
    return None
def section_profile(mesh,A,B,samples=500):
    e1,n1=A;e2,n2=B;length=math.hypot(e2-e1,n2-n1)
    if length<=1e-9:raise ValueError('A e B não podem ser iguais.')
    prof=[];samples=max(20,min(int(samples),1000))
    for i in range(samples+1):
        t=i/samples;e=e1+(e2-e1)*t;n=n1+(n2-n1)*t;z=z_at(mesh,e,n)
        if z is not None:prof.append({'distance':length*t,'east':e,'north':n,'z':z})
    return {'length':length,'profile':prof}
def geojson_surface(mesh):
    features=[];vs=mesh['vertices']
    for f in mesh['faces']:
        coords=[]
        for idx in f+[f[0]]:
            lon,lat=TRANSFORM.transform(vs[idx]['east'],vs[idx]['north']);coords.append([lon,lat])
        features.append({'type':'Feature','properties':{'zavg':sum(vs[i]['z'] for i in f)/3},'geometry':{'type':'Polygon','coordinates':[coords]}})
    for ii,v in enumerate(vs):
        lon,lat=TRANSFORM.transform(v['east'],v['north']);features.append({'type':'Feature','properties':{'index':ii,'id':v['id'],'description':v['description'],'z':v['z'],'east':v['east'],'north':v['north']},'geometry':{'type':'Point','coordinates':[lon,lat]}})
    return {'type':'FeatureCollection','features':features}

def base_path(sau):
    if sau not in ALLOWED: raise ValueError('SAU inválida')
    folder=DATA/sau; folder.mkdir(parents=True, exist_ok=True)
    return folder/'base.json'

def get_base(sau):
    p=base_path(sau)
    if not p.exists(): return None
    return json.loads(p.read_text())

def save_base(sau,name,filename,source,pts,mesh):
    payload={'id':'BASE','name':name or 'Projeto Final - Base','filename':filename,'source':source,'created':__import__('datetime').datetime.now().isoformat(timespec='seconds'),'points_raw':pts,'mesh':mesh}
    base_path(sau).write_text(json.dumps(payload,ensure_ascii=False))
    return payload

def load_mesh_for_ref(sau,ref):
    if ref in (None,'','__BASE__','BASE'):
        b=get_base(sau)
        if not b: raise ValueError('Projeto final (Base) ainda não foi cadastrado para esta SAU.')
        return b['mesh']
    return load_surface(sau,ref)['mesh']

def edge_key(a,b): return (a,b) if a<b else (b,a)


def orient2d(a,b,c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])

def segments_proper(a,b,c,d,eps=1e-9):
    o1=orient2d(a,b,c); o2=orient2d(a,b,d); o3=orient2d(c,d,a); o4=orient2d(c,d,b)
    return ((o1>eps and o2<-eps) or (o1<-eps and o2>eps)) and ((o3>eps and o4<-eps) or (o3<-eps and o4>eps))

def mesh_edges(mesh):
    out=set()
    for f in mesh['faces']:
        for i in range(3): out.add(edge_key(f[i],f[(i+1)%3]))
    return out

def recover_breakline(mesh, chain, max_iter=20000):
    chain=[int(x) for x in chain]
    if len(chain)<2: raise ValueError('Uma breakline precisa de pelo menos 2 vértices.')
    n=len(mesh['vertices'])
    if any(i<0 or i>=n for i in chain): raise ValueError('Breakline contém vértice inválido.')
    chain2=[]
    for i in chain:
        if not chain2 or chain2[-1]!=i: chain2.append(i)
    for u,v in zip(chain2,chain2[1:]):
        if u==v: continue
        A=(mesh['vertices'][u]['east'],mesh['vertices'][u]['north']); B=(mesh['vertices'][v]['east'],mesh['vertices'][v]['north'])
        for _ in range(max_iter):
            edges=mesh_edges(mesh)
            target=edge_key(u,v)
            if target in edges: break
            crossings=[]
            for a,b in edges:
                if u in (a,b) or v in (a,b): continue
                C=(mesh['vertices'][a]['east'],mesh['vertices'][a]['north']); D=(mesh['vertices'][b]['east'],mesh['vertices'][b]['north'])
                if segments_proper(A,B,C,D): crossings.append((a,b))
            if not crossings:
                raise ValueError(f'Não foi possível encaixar a breakline entre os vértices {mesh["vertices"][u]["id"]} e {mesh["vertices"][v]["id"]}.')
            flipped=False
            for a,b in crossings:
                adj=[];keep=[]
                for f in mesh['faces']:
                    if a in f and b in f: adj.append(f)
                    else: keep.append(f)
                if len(adj)!=2: continue
                c=next(x for x in adj[0] if x not in (a,b)); d=next(x for x in adj[1] if x not in (a,b))
                P=[(mesh['vertices'][i]['east'],mesh['vertices'][i]['north']) for i in (a,b,c,d)]
                if orient2d(P[0],P[1],P[2])*orient2d(P[0],P[1],P[3])>=0: continue
                if orient2d(P[2],P[3],P[0])*orient2d(P[2],P[3],P[1])>=0: continue
                newe=edge_key(c,d)
                if newe in edges and newe!=target: continue
                mesh['faces']=keep+[[c,d,a],[d,c,b]]
                flipped=True; break
            if not flipped:
                raise ValueError(f'Não foi possível recuperar a breakline entre os vértices {mesh["vertices"][u]["id"]} e {mesh["vertices"][v]["id"]}.')
        else:
            raise ValueError('Limite de iterações atingido ao recuperar a breakline.')
    mesh['triangles']=len(mesh['faces'])
    mesh['breaklines']=mesh.get('breaklines',[])+[chain2]
    return mesh

def remove_edge_from_mesh(mesh,a,b):
    a,b=int(a),int(b)
    faces=mesh['faces']
    adj=[];keep=[]
    for f in faces:
        if a in f and b in f: adj.append(f)
        else: keep.append(f)
    if not adj:
        raise ValueError('Aresta não encontrada.')
    if len(adj)==1:
        newfaces=keep
        action='boundary_hole'
    else:
        # Interior edge: replace the selected diagonal by the opposite diagonal.
        others=[]
        for f in adj:
            others.extend([x for x in f if x not in (a,b)])
        if len(others)!=2 or others[0]==others[1]:
            raise ValueError('Não foi possível identificar os triângulos adjacentes.')
        c,d=others
        newfaces=keep+[[c,d,a],[d,c,b]]
        action='flip_to_opposite_diagonal'
    out=dict(mesh);out['faces']=newfaces;out['triangles']=len(newfaces)
    return out,action



def flip_edge_in_mesh(mesh,a,b):
    a,b=int(a),int(b)
    adj=[];keep=[]
    for f in mesh['faces']:
        if a in f and b in f: adj.append(f)
        else: keep.append(f)
    if len(adj)!=2:
        raise ValueError('A diagonal interna com dois triângulos adjacentes é necessária.')
    others=[]
    for f in adj:
        others.append(next(x for x in f if x not in (a,b)))
    c,d=others
    if c==d: raise ValueError('Não foi possível identificar os vértices opostos.')
    # Ensure the new diagonal c-d replaces a-b.
    newfaces=keep+[[c,d,a],[d,c,b]]
    out=dict(mesh);out['faces']=newfaces;out['triangles']=len(newfaces)
    return out

def reset_mesh_from_points(payload):
    pts=payload.get('points_raw') or []
    return build_tin(pts)

try:
    from shapely.geometry import Polygon, Point
    from shapely.strtree import STRtree
    HAVE_SHAPELY=True
except Exception:
    HAVE_SHAPELY=False

def triangle_plane_z(v1,v2,v3,x,y):
    x1,y1,z1=v1['east'],v1['north'],v1['z']; x2,y2,z2=v2['east'],v2['north'],v2['z']; x3,y3,z3=v3['east'],v3['north'],v3['z']
    den=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)
    if abs(den)<1e-12: return None
    return z1 + ((z2-z1)*(y3-y1)-(z3-z1)*(y2-y1))*(x-x1)/den + (-(z2-z1)*(x3-x1)+(z3-z1)*(x2-x1))*(y-y1)/den

def compare_tins(mesh_a,mesh_b):
    if not HAVE_SHAPELY: raise ValueError('Shapely não está instalado no ambiente.')
    va,vb=mesh_a['vertices'],mesh_b['vertices']
    polys=[]; valid_b=[]
    for j,f in enumerate(mesh_b['faces']):
        coords=[(vb[i]['east'],vb[i]['north']) for i in f]
        poly=Polygon(coords)
        if poly.is_valid and poly.area>0:
            polys.append(poly);valid_b.append(j)
    if not polys: raise ValueError('A superfície de referência não possui triângulos válidos.')
    tree=STRtree(polys)
    cut=fill=area=wsum=0.0
    va_idx=list(range(len(va)))
    for f in mesh_a['faces']:
        ca=[(va[i]['east'],va[i]['north']) for i in f]
        pa=Polygon(ca)
        if not pa.is_valid or pa.area<=0: continue
        for hit in tree.query(pa, predicate='intersects'):
            # shapely 2 returns integer indices; support older versions returning geometry objects.
            if isinstance(hit,(int,np.integer)):
                j=int(hit); pb=polys[j]; fb=mesh_b['faces'][valid_b[j]]
            else:
                try:j=polys.index(hit)
                except ValueError: continue
                pb=hit; fb=mesh_b['faces'][valid_b[j]]
            inter=pa.intersection(pb)
            if inter.is_empty: continue
            geoms=list(getattr(inter,'geoms',[])) if inter.geom_type=='GeometryCollection' else [inter]
            for g in geoms:
                if g.geom_type not in ('Polygon','MultiPolygon') or g.area<=1e-9: continue
                parts=list(g.geoms) if g.geom_type=='MultiPolygon' else [g]
                for part in parts:
                    area_i=part.area
                    c=part.centroid
                    za=triangle_plane_z(*(va[i] for i in f),c.x,c.y)
                    zb=triangle_plane_z(*(vb[i] for i in fb),c.x,c.y)
                    if za is None or zb is None: continue
                    diff=za-zb
                    vol=area_i*diff
                    area+=area_i;wsum+=vol
                    if diff>=0: cut+=vol
                    else: fill-=vol
    return {'method':'TIN × TIN','cut':cut,'fill':fill,'net':cut-fill,'area':area,'mean_diff':(wsum/area if area else 0.0)}
BASE_CAD_FILES={'SAU4':'base_cad.json'}
def cad_path(sau):
    if sau not in ALLOWED: raise ValueError('SAU inválida')
    folder=DATA/sau; folder.mkdir(parents=True,exist_ok=True)
    return folder/BASE_CAD_FILES.get(sau,'base_cad.json')

def load_base_cad(sau):
    p=cad_path(sau)
    if not p.exists(): return None
    return json.loads(p.read_text())

def cad_geojson(cad):
    features=[]
    for g in cad.get('geometry',[]):
        pts=g.get('points',[])
        if not pts: continue
        layer=g.get('layer','0')
        typ=g.get('type')
        if typ=='POINT':
            x,y=pts[0][0],pts[0][1]; lon,lat=TRANSFORM.transform(x,y)
            features.append({'type':'Feature','properties':{'layer':layer,'cad_type':typ,'z':pts[0][2]},'geometry':{'type':'Point','coordinates':[lon,lat]}})
        elif len(pts)>=2:
            coords=[]
            for x,y,z in pts:
                lon,lat=TRANSFORM.transform(x,y); coords.append([lon,lat])
            if g.get('closed') and coords[0]!=coords[-1]: coords.append(coords[0])
            features.append({'type':'Feature','properties':{'layer':layer,'cad_type':typ},'geometry':{'type':'LineString','coordinates':coords}})
    return {'type':'FeatureCollection','features':features,'source':cad.get('source',''),'sau':cad.get('sau','')}

def segment_intersection_param(a,b,c,d,eps=1e-10):
    ax,ay=a; bx,by=b; cx,cy=c; dx,dy=d
    rx,ry=bx-ax,by-ay; sx,sy=dx-cx,dy-cy
    den=rx*sy-ry*sx
    if abs(den)<eps:return None
    qpx,qpy=cx-ax,cy-ay
    t=(qpx*sy-qpy*sx)/den; u=(qpx*ry-qpy*rx)/den
    if -eps<=t<=1+eps and -eps<=u<=1+eps:return max(0.0,min(1.0,t))
    return None

def section_cad_intersections(cad,A,B):
    if not cad:return []
    ae,an=A;be,bn=B;out=[];dx=be-ae;dy=bn-an;L=math.hypot(dx,dy)
    if L<=1e-9:return out
    for g in cad.get('geometry',[]):
        pts=g.get('points',[])
        if len(pts)<2:continue
        for i in range(len(pts)-1):
            p1,p2=pts[i],pts[i+1]
            t=segment_intersection_param((ae,an),(be,bn),(p1[0],p1[1]),(p2[0],p2[1]))
            if t is None:continue
            z=None
            z1,z2=p1[2],p2[2]
            if abs(z1)>1e-9 or abs(z2)>1e-9:
                z=z1+(z2-z1)*(t if math.hypot(p2[0]-p1[0],p2[1]-p1[1])>1e-12 else 0)
            out.append({'distance':L*t,'east':ae+dx*t,'north':an+dy*t,'z':z,'layer':g.get('layer','0'),'cad_type':g.get('type','')})
    return out

class Handler(BaseHTTPRequestHandler):
    def _send(self,code,data,ctype='application/json; charset=utf-8'):
        body=data if isinstance(data,bytes) else data.encode('utf-8');self.send_response(code);self.send_header('Content-Type',ctype);self.send_header('Content-Length',str(len(body)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(body)
    def _json(self,code,obj):self._send(code,json.dumps(obj,ensure_ascii=False))
    def do_GET(self):
        u=urlparse(self.path)
        if u.path=='/':return self._send(200,(STATIC/'index.html').read_bytes(),'text/html; charset=utf-8')
        if u.path=='/health':return self._json(200,{'ok':True})
        if u.path.startswith('/api/surfaces/'):
            sau=u.path.split('/')[3].upper();return self._json(200,list_surfaces(sau)) if sau in ALLOWED else self._json(400,{'error':'SAU inválida'})
        if u.path.startswith('/api/surface_geo/'):
            p=u.path.split('/');sau=p[3].upper();sid=p[4]
            try:return self._json(200,geojson_surface(load_surface(sau,sid)['mesh']))
            except Exception as e:return self._json(404,{'error':str(e)})
        if u.path.startswith('/api/surface/'):
            p=u.path.split('/');sau=p[3].upper();sid=p[4]
            try:return self._json(200,load_surface(sau,sid))
            except Exception as e:return self._json(404,{'error':str(e)})
        if u.path.startswith('/api/base_cad/'):
            sau=u.path.split('/')[3].upper()
            try:
                c=load_base_cad(sau)
                return self._json(200,cad_geojson(c)) if c else self._json(404,{'error':'Base CAD não cadastrada.'})
            except Exception as e:return self._json(400,{'error':str(e)})
        if u.path.startswith('/api/base/'):
            sau=u.path.split('/')[3].upper()
            try:
                b=get_base(sau)
                return self._json(200,b) if b else self._json(404,{'error':'Base não cadastrada.'})
            except Exception as e:return self._json(400,{'error':str(e)})
        return self._send(404,b'Not found','text/plain; charset=utf-8')
    def _multipart(self,body,ctype):
        m=re.search(r'boundary=(?:"([^"]+)"|([^;]+))',ctype)
        if not m:raise ValueError('Boundary ausente.')
        boundary=(m.group(1) or m.group(2)).encode();parts=body.split(b'--'+boundary);fields={};file_name=file_bytes=None
        for part in parts:
            if b'Content-Disposition' not in part:continue
            head,sep,content=part.partition(b'\r\n\r\n')
            if not sep:head,sep,content=part.partition(b'\n\n')
            if not sep:continue
            fm=re.search(br'filename="([^"]*)"',head);nm=re.search(br'name="([^"]+)"',head);content=content.rstrip(b'\r\n')
            if fm:file_name=fm.group(1).decode('utf-8','ignore');file_bytes=content
            elif nm:fields[nm.group(1).decode('utf-8','ignore')]=content.decode('utf-8','ignore')
        return fields,file_name,file_bytes
    def do_POST(self):
        u=urlparse(self.path);body=self.rfile.read(int(self.headers.get('Content-Length','0')))
        try:
            if u.path=='/api/upload':
                fields,filename,file_bytes=self._multipart(body,self.headers.get('Content-Type',''))
                if not file_bytes or not filename:raise ValueError('Nenhum arquivo recebido.')
                sau=fields.get('sau','SAU3').upper();name=fields.get('name','').strip() or Path(filename).stem
                if sau not in ALLOWED:raise ValueError('SAU inválida.')
                pts,faces,source=parse_file(filename,file_bytes);mesh=build_tin(pts,faces);return self._json(200,{'ok':True,'surface':save_surface(sau,name,filename,source,pts,mesh)})
            if u.path=='/api/base_upload':
                fields,filename,file_bytes=self._multipart(body,self.headers.get('Content-Type',''))
                if not file_bytes or not filename: raise ValueError('Nenhum arquivo recebido.')
                sau=fields.get('sau','SAU3').upper(); name=fields.get('name','').strip() or Path(filename).stem
                if sau not in ALLOWED: raise ValueError('SAU inválida.')
                pts,faces,source=parse_file(filename,file_bytes); mesh=build_tin(pts,faces); return self._json(200,{'ok':True,'base':save_base(sau,name,filename,source,pts,mesh)})
            if u.path=='/api/base_delete':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper();
                p=base_path(sau)
                if not p.exists(): raise ValueError('Base não encontrada.')
                p.unlink(); return self._json(200,{'ok':True})
            if u.path=='/api/breakline':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper(); sid=d.get('id'); chain=d.get('chain') or []
                if sau not in ALLOWED or not sid: raise ValueError('Superfície inválida.')
                payload=load_surface(sau,sid); payload['mesh']=recover_breakline(payload['mesh'],chain); save_payload(sau,payload); return self._json(200,{'ok':True,'surface':payload})
            if u.path=='/api/flip_edge':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper(); sid=d.get('id'); edge=d.get('edge') or []
                if sau not in ALLOWED or not sid or len(edge)!=2: raise ValueError('Aresta inválida.')
                payload=load_surface(sau,sid); payload['mesh']=flip_edge_in_mesh(payload['mesh'],edge[0],edge[1]); save_payload(sau,payload); return self._json(200,{'ok':True,'surface':payload})
            if u.path=='/api/edit_edge':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper(); sid=d.get('id'); edge=d.get('edge') or []
                if sau not in ALLOWED or not sid or len(edge)!=2: raise ValueError('Aresta inválida.')
                payload=load_surface(sau,sid); mesh,action=remove_edge_from_mesh(payload['mesh'],edge[0],edge[1]); payload['mesh']=mesh; save_payload(sau,payload); return self._json(200,{'ok':True,'action':action,'surface':payload})
            if u.path=='/api/reset_tin':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper(); sid=d.get('id')
                payload=load_surface(sau,sid); payload['mesh']=reset_mesh_from_points(payload); save_payload(sau,payload); return self._json(200,{'ok':True,'surface':payload})
            if u.path=='/api/compare':
                d=json.loads(body or b'{}'); sau=d.get('sau','').upper(); sid=d.get('id'); ref=d.get('reference_id')
                payload=load_surface(sau,sid); refmesh=load_mesh_for_ref(sau,ref); return self._json(200,compare_tins(payload['mesh'],refmesh))
            if u.path=='/api/delete':
                d=json.loads(body or b'{}');sau=d.get('sau','').upper();sid=d.get('id','');p=DATA/sau/f'{sid}.json'
                if sau not in ALLOWED or not p.exists():raise ValueError('Superfície não encontrada.')
                p.unlink();return self._json(200,{'ok':True})
            if u.path=='/api/section':
                d=json.loads(body or b'{}');sau=d.get('sau','').upper();sid=d.get('id');A=(fnum(d['A']['east']),fnum(d['A']['north']));B=(fnum(d['B']['east']),fnum(d['B']['north']));res=section_profile(load_surface(sau,sid)['mesh'],A,B,d.get('samples',500));ref=d.get('reference_id')
                if ref:res['reference']=section_profile(load_mesh_for_ref(sau,ref),A,B,d.get('samples',500))
                cad=load_base_cad(sau)
                if cad:res['base_cad']=section_cad_intersections(cad,A,B)
                return self._json(200,res)
            raise ValueError('Endpoint não encontrado.')
        except Exception as e:return self._json(400,{'error':str(e)})
if __name__=='__main__':
    for s in ALLOWED:(DATA/s).mkdir(parents=True,exist_ok=True)
    print(f'Portal de Superfícies rodando em http://{HOST}:{PORT}')
    ThreadingHTTPServer((HOST,PORT),Handler).serve_forever()

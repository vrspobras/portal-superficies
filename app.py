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
    for v in vs:
        lon,lat=TRANSFORM.transform(v['east'],v['north']);features.append({'type':'Feature','properties':{'id':v['id'],'description':v['description'],'z':v['z'],'east':v['east'],'north':v['north']},'geometry':{'type':'Point','coordinates':[lon,lat]}})
    return {'type':'FeatureCollection','features':features}
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
            if u.path=='/api/delete':
                d=json.loads(body or b'{}');sau=d.get('sau','').upper();sid=d.get('id','');p=DATA/sau/f'{sid}.json'
                if sau not in ALLOWED or not p.exists():raise ValueError('Superfície não encontrada.')
                p.unlink();return self._json(200,{'ok':True})
            if u.path=='/api/section':
                d=json.loads(body or b'{}');sau=d.get('sau','').upper();sid=d.get('id');A=(fnum(d['A']['east']),fnum(d['A']['north']));B=(fnum(d['B']['east']),fnum(d['B']['north']));res=section_profile(load_surface(sau,sid)['mesh'],A,B,d.get('samples',500));ref=d.get('reference_id')
                if ref:res['reference']=section_profile(load_surface(sau,ref)['mesh'],A,B,d.get('samples',500))
                return self._json(200,res)
            raise ValueError('Endpoint não encontrado.')
        except Exception as e:return self._json(400,{'error':str(e)})
if __name__=='__main__':
    for s in ALLOWED:(DATA/s).mkdir(parents=True,exist_ok=True)
    print(f'Portal de Superfícies rodando em http://{HOST}:{PORT}')
    ThreadingHTTPServer((HOST,PORT),Handler).serve_forever()

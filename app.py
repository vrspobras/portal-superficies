from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse
import json, os, uuid, re, csv, io
from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
STATIC = ROOT / 'static'
HOST, PORT = '0.0.0.0', 8000


def finite(v):
    return np.isfinite(float(v))


def parse_txt(text):
    # Normalize BOM/newlines and detect delimiter.
    text = text.lstrip('\ufeff').replace('\r\n','\n').replace('\r','\n')
    lines = [ln.strip() for ln in text.split('\n') if ln.strip()]
    if not lines:
        raise ValueError('Arquivo vazio.')
    delim = ','
    # Header-aware parsing.
    first = lines[0]
    if ';' in first and first.count(';') > first.count(','):
        delim = ';'
    if '\t' in first and first.count('\t') > first.count(delim):
        delim = '\t'
    raw = [next(csv.reader([ln], delimiter=delim)) for ln in lines]

    header = [c.strip().lower() for c in raw[0]]
    # Format A: Name,N,E,Z,...
    if all(x in header for x in ['n','e','z']):
        i_n, i_e, i_z = header.index('n'), header.index('e'), header.index('z')
        i_name = header.index('name') if 'name' in header else None
        pts=[]
        for row in raw[1:]:
            if len(row) <= max(i_n,i_e,i_z):
                continue
            name = row[i_name].strip() if i_name is not None and i_name < len(row) else str(len(pts)+1)
            # Base line must not become a point.
            if name.lower().startswith('base:'):
                continue
            try:
                n,e,z = map(float, [row[i_n], row[i_e], row[i_z]])
            except Exception:
                continue
            if np.isfinite(n) and np.isfinite(e) and np.isfinite(z):
                pts.append({'id': name, 'description': name, 'north': n, 'east': e, 'z': z})
        return pts

    # Format B: Numero do ponto, Descricao, Norte, Este, Cota
    hnorm = [re.sub(r'[^a-z0-9]+','',c) for c in header]
    candidates = {
        'id': ['numerodoponto','numero','ponto','name','nome'],
        'desc': ['descricao','description','desc'],
        'north': ['norte','north','n'],
        'east': ['este','east','e'],
        'z': ['cota','elevacao','elevation','z']
    }
    idx={}
    for key,names in candidates.items():
        for name in names:
            if name in hnorm:
                idx[key]=hnorm.index(name); break
    if all(k in idx for k in ['id','north','east','z']):
        pts=[]
        for row in raw[1:]:
            if len(row) <= max(idx.values()): continue
            try:
                n,e,z = float(row[idx['north']]), float(row[idx['east']]), float(row[idx['z']])
            except Exception:
                continue
            ident = row[idx['id']].strip()
            desc = row[idx['desc']].strip() if 'desc' in idx and idx['desc'] < len(row) else ident
            if finite(n) and finite(e) and finite(z):
                pts.append({'id': ident, 'description': desc, 'north': n, 'east': e, 'z': z})
        return pts

    # Fallbacks based on numeric column patterns.
    pts=[]
    for row in raw:
        if len(row) >= 5:
            # Known original format: index, number, N, E, Z, ..., ..., description
            try:
                int(float(row[0])); ident=row[1].strip(); n=float(row[2]); e=float(row[3]); z=float(row[4])
                if np.isfinite(n) and np.isfinite(e) and np.isfinite(z):
                    desc = row[-1].strip() if row else ident
                    pts.append({'id': ident, 'description': desc, 'north': n, 'east': e, 'z': z})
            except Exception: pass
    if len(pts) >= 3:
        return pts
    raise ValueError('Não foi possível identificar as colunas de ponto, Norte, Este e Cota.')


def build_tin(points):
    if len(points) < 3:
        raise ValueError('São necessários pelo menos 3 pontos válidos.')
    # Remove exact XY duplicates, retain first point.
    unique=[]; seen=set()
    for p in points:
        key=(round(p['east'], 8), round(p['north'], 8))
        if key in seen: continue
        seen.add(key); unique.append(p)
    if len(unique) < 3:
        raise ValueError('São necessários pelo menos 3 pontos com coordenadas planimétricas distintas.')
    xy=np.array([[p['east'],p['north']] for p in unique], dtype=float)
    # localize coordinates for numerical stability
    origin=xy.mean(axis=0)
    local=xy-origin
    try:
        tri=Delaunay(local, qhull_options='QJ')
    except Exception as e:
        raise ValueError(f'Não foi possível formar o TIN: {e}')
    faces=tri.simplices.astype(int).tolist()
    # normalized coordinates kept in local form for browser convenience
    vertices=[]
    zvals=[]
    for p in unique:
        vertices.append({'id':p['id'],'description':p['description'],'east':p['east'],'north':p['north'],'z':p['z'], 'x':p['east']-origin[0], 'y':p['north']-origin[1]})
        zvals.append(p['z'])
    return {'vertices':vertices,'faces':faces,'origin':{'east':origin[0],'north':origin[1]},'zmin':float(min(zvals)),'zmax':float(max(zvals)),'points':len(vertices),'triangles':len(faces)}


def project_json(mesh):
    # Strip redundant fields? keep as is
    return mesh


def list_surfaces(sau):
    folder=DATA/sau
    folder.mkdir(parents=True,exist_ok=True)
    out=[]
    for f in sorted(folder.glob('*.json'), key=lambda p:p.stat().st_mtime, reverse=True):
        try:
            d=json.loads(f.read_text(encoding='utf-8'))
            out.append({'id':d['id'],'name':d['name'],'filename':d.get('filename',''),'created':d.get('created',''),'points':d['mesh']['points'],'triangles':d['mesh']['triangles']})
        except Exception: pass
    return out


def save_surface(sau, name, filename, pts, mesh):
    folder=DATA/sau; folder.mkdir(parents=True,exist_ok=True)
    sid=uuid.uuid4().hex[:10]
    payload={'id':sid,'name':name,'filename':filename,'created':__import__('datetime').datetime.now().isoformat(timespec='seconds'),'points_raw':pts,'mesh':mesh}
    (folder/f'{sid}.json').write_text(json.dumps(payload, ensure_ascii=False), encoding='utf-8')
    return payload


def load_surface(sau,sid):
    p=DATA/sau/f'{sid}.json'
    if not p.exists(): raise FileNotFoundError(sid)
    return json.loads(p.read_text(encoding='utf-8'))


class Handler(BaseHTTPRequestHandler):
    server_version='PortalSuperficies/1.0'
    def _send(self, code, data, ctype='application/json; charset=utf-8'):
        body = data if isinstance(data,bytes) else data.encode('utf-8')
        self.send_response(code); self.send_header('Content-Type',ctype); self.send_header('Content-Length',str(len(body))); self.end_headers(); self.wfile.write(body)
    def do_GET(self):
        u=urlparse(self.path)
        if u.path=='/':
            self._send(200,(STATIC/'index.html').read_bytes(),'text/html; charset=utf-8'); return
        if u.path.startswith('/api/surfaces/'):
            sau=u.path.split('/')[3].upper()
            if sau not in ('SAU3','SAU4'): self._send(400,json.dumps({'error':'SAU inválida'})); return
            self._send(200,json.dumps(list_surfaces(sau),ensure_ascii=False)); return
        if u.path.startswith('/api/surface/'):
            parts=u.path.split('/')
            if len(parts)>=4:
                sau=parts[3].upper(); sid=parts[4] if len(parts)>4 else ''
                try: self._send(200,json.dumps(load_surface(sau,sid),ensure_ascii=False))
                except Exception: self._send(404,json.dumps({'error':'Superfície não encontrada'}))
                return
        if u.path=='/api/test':
            self._send(200,json.dumps({'ok':True}));return
        self._send(404,b'Not found','text/plain')

    def do_POST(self):
        u=urlparse(self.path)
        if u.path=='/api/upload':
            ctype=self.headers.get('Content-Type','')
            length=int(self.headers.get('Content-Length','0'))
            body=self.rfile.read(length)
            # Simple multipart parser without third party libs.
            if 'multipart/form-data' not in ctype:
                self._send(400,json.dumps({'error':'Envio deve ser multipart/form-data'})); return
            m=re.search(r'boundary=(?:"([^"]+)"|([^;]+))',ctype)
            if not m: self._send(400,json.dumps({'error':'Boundary ausente'})); return
            boundary=(m.group(1) or m.group(2)).encode()
            chunks=body.split(b'--'+boundary)
            filename='arquivo.txt'; filebytes=None; name='Nova superfície'; sau='SAU3'
            for ch in chunks:
                if b'filename=' not in ch: continue
                head, sep, content=ch.partition(b'\r\n\r\n')
                if not sep: head,sep,content=ch.partition(b'\n\n')
                if not sep: continue
                fm=re.search(br'filename="([^"]*)"',head)
                if fm: filename=fm.group(1).decode('utf-8','ignore')
                content=content.rsplit(b'\r\n',1)[0] if content.endswith(b'\r\n') else content.rsplit(b'\n',1)[0] if content.endswith(b'\n') else content
                filebytes=content; break
            # fields are parsed from all parts
            for ch in chunks:
                head,sep,content=ch.partition(b'\r\n\r\n')
                if not sep: continue
                if b'name="name"' in head:
                    name=content.rstrip(b'\r\n-').decode('utf-8','ignore')
                if b'name="sau"' in head:
                    sau=content.rstrip(b'\r\n-').decode('utf-8','ignore').upper()
            if filebytes is None: self._send(400,json.dumps({'error':'Nenhum arquivo recebido'}));return
            try: text=filebytes.decode('utf-8-sig',errors='replace')
            except Exception: text=filebytes.decode('latin1',errors='replace')
            try:
                pts=parse_txt(text); mesh=build_tin(pts); payload=save_surface(sau,name or 'Nova superfície',filename,pts,mesh)
                self._send(200,json.dumps({'ok':True,'surface':payload},ensure_ascii=False))
            except Exception as e:
                self._send(400,json.dumps({'error':str(e)},ensure_ascii=False))
            return
        if u.path=='/api/delete':
            length=int(self.headers.get('Content-Length','0')); data=json.loads(self.rfile.read(length) or '{}')
            sau=data.get('sau','').upper(); sid=data.get('id','')
            if sau not in ('SAU3','SAU4'): self._send(400,json.dumps({'error':'SAU inválida'}));return
            p=DATA/sau/f'{sid}.json'
            if not p.exists(): self._send(404,json.dumps({'error':'Não encontrada'}));return
            p.unlink(); self._send(200,json.dumps({'ok':True})); return
        self._send(404,json.dumps({'error':'Endpoint não encontrado'}))

if __name__=='__main__':
    print(f'Portal de Superfícies rodando em http://{HOST}:{PORT}')
    ThreadingHTTPServer((HOST,PORT),Handler).serve_forever()

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse
from pathlib import Path
import json, os, uuid, re, csv, io, math
from xml.etree import ElementTree as ET
import numpy as np
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
STATIC = ROOT / 'static'
HOST, PORT = '0.0.0.0', 8000
ALLOWED = {'SAU3', 'SAU4'}


def clean_tag(tag):
    return tag.rsplit('}', 1)[-1].lower()


def fnum(v):
    return float(str(v).strip().replace(',', '.'))


def parse_delimited(text):
    text = text.lstrip('\ufeff').replace('\r\n', '\n').replace('\r', '\n')
    lines = [x.strip() for x in text.split('\n') if x.strip()]
    if not lines:
        raise ValueError('Arquivo vazio.')
    first = lines[0]
    delimiter = ';' if first.count(';') > first.count(',') and first.count(';') >= first.count('\t') else '\t' if first.count('\t') > first.count(',') else ','
    rows = [next(csv.reader([ln], delimiter=delimiter)) for ln in lines]
    header = [re.sub(r'[^a-z0-9]+', '', c.strip().lower()) for c in rows[0]]
    # Name,N,E,Z,...
    if all(k in header for k in ('n', 'e', 'z')):
        i_n, i_e, i_z = header.index('n'), header.index('e'), header.index('z')
        i_name = header.index('name') if 'name' in header else None
        pts = []
        for row in rows[1:]:
            if len(row) <= max(i_n, i_e, i_z):
                continue
            name = row[i_name].strip() if i_name is not None and i_name < len(row) else str(len(pts) + 1)
            if name.lower().startswith('base:'):
                continue
            try:
                n, e, z = fnum(row[i_n]), fnum(row[i_e]), fnum(row[i_z])
            except Exception:
                continue
            if all(np.isfinite(v) for v in (n, e, z)):
                pts.append({'id': name, 'description': name, 'north': n, 'east': e, 'z': z})
        if len(pts) >= 3:
            return pts

    # Numero do ponto, Descricao, Norte, Este, Cota
    cand = {
        'id': ['numerodoponto', 'numero', 'ponto', 'name', 'nome'],
        'desc': ['descricao', 'description', 'desc'],
        'north': ['norte', 'north', 'n'],
        'east': ['este', 'east', 'e'],
        'z': ['cota', 'elevacao', 'elevation', 'z'],
    }
    idx = {}
    for key, names in cand.items():
        for name in names:
            if name in header:
                idx[key] = header.index(name)
                break
    if all(k in idx for k in ('id', 'north', 'east', 'z')):
        pts = []
        for row in rows[1:]:
            if len(row) <= max(idx.values()):
                continue
            try:
                n, e, z = fnum(row[idx['north']]), fnum(row[idx['east']]), fnum(row[idx['z']])
            except Exception:
                continue
            if not all(np.isfinite(v) for v in (n, e, z)):
                continue
            ident = row[idx['id']].strip() or str(len(pts) + 1)
            desc = row[idx['desc']].strip() if 'desc' in idx and idx['desc'] < len(row) else ident
            pts.append({'id': ident, 'description': desc, 'north': n, 'east': e, 'z': z})
        if len(pts) >= 3:
            return pts

    # Original pattern: index, number, N, E, Z, ..., description
    pts = []
    for row in rows:
        if len(row) < 5:
            continue
        try:
            float(row[0]); ident = row[1].strip(); n, e, z = fnum(row[2]), fnum(row[3]), fnum(row[4])
        except Exception:
            continue
        if all(np.isfinite(v) for v in (n, e, z)):
            pts.append({'id': ident or str(len(pts) + 1), 'description': row[-1].strip() if row else ident, 'north': n, 'east': e, 'z': z})
    if len(pts) >= 3:
        return pts
    raise ValueError('Não foi possível identificar pelo menos 3 pontos com Norte, Este e Cota.')


def parse_landxml(text):
    try:
        root = ET.fromstring(text)
    except ET.ParseError as e:
        raise ValueError(f'XML inválido: {e}')
    points = []
    id_to_idx = {}
    faces_raw = []
    # Common LandXML CgPoints
    for elem in root.iter():
        tag = clean_tag(elem.tag)
        if tag in ('cgpoint', 'p'):
            txt = (elem.text or '').strip()
            vals = re.split(r'[\s,;]+', txt)
            if len(vals) < 3:
                continue
            try:
                a, b, z = fnum(vals[0]), fnum(vals[1]), fnum(vals[2])
            except Exception:
                continue
            # LandXML P/CgPoint coordinates are E N Z.
            ident = str(elem.attrib.get('name') or elem.attrib.get('id') or len(points) + 1)
            desc = str(elem.attrib.get('desc') or elem.attrib.get('description') or ident)
            if clean_tag(elem.tag) == 'cgpoint' and ident.lower().startswith('base'):
                continue
            if ident in id_to_idx:
                continue
            id_to_idx[ident] = len(points)
            points.append({'id': ident, 'description': desc, 'north': b, 'east': a, 'z': z})
        elif tag == 'face':
            vals = re.findall(r'-?\d+', (elem.text or ''))
            if len(vals) >= 3:
                faces_raw.append([vals[0], vals[1], vals[2]])

    # If the XML has faces, translate point ids to indexes.
    faces = []
    for f in faces_raw:
        try:
            idxs = [id_to_idx[x] for x in f]
        except KeyError:
            # Some simple XMLs use zero-based integer vertex indexes.
            try:
                idxs = [int(x) for x in f]
                if min(idxs) < 0 or max(idxs) >= len(points):
                    continue
            except Exception:
                continue
        if len(set(idxs)) == 3:
            faces.append(idxs)
    if len(points) < 3:
        # Support generic <point> with attributes.
        for elem in root.iter():
            if clean_tag(elem.tag) != 'point':
                continue
            at = elem.attrib
            try:
                e = fnum(at.get('e') or at.get('x') or at.get('este'))
                n = fnum(at.get('n') or at.get('y') or at.get('norte'))
                z = fnum(at.get('z') or at.get('cota') or at.get('elev'))
            except Exception:
                continue
            ident = at.get('id') or at.get('name') or str(len(points) + 1)
            if ident in id_to_idx:
                continue
            points.append({'id': ident, 'description': at.get('desc') or ident, 'north': n, 'east': e, 'z': z})
            id_to_idx[ident] = len(points) - 1
    if len(points) < 3:
        raise ValueError('O XML não contém pelo menos 3 pontos reconhecíveis.')
    return points, faces


def build_tin(points, faces_override=None):
    unique, seen = [], set()
    for p in points:
        key = (round(float(p['east']), 8), round(float(p['north']), 8))
        if key in seen:
            continue
        seen.add(key)
        unique.append(p)
    if len(unique) < 3:
        raise ValueError('São necessários pelo menos 3 pontos planimetricamente distintos.')
    xy = np.array([[p['east'], p['north']] for p in unique], float)
    origin = xy.mean(axis=0)
    local = xy - origin
    if faces_override:
        faces = [list(map(int, f)) for f in faces_override if len(f) == 3 and max(f) < len(unique)]
        if not faces:
            raise ValueError('O XML possui faces, mas nenhuma face utilizável; será necessário recalcular o TIN.')
    else:
        try:
            faces = Delaunay(local, qhull_options='QJ').simplices.astype(int).tolist()
        except Exception as e:
            raise ValueError(f'Falha na triangulação TIN: {e}')
    verts = []
    for p in unique:
        verts.append({
            'id': p['id'], 'description': p.get('description', p['id']),
            'east': float(p['east']), 'north': float(p['north']), 'z': float(p['z']),
            'x': float(p['east'] - origin[0]), 'y': float(p['north'] - origin[1])
        })
    return {'vertices': verts, 'faces': faces, 'origin': {'east': float(origin[0]), 'north': float(origin[1])},
            'zmin': float(min(v['z'] for v in verts)), 'zmax': float(max(v['z'] for v in verts)),
            'points': len(verts), 'triangles': len(faces)}


def parse_file(filename, data):
    lower = filename.lower()
    text = data.decode('utf-8-sig', errors='replace')
    if lower.endswith('.xml') or lower.endswith('.landxml'):
        pts, faces = parse_landxml(text)
        return pts, faces, 'LandXML'
    pts = parse_delimited(text)
    return pts, None, 'TXT/CSV'


def list_surfaces(sau):
    folder = DATA / sau
    folder.mkdir(parents=True, exist_ok=True)
    out = []
    for f in sorted(folder.glob('*.json'), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            d = json.loads(f.read_text(encoding='utf-8'))
            m = d['mesh']
            out.append({'id': d['id'], 'name': d['name'], 'filename': d.get('filename', ''), 'created': d.get('created', ''),
                        'source': d.get('source', ''), 'points': m['points'], 'triangles': m['triangles']})
        except Exception:
            continue
    return out


def save_surface(sau, name, filename, source, pts, mesh):
    sid = uuid.uuid4().hex[:10]
    payload = {'id': sid, 'name': name or 'Nova superfície', 'filename': filename, 'source': source,
               'created': __import__('datetime').datetime.now().isoformat(timespec='seconds'),
               'points_raw': pts, 'mesh': mesh}
    (DATA / sau / f'{sid}.json').write_text(json.dumps(payload, ensure_ascii=False), encoding='utf-8')
    return payload


def load_surface(sau, sid):
    if sau not in ALLOWED:
        raise ValueError('SAU inválida')
    p = DATA / sau / f'{sid}.json'
    if not p.exists():
        raise FileNotFoundError(sid)
    return json.loads(p.read_text(encoding='utf-8'))


def save_payload(sau, payload):
    (DATA / sau / f"{payload['id']}.json").write_text(json.dumps(payload, ensure_ascii=False), encoding='utf-8')


def point_in_triangle(pt, a, b, c, eps=1e-9):
    px, py = pt
    v0x, v0y = c[0]-a[0], c[1]-a[1]
    v1x, v1y = b[0]-a[0], b[1]-a[1]
    v2x, v2y = px-a[0], py-a[1]
    den = v0x*v1y - v1x*v0y
    if abs(den) < eps:
        return None
    u = (v2x*v1y - v1x*v2y) / den
    v = (v0x*v2y - v2x*v0y) / den
    if u >= -eps and v >= -eps and u + v <= 1 + eps:
        w = 1 - u - v
        return w*a[2] + v*b[2] + u*c[2]
    return None


def z_at(mesh, e, n):
    vs = mesh['vertices']; faces = mesh['faces']
    for f in faces:
        a0, b0, c0 = (vs[f[i]] for i in range(3))
        val = point_in_triangle((e, n), (a0['east'], a0['north'], a0['z']),
                                (b0['east'], b0['north'], b0['z']),
                                (c0['east'], c0['north'], c0['z']))
        if val is not None:
            return float(val)
    return None


def section_profile(mesh, A, B, samples=500):
    e1, n1 = A; e2, n2 = B
    length = math.hypot(e2-e1, n2-n1)
    if length <= 1e-9:
        raise ValueError('Os pontos A e B da seção não podem ser iguais.')
    samples = max(20, min(int(samples), 1000))
    prof = []
    for i in range(samples + 1):
        t = i / samples
        e = e1 + (e2-e1)*t
        n = n1 + (n2-n1)*t
        z = z_at(mesh, e, n)
        if z is not None:
            prof.append({'distance': length*t, 'east': e, 'north': n, 'z': z})
    return {'length': length, 'profile': prof}


def edge_key(a, b):
    return tuple(sorted((int(a), int(b))))


def edge_edit(mesh, a, b):
    key = edge_key(a, b)
    adjacent = []
    for i, f in enumerate(mesh['faces']):
        edges = (edge_key(f[0], f[1]), edge_key(f[1], f[2]), edge_key(f[2], f[0]))
        if key in edges:
            adjacent.append(i)
    if not adjacent:
        raise ValueError('Aresta não encontrada.')
    faces = [f[:] for f in mesh['faces']]
    if len(adjacent) == 1:
        # Boundary: delete the triangle on that boundary edge.
        faces.pop(adjacent[0])
        msg = 'Aresta de borda removida; o triângulo adjacente foi retirado.'
    else:
        # Interior: remove the two triangles and reconnect the quadrilateral by the opposite diagonal.
        i1, i2 = adjacent[:2]
        f1, f2 = faces[i1], faces[i2]
        opp1 = next(v for v in f1 if v not in key)
        opp2 = next(v for v in f2 if v not in key)
        new_edge = edge_key(opp1, opp2)
        # Do not create a degenerate triangle.
        vs = mesh['vertices']
        def area2(f):
            p, q, r = (vs[i] for i in f)
            return (q['east']-p['east'])*(r['north']-p['north']) - (q['north']-p['north'])*(r['east']-p['east'])
        nf1 = [opp1, opp2, key[0]]; nf2 = [opp1, opp2, key[1]]
        if abs(area2(nf1)) < 1e-8 or abs(area2(nf2)) < 1e-8:
            raise ValueError('A diagonal alternativa degeneraria um triângulo.')
        for idx in sorted((i1, i2), reverse=True):
            faces.pop(idx)
        faces.extend((nf1, nf2))
        msg = f'Aresta {a}–{b} removida e região re-triangulada pela diagonal {new_edge[0]}–{new_edge[1]}.'
    mesh = dict(mesh)
    mesh['faces'] = faces
    mesh['triangles'] = len(faces)
    return mesh, msg


class Handler(BaseHTTPRequestHandler):
    server_version = 'PortalSuperficies/2.0'
    def _send(self, code, data, ctype='application/json; charset=utf-8'):
        body = data if isinstance(data, bytes) else data.encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code, obj):
        self._send(code, json.dumps(obj, ensure_ascii=False))

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == '/':
            self._send(200, (STATIC / 'index.html').read_bytes(), 'text/html; charset=utf-8'); return
        if u.path.startswith('/api/surfaces/'):
            sau = u.path.split('/')[3].upper()
            if sau not in ALLOWED: self._json(400, {'error': 'SAU inválida'}); return
            self._json(200, list_surfaces(sau)); return
        if u.path.startswith('/api/surface/'):
            parts = u.path.split('/')
            sau = parts[3].upper() if len(parts) > 3 else ''
            sid = parts[4] if len(parts) > 4 else ''
            try: self._json(200, load_surface(sau, sid))
            except Exception as e: self._json(404, {'error': str(e)}); return
        if u.path == '/api/test': self._json(200, {'ok': True}); return
        self._send(404, b'Not found', 'text/plain; charset=utf-8')

    def _multipart(self, body, ctype):
        m = re.search(r'boundary=(?:"([^"]+)"|([^;]+))', ctype)
        if not m: raise ValueError('Boundary ausente.')
        boundary = (m.group(1) or m.group(2)).encode()
        parts = body.split(b'--' + boundary)
        fields, file_name, file_bytes = {}, None, None
        for part in parts:
            if b'Content-Disposition' not in part: continue
            head, sep, content = part.partition(b'\r\n\r\n')
            if not sep: head, sep, content = part.partition(b'\n\n')
            if not sep: continue
            fm = re.search(br'filename="([^"]*)"', head)
            nm = re.search(br'name="([^"]+)"', head)
            content = content.rstrip(b'\r\n')
            if fm:
                file_name = fm.group(1).decode('utf-8', 'ignore')
                file_bytes = content
            elif nm:
                fields[nm.group(1).decode('utf-8', 'ignore')] = content.decode('utf-8', 'ignore')
        return fields, file_name, file_bytes

    def do_POST(self):
        u = urlparse(self.path)
        length = int(self.headers.get('Content-Length', '0'))
        body = self.rfile.read(length)
        try:
            if u.path == '/api/upload':
                ctype = self.headers.get('Content-Type', '')
                fields, filename, file_bytes = self._multipart(body, ctype)
                if not file_bytes or not filename: raise ValueError('Nenhum arquivo recebido.')
                sau = fields.get('sau', 'SAU3').upper()
                if sau not in ALLOWED: raise ValueError('SAU inválida.')
                name = fields.get('name', '').strip() or Path(filename).stem
                pts, faces, source = parse_file(filename, file_bytes)
                mesh = build_tin(pts, faces)
                payload = save_surface(sau, name, filename, source, pts, mesh)
                self._json(200, {'ok': True, 'surface': payload}); return

            if u.path == '/api/delete':
                data = json.loads(body or b'{}'); sau = data.get('sau', '').upper(); sid = data.get('id', '')
                p = DATA / sau / f'{sid}.json'
                if sau not in ALLOWED or not p.exists(): raise ValueError('Superfície não encontrada.')
                p.unlink(); self._json(200, {'ok': True}); return

            if u.path == '/api/surface/update':
                data = json.loads(body or b'{}'); sau = data.get('sau', '').upper(); payload = data.get('surface')
                if sau not in ALLOWED or not payload: raise ValueError('Dados inválidos.')
                sid = payload.get('id', '')
                old = load_surface(sau, sid)
                payload['name'] = old['name']; payload['filename'] = old.get('filename', ''); payload['source'] = old.get('source', '')
                payload['points_raw'] = old.get('points_raw', []); payload['created'] = old.get('created', payload.get('created', ''))
                save_payload(sau, payload)
                self._json(200, {'ok': True}); return

            if u.path == '/api/edge/delete':
                data = json.loads(body or b'{}'); sau = data.get('sau', '').upper(); sid = data.get('id', ''); a, b = int(data['a']), int(data['b'])
                surface = load_surface(sau, sid)
                mesh, msg = edge_edit(surface['mesh'], a, b)
                surface['mesh'] = mesh
                save_payload(sau, surface)
                self._json(200, {'ok': True, 'message': msg, 'surface': surface}); return

            if u.path == '/api/section':
                data = json.loads(body or b'{}'); sau = data.get('sau', '').upper(); sid = data.get('id', '')
                A = (fnum(data['A']['east']), fnum(data['A']['north'])); B = (fnum(data['B']['east']), fnum(data['B']['north']))
                surface = load_surface(sau, sid)
                result = section_profile(surface['mesh'], A, B, data.get('samples', 500))
                ref_id = data.get('reference_id')
                if ref_id:
                    ref = load_surface(sau, ref_id)
                    result['reference'] = section_profile(ref['mesh'], A, B, data.get('samples', 500))
                    # Difference profile where both are present.
                    if result['reference']['profile']:
                        ref_by_i = result['reference']['profile']
                        diffs=[]
                        for p in result['profile']:
                            nearest = min(ref_by_i, key=lambda q: abs(q['distance']-p['distance']))
                            diffs.append({'distance': p['distance'], 'delta': p['z']-nearest['z'], 'current_z': p['z'], 'reference_z': nearest['z']})
                        result['differences'] = diffs
                result['A'] = {'east': A[0], 'north': A[1]}; result['B'] = {'east': B[0], 'north': B[1]}
                self._json(200, result); return

            raise ValueError('Endpoint não encontrado.')
        except Exception as e:
            self._json(400, {'error': str(e)})


if __name__ == '__main__':
    for s in ALLOWED: (DATA / s).mkdir(parents=True, exist_ok=True)
    print(f'Portal de Superfícies rodando em http://{HOST}:{PORT}')
    ThreadingHTTPServer((HOST, PORT), Handler).serve_forever()

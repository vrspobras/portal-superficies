# Portal de Superfícies — Codespaces V7

Portal local no GitHub Codespaces para SAU3 e SAU4.

## Novidades desta versão
- Base CAD completa da SAU4 a partir do DXF fornecido pelo usuário, usando a geometria do Model expandida em blocos.
- Camada Base CAD sobre a planta/satélite.
- Seções continuam cortando o TIN e agora recebem interseções da Base CAD quando as entidades têm Z útil.
- Importação TXT/CSV e LandXML.
- TIN, edição de arestas, breakline, comparação e projeto BASE.

## Executar
```bash
pip install -r requirements.txt
python app.py
```
Abra a porta 8000.

## Referencial
SIRGAS 2000 / UTM 22S (EPSG:31982).

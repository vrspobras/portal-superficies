# Portal de Superfícies — Codespaces V3

Portal local/remoto para SAU3 e SAU4 com TIN, TXT/CSV, LandXML e visualização em planta com imagem de satélite.

## Executar
```bash
pip install -r requirements.txt
python app.py
```
Abra a porta 8000 no Codespaces.

## Referencial
Os dados de engenharia são tratados como SIRGAS 2000 / UTM 22S (EPSG:31982). A planta web converte a visualização para latitude/longitude apenas para posicionamento cartográfico.

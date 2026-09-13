# Portal de Superfícies — Codespaces V8

Base do portal para SAU3 e SAU4.

Correções desta versão:
- Ao gerar/abrir uma superfície, a visualização 3D aparece imediatamente.
- Retorno 3D ↔ Planta sem travamento.
- Menu Projetos fica sobre o mapa e pode ser aberto na planta.
- Limite de zoom do satélite para evitar "Map data not yet available".
- Navegação por menus continua disponível em todas as telas.
- Base CAD da SAU4 mantida.

Execução no Codespaces:

```bash
pip install -r requirements.txt
python app.py
```

Abra a porta 8000.

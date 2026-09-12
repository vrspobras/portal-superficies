PORTAL DE SUPERFÍCIES - PYTHON

1) Extraia esta pasta.
2) Dê duplo clique em start.bat.
3) Abra o navegador em http://127.0.0.1:8000

Requisitos:
- Python 3
- numpy
- scipy

Se numpy/scipy não estiverem instalados:
    python -m pip install --user numpy scipy

O portal usa Python no servidor para:
- ler TXT/CSV;
- identificar Norte, Este e Cota;
- remover duplicidades XY;
- calcular a triangulação Delaunay/TIN;
- salvar as superfícies por SAU3/SAU4.

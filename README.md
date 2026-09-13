# Portal de Superfícies — Codespaces V4

Portal de engenharia para SAU3 e SAU4 com processamento de TIN no Python do GitHub Codespaces.

## Executar
```bash
pip install -r requirements.txt
python app.py
```
Abra a porta 8000 no Codespaces.

## Importação
- TXT/CSV de pontos
- XML/LandXML
- TIN automático por Delaunay
- LandXML com faces aproveita a triangulação quando disponível

## Edição TIN
Ative **Cursor CAD** na visualização 3D. O cursor funciona com mira/crosshair, destaca a aresta sob o mouse e permite selecioná-la. **Remover aresta** atualiza a malha no servidor. Em aresta interna, a malha troca para a diagonal oposta; em aresta de contorno, o triângulo adjacente é removido, criando um vazio.

**Restaurar TIN** reconstrói a malha a partir dos pontos originais.

## Projeto final — Base
Cada SAU possui um cadastro separado de **Projeto final · BASE**. A base aceita TXT/CSV ou XML/LandXML e fica salva como referência do projeto.

A base aparece como referência para:
- comparação TIN × TIN de corte/aterro;
- seções A–B comparativas.

## Comparação
A comparação com a BASE usa interseções reais entre os triângulos dos dois TINs e integra a diferença de cota sobre a área comum.

## Referencial
Os dados de engenharia são tratados como **SIRGAS 2000 / UTM 22S — EPSG:31982**. A conversão para latitude/longitude é usada somente para posicionamento cartográfico na planta.

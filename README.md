# Portal de Superfícies — Codespaces V6

Portal de engenharia para SAU3 e SAU4 com processamento TIN no Python do GitHub Codespaces.

## Executar
```bash
python app.py
```
Abra a porta 8000 no Codespaces.

## Interface
- Tela inicial para escolher SAU3 ou SAU4.
- Área de mapa/superfície ocupando a tela.
- Menus para Projetos, TIN, Comparar, Seções e Projeto BASE.
- Alternância 3D / Planta, Orbitar / PAN e Cursor CAD.

## Importação
- TXT/CSV de pontos.
- XML/LandXML.
- LandXML com faces existentes aproveita as faces quando compatíveis.
- TIN automático por Delaunay.

## Edição do TIN
- Cursor CAD com destaque de aresta ao aproximar o mouse.
- Remover aresta sem confirmação intermediária.
- Trocar diagonal.
- Re-Triangular a partir dos pontos originais.
- Breakline: selecione os vértices em sequência e clique em **Aplicar BL**. O backend recupera a linha na malha por operações de arestas, fazendo o TIN respeitar a linha de quebra.
- A breakline fica registrada na evolução da superfície.

## Projeto final — BASE
Cada SAU possui uma BASE separada, aceita TXT/CSV ou XML/LandXML e fica disponível para comparação e seções.

## Comparação
A comparação com a BASE usa interseções reais entre triângulos dos dois TINs e integra a diferença de cota sobre a área comum.

## Seções
- Definição A/B por coordenadas.
- Definição A/B na planta.
- Perfil da superfície.
- Perfil comparativo com outra evolução ou com a BASE.

## Referencial
Dados de engenharia: **SIRGAS 2000 / UTM 22S — EPSG:31982**. A transformação para latitude/longitude é usada somente para posicionamento cartográfico da planta.

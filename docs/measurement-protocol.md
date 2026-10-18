# Ângulos de repouso e de liberação/deslizamento

## O que o projeto mede

Os nomes das respostas são legados. `discharge_deg` identifica a inclinação média da superfície livre após a remoção da base no ensaio `Ang Deslizamento`; não deve ser interpretada automaticamente como ângulo universal de descarga de silo nem como ângulo de atrito interno. `repose_deg` pertence ao ensaio distinto `Ang Repouso`. Os dois valores não são intercambiáveis: decks, condições de contorno, sequência temporal e método de extração são diferentes.

A análise desta versão mantém os dois ensaios, a geometria e os decks. No deslizamento, a seleção usa o método de perfil já existente e aplica o mesmo limite inferior do repouso: por face, inicia o ajuste a 2 diâmetros médios acima da menor altura de centro de partícula naquela seção. Isso elimina as camadas basais da seleção, mas não troca o seletor do ensaio de deslizamento pelo seletor radial de repouso. Em cada lado, o perfil mantém o ajuste único quando não há reversão sustentada; quando a reversão sustentada é detectada, são usados até dois segmentos, com ângulos/linhas reportados separadamente e suporte mínimo exigido para ambos. Não remove reversões sustentadas como outliers. No repouso, preservam-se os quatro cortes diametrais/oito faces, mas o seletor e o ajuste são revisados conforme a seção específica abaixo.

## Cortes de deslizamento e coordenadas físicas

O estimador mede três secções transversais cujos centros físicos são `-0.100 m`, `0.000 m` e `+0.100 m`. O eixo horizontal do perfil é Z, o eixo vertical é Y e o corte é selecionado ao longo de X. As secções são retornadas na ordem `[-0.1, 0.0, +0.1] m`; os nomes dos campos do CSV, porém, são semânticos:

- `discharge_xminus01_*` → centro físico `-0.100 m`;
- `discharge_x0_*` → centro físico `0.000 m`;
- `discharge_xplus01_*` → centro físico `+0.100 m`.

A persistência dos resultados associa cada secção ao seu centro físico explícito; não depende da ordem da lista. As colunas `*_center_physical_m`, `*_band_min_m` e `*_band_max_m` permitem conferir a identidade do corte e os limites efetivos usados. `*_band_clipped` indica que a banda nominal foi truncada pelos limites da malha escalada.

A seleção de partículas usa a interseção entre a banda nominal do corte e os limites da malha naquele eixo. A banda efetiva usada pelo estimador é a mesma banda usada para desenhar o fundo de partículas da figura oficial. Assim, pontos fora da malha não podem influenciar o ajuste nem aparecer como se pertencessem ao corte. Um corte truncado continua explicitamente marcado; a figura não oculta essa limitação.

Para cada lado de cada corte, o limite inferior do perfil é `altura mínima da face + 2 × diâmetro mediano` (equivalente a quatro raios medianos), seguindo o valor `minimum_height_particle_diameters=2.0` do repouso. A seleção por altura é aplicada antes de se escolher a amostra da banda; mantém-se o seletor de deslizamento mais próximo à linha central. O CSV grava esse fator e as alturas Y efetivas por lado. Na figura oficial, as linhas tracejadas coloridas indicam os limites; pontos abaixo delas são contexto da nuvem, não amostras do ajuste. Se a medição for rejeitada, a figura mantém os três cortes e as linhas de limite, mostra `NO VALID ANGLE / NO FIT` e não inventa ângulos nem linhas de ajuste.

Quando há reversão sustentada no perfil, o algoritmo pode dividir cada lado em duas partes e ajustar uma reta por parte. Os ângulos e suportes por segmento ficam separados no CSV e na figura; a média de pesos iguais só é devolvida quando ambos os segmentos atingem o suporte mínimo. Se qualquer parte não tiver suporte suficiente, o corte é rejeitado, sem forçar um ângulo.

## Revisão do seletor de repouso

Esta revisão modifica a extração/representação de repouso e seu diagnóstico. Na descarga, acrescenta somente círculos com o raio individual lido do VTK, projetados na mesma banda espessa já usada pelo estimador; não troca a banda por um plano de espessura zero. O raio em pixels é `raio_VTK × escala_igual_XY`, sem fator de ampliação arbitrário. A base translúcida segue fora do ajuste. O estimador de descarga, seus três cortes, os decks e os parâmetros físicos permanecem inalterados. A prévia de repouso foi aprovada em 17/10/2026, com autorização para nova campanha de 12 trials independentes de um rank.

- As quatro seções diametrais continuam fixas em 0°, 45°, 90° e 135°, passando por `(x,z)=(0,0)`. Uma esfera só participa se interceptar o plano. Seu raio na seção é `sqrt(r_particula² - distancia_normal²)`; o raio original 3D, não o raio cortado, define a resolução vertical e o limite inferior.
- Por face, centros abaixo de `min(Y_centro) + 4 * raio_mediano_3D` são inelegíveis. A base translúcida é referência geométrica, nunca suporte do ajuste.
- A superfície é amostrada em alturas regulares com passo padrão `0.25 * raio_mediano_3D`. Cada esfera que cruza o plano e a altura fornece a coordenada horizontal da sua superfície. Escolhe-se a superfície horizontalmente mais próxima ao eixo lateral externo da seção. Não se substitui uma superfície ausente por uma partícula fictícia. Esferas que cruzam o eixo central podem contribuir para ambas as faces; suas calotas também são consideradas.
- Os eixos laterais são os limites externos das superfícies das esferas daquela seção, não uma única partícula da primeira faixa e não âncoras do ajuste. A origem desse eixo é registrada em `axis_source`. `profile_selection=base_edge_nearest` é um identificador legado; os campos `height_sampling_method`, `repose_axis_reference` e `profile_method` identificam a implementação efetiva.
- A reta tem intercepto livre. Um consenso determinístico de pares separa amostras incompatíveis da face predominante; o refit e o R² usam pesos dos intervalos físicos de altura. Amostras repetidas na mesma altura dividem seu peso. O número de partículas distintas é verificado separadamente e não é confundido com o número de alturas amostradas.
- Políticas declaradas desta revisão: tolerância ortogonal padrão de `1.75 * raio_mediano_3D`, pelo menos oito amostras ajustadas, cinco partículas distintas e R² ponderado dos inliers ≥ 0.90. Não são limiares universais de literatura. Lacunas verticais maiores que dois diâmetros medianos ou desvios sustentados de múltiplas partículas invalidam a face, em vez de serem disfarçados como outliers isolados.
- A figura mantém o formato 2×2. A nuvem cinza e os círculos pertencem à mesma seção física usada pelo seletor, e não à projeção de toda a pilha. Verde indica superfície escolhida por altura; X laranja indica amostras excluídas, não pontos ainda incluídos com peso reduzido. A reta fica restrita ao suporte acima do corte.
- Uma face inválida não recebe escalar nem reta ajustada. Se qualquer face for inválida, não existe média global utilizável no objetivo. O diagnóstico carrega as mesmas seções, eixos, limites e amostras. Falhas anteriores à geração desses metadados têm apenas um diagnóstico reduzido, sem inventar seleção ou limites.

Os testes incluem um ponto interno na primeira faixa, inversão de resíduos, PSD, esfera central, calota superior, multiplicidade, lacunas/desvios sustentados, um perfil em coluna sem face inclinada e uma regressão derivada de superfície histórica para verificar estabilidade ao refinar a malha de alturas.

## Validade, agregação e falhas

- Os ângulos por corte e por face ficam disponíveis antes das médias. A agregação existente e os campos de resposta são mantidos.
- Cortes sem suporte geométrico suficiente continuam sendo rejeitados; o código não inventa um ângulo para completar a média.
- A figura oficial de descarga permanece no formato de três cortes, inclusive quando não há ângulo válido.
- A alteração de limites de seleção pode mudar o resultado se existirem partículas fora dos limites da malha na banda nominal; por isso, resultados de campanhas anteriores não devem ser considerados equivalentes sem reprocessar os dumps originais com a versão identificada do método.
- Nenhum deck, parâmetro de material, alvo, limite DOE ou duração física é alterado por esta correção.

## Verificações automatizadas

Os testes cobrem: associação dos campos aos centros físicos, persistência de centros/bandas/clipping, exclusão de partículas fora dos limites escalados, invariância dos ângulos quando apenas pontos exteriores à malha são acrescentados, e consistência entre os pontos selecionados pelo ajuste e os exibidos na figura. A suíte completa e o preflight do controlador devem passar antes de iniciar uma campanha.

## Proveniência da medição

Cada campanha deve registrar a revisão do código, os decks e a configuração utilizados. Resultados, dumps, logs, figuras por trial e demais artefatos da campanha permanecem fora do repositório público, salvo autorização específica. Uma correção do estimador não reclassifica retroativamente resultados antigos; os artefatos precisam ser reprocessados em uma cópia separada, preservando originais.

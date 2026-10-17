# Ângulos de repouso e de liberação/deslizamento

## O que o projeto mede

Os nomes das respostas são legados. `discharge_deg` identifica a inclinação média da superfície livre após a remoção da base no ensaio `Ang Deslizamento`; não deve ser interpretada automaticamente como ângulo universal de descarga de silo nem como ângulo de atrito interno. `repose_deg` pertence ao ensaio distinto `Ang Repouso`. Os dois valores não são intercambiáveis: decks, condições de contorno, sequência temporal e método de extração são diferentes.

A análise desta versão mantém os dois ensaios, a geometria e os decks. No deslizamento, preserva o estimador aprovado de até dois segmentos sustentados (inferior/superior) e a média dos segmentos quando ambos têm suporte. Não remove reversões sustentadas como outliers. No repouso, mantém o estimador e a composição dos quatro cortes diametrais/oito faces.

## Cortes de deslizamento e coordenadas físicas

O estimador mede três secções transversais cujos centros físicos são `-0.100 m`, `0.000 m` e `+0.100 m`. O eixo horizontal do perfil é Z, o eixo vertical é Y e o corte é selecionado ao longo de X. As secções são retornadas na ordem `[-0.1, 0.0, +0.1] m`; os nomes dos campos do CSV, porém, são semânticos:

- `discharge_xminus01_*` → centro físico `-0.100 m`;
- `discharge_x0_*` → centro físico `0.000 m`;
- `discharge_xplus01_*` → centro físico `+0.100 m`.

A persistência dos resultados associa cada secção ao seu centro físico explícito; não depende da ordem da lista. As colunas `*_center_physical_m`, `*_band_min_m` e `*_band_max_m` permitem conferir a identidade do corte e os limites efetivos usados. `*_band_clipped` indica que a banda nominal foi truncada pelos limites da malha escalada.

A seleção de partículas usa a interseção entre a banda nominal do corte e os limites da malha naquele eixo. A banda efetiva usada pelo estimador é a mesma banda usada para desenhar o fundo de partículas da figura oficial. Assim, pontos fora da malha não podem influenciar o ajuste nem aparecer como se pertencessem ao corte. Um corte truncado continua explicitamente marcado; a figura não oculta essa limitação.

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

# Revisão do repouso e raios da descarga — 17/10/2026

## Escopo aprovado

- Quatro cortes diametrais, oito faces; seção real das esferas, amostragem regular da superfície, intercepto livre, consenso robusto e diagnóstico sem ângulo em faces sem suporte.
- A imagem de descarga passa a desenhar o raio individual do VTK na mesma escala X/Y do perfil. Trata-se de projeção na banda efetiva original, não de troca do seletor ou do ajuste aprovado.
- Não foram alterados decks, duração, parâmetros, alvos ou limites da campanha.
- Campanhas e dumps anteriores permanecem preservados e fora do Git. Scripts de auditoria locais não são dependências de execução.

## Checagem dos Trials 07 e 12

Fonte: campanha histórica `calibration_runs_20261017_163014`, decks e logs de `Ang Repouso`, parâmetros por trial e snapshot `post/repouso_10s.vtk`.

- O nome `repouso_10s.vtk` é legado: o snapshot final é escrito no passo 130000 com dt=0,0001 s, isto é, aos **13 s totais**. São 1 s de alimentação, 2 s de espera e 10 s de subida do tubo. Não existe espera adicional depois da subida. Não se deve rotular este snapshot como 10 s totais.
- Trial 07: CED11 = 11264,121845217265; sf11 = 0,9142933221924179; rf11 = 0,8483472063817843.
- Trial 12: CED11 = 5805461,533862077; sf11 = 0,35543555060020193; rf11 = 1,0899367080251015.
- A coesão do Trial 12 é 515,394 vezes a do 07. Atritos e sementes também diferem: isso não isola causalidade.
- O STL do tubo tem Y inicial de aproximadamente 0 a 0,5 m; escala 37,5 e subida de 6,25 m colocam sua borda inferior em Y=6,25 m na captura.
- O topo das superfícies das esferas fica em Y=4,12222 m no Trial 12 e 2,27904 m no 07. Não há sobreposição vertical com o tubo na captura final.
- O Trial 12 é um diagnóstico de perfil sem faces livres confiáveis; a ausência de confinamento e a coesão elevada são compatíveis com a coluna, não uma prova de que coesão é sua única causa. O instante final não comprova acomodação completa.

## Verificações

O protocolo exige testes completos em WSL, preflight no projeto promovido e verificação dos processos/logs após lançamento. A revisão acrescenta regressões de raios individuais e propagação pelo controlador. Comparação real de descarga nos Trials 01, 02, 07 e 12 conserva integralmente resultados e motivos de rejeição do estimador anterior; por exemplo, o Trial 07 mantém média de 43,904858368641975°. O AST das funções de medição não relacionadas ao repouso permanece idêntico.

## Referências metodológicas

- Müller, Fimbinger e Brand (2021): contorno externo e regressão — https://pure.unileoben.ac.at/ws/portalfiles/portal/7308811/1_s2.0_S0032591021000176_main.pdf
- Consenso robusto como referência, sem afirmar uso literal de scikit-learn — https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RANSACRegressor.html

O algoritmo usa consenso determinístico próprio. Passo e tolerâncias são políticas declaradas do projeto, não consenso universal da literatura. Métricas numéricas não substituem inspeção visual e validação física.

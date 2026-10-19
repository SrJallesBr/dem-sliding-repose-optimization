# Otimização — repouso e abertura sequenciais (Jenike)

## Escopo e estado

Candidato de integração para o caso fornecido por Henrique em `Ang Repouso e Desgarga/in.repousodesgarga` (grafia original preservada). O deck e as quatro malhas **não são alterados** pelo controlador no projeto-fonte. Cada trial recebe uma cópia limpa: somente o deck e `meshes/`; dumps/logs/restarts históricos não são copiados.

- Um único solver por trial; não existem mais dois ensaios independentes.
- `post/repouso_final.vtk`: após formar a pilha, em **13 s**.
- Duas metades da base se movem em sentidos opostos de Z de **13 a 16 s**.
- `post/descarga_final.vtk`: após mais 7 s, em **23 s**.
- O vão final calculado das malhas e velocidades é **0,150 m**.
- O usuário confirmou: medir as **faces internas voltadas ao vão**; alvos **40°** no repouso e **73,4°** após abertura.

**Prévia aprovada por Henrique em 18/10/2026. Não há campanha iniciada ou autorizada neste salvamento.** O flag `measurement_approved=false` foi preservado nesta cópia para impedir lançamento involuntário; a aprovação visual está registrada em `APPROVAL.json`. Os alvos já estão confirmados. O controlador antigo é mantido como base de algoritmos e compatibilidade de regressões; a entrada deste projeto é `optimizer.opening_controller`.

## Medições

### Repouso — preservado

Mesmo backend aprovado: quatro seções diametrais, oito faces, partículas/esfera-plano, intercepto livre e validação geométrica. O conteúdo de `repose_profile.py` e a AST da função de repouso foram comparados com o projeto anterior e são idênticos. Não se reutiliza um dump de outra simulação: mede-se o estado aos 13 s deste próprio trial.

### Abertura — um corte central, não três

Uma seção central **ZY em X=0**, mantendo a meia-largura de observação da descarga anterior: **0,075 m**. Trata-se da projeção das partículas dessa faixa, com seus raios VTK completos; **não** é uma interseção matemática de espessura zero. O tipo da observação e a espessura estão explícitos no JSON, CSV e PNG. Preservar a faixa anterior evita confundir uma mudança no número de cortes com uma mudança silenciosa na espessura.

O contorno da interface material-ar voltada ao vão é obtido em passos de **0,25 raio mediano**, pela interseção de linhas horizontais com os círculos projetados. Nesta revisão, o corte inferior da descarga foi removido (`minimum_height_particle_diameters=0`): o perfil cobre desde os centros mais baixos até a crista, sem descartar a faixa basal. Centros abaixo da superfície física de suporte permanecem inelegíveis; isso não é o antigo corte de dois diâmetros. O contorno inteiro é candidato ao consenso; pontos incompatíveis da crista ficam registrados e visíveis em laranja, com peso zero no fit, sem excluir indiscriminadamente o topo. A medição de repouso não mudou. Raios individuais corrigem a superfície, não apenas os centros das partículas.

O ajuste usa o motor de consenso geométrico já testado no repouso, **sem alterar seus limiares**: intercepto livre, cinco partículas distintas, extensão física, resíduos ortogonais e continuidade. Pontos excluídos são identificados como excluídos, não como Huber. Uma face inválida não recebe escalar; ausência de resposta não desaparece do denominador. A média é calculada apenas quando **ambas** as faces internas são válidas.

Isso reutiliza componentes testados, mas não é o encaixe literal do antigo fit ancorado e não é uma alegada implementação integral de Müller ou certificação do ensaio. O encaixe direto da antiga amostragem por centros/diâmetro foi testado e rejeitado no baseline real. A alternativa de plano de espessura zero também foi preservada como diagnóstico de desenvolvimento; ela não é a observação ativa.

## Campanha preparada

Configuração herdada, sem estreitar limites materiais: `sf11=0,2..2,0`, `rf11=0,2..2,0`, `CED11=0..10.000.000`; iniciais do deck `0,5 / 0,8 / 1.000`. 12 trials concorrentes, 1 rank cada, reposição imediata de vaga. Plano-base: 45 DOE, 48 refinamento, 15 validações independentes; o controlador existente pode gerar contingência se faltar DOE válido.

A redução automática de limite de coesão está desativada nesta integração: rejeição geométrica não prova causalidade da coesão e o protocolo novo não tem um teto empírico validado. Os limites confirmados não são alterados. Nenhuma campanha antiga é importada/reaproveitada.

## Verificar sem iniciar simulação

Em WSL, na pasta deste projeto:

```bash
python3 -B -m unittest discover -s tests -v
python3 -B -m optimizer.opening_controller --preflight --parallel-trials 12 --run-root verification_preflight
python3 -B -m optimizer.opening_controller --prepare-only --run-root verification_prepare
```

`--prepare-only` cria uma cópia de trial e reescreve os parâmetros/sementes, mas **não** chama LIGGGHTS. Não reutilize uma pasta de trial existente.

`verificar_otimizacao.cmd` executa apenas preflight em uma janela visível. `run_visible.cmd` é a entrada de uma **campanha nova**, não uma retomada. Ele está bloqueado pela aprovação de medição pendente; não editar o flag para contornar a revisão visual.

## Artefatos por trial

- `sequential.log` e log original do solver dentro do caso.
- `parameters.json`: parâmetros, sementes, hashes do input fonte/copiado e tempos.
- Dois dumps no mesmo `post/`.
- `repose_measurement.json`, `discharge_measurement.json`: metadados, validade e hash do dump.
- Figuras de repouso 2×2 e abertura em um único painel, SVG + PNG.
- `results.csv`: apenas pares completos, sem colunas de cortes ±0,1 m.
- `partial_responses.json`: resposta preservada quando só a outra foi rejeitada.
- `status.csv`, `diagnostics.csv`, `status_execucao.csv`: falhas de solver/artefato separadas das rejeições geométricas.

## Validação real e limites

Foi executada uma cópia isolada do deck original, sem alterar os parâmetros, física, tempos ou malhas: final em 23 s, 3.657 partículas restantes, saída zero. Os resultados foram então reprocessados pelo controlador para verificar os dois cálculos, figuras, CSV e score; isso **não foi uma campanha de calibração**.

Baseline `sf11=0,5`, `rf11=0,8`, `CED11=1.000`:
- Repouso: **30,6867949661°**.
- Abertura: esquerda **49,1639103925°**, direita **50,8098792032°**, média **49,9868947979°**; suportes de 14 e 17 partículas distintas. Valores da revisão sem corte inferior, reprocessando o mesmo dump salvo.

Auditoria desta revisão sem corte inferior, mesma faixa ±0,075 m: passos de 0,125/0,25/0,5 raio mantiveram ambas as faces mensuráveis e deram médias de 49,826774° a 49,986895° (amplitude 0,160121°). Escalas globais 0,1 e 10 preservaram o ângulo. Não é intervalo de confiança nem ensaio adicional. A faixa física continua fixa; não foi escolhida pelo alvo. A imagem usa o padrão gráfico do repouso, adaptado para um único painel ZY com as duas faces internas e sem linha de cutoff inferior.

Esses valores não atingem os objetivos e não representam uma otimização concluída. Um único baseline não demonstra cobertura de 90%, convergência, estabilidade temporal ou robustez em toda a faixa de parâmetros. A prévia do observável sem corte inferior foi aprovada; o início de uma campanha continua dependendo de autorização explícita.

O dump de descarga originalmente fornecido tinha 611.241 bytes, todos zero, e foi preservado. A verificação produziu um dump novo fora da fonte. O primeiro teste do executável com `-screen none` sofreu SIGSEGV; a repetição isolada com saída normal concluiu com código zero. O runner desta integração não usa `-screen none`. A causa do arquivo histórico nulo não foi estabelecida.

Backup integral da fonte e evidências desta integração: `C:/DEM/_opening_integration_20261018/`. Os testes de pipeline usam respostas analíticas explicitamente sintéticas apenas para exercitar DOE/refino/validação; não são resultados DEM.

# Otimização DEM — deslizamento + repouso

Pacote limpo para iniciar campanha nova. Não contém a campanha interrompida, resultados, dumps, logs de trials ou configurações legadas. A extração não inicia nada.

## Campanha
- Objetivos: deslizamento 73,4° e repouso 40° (não são resultados garantidos).
- Limites: sf11 0,2–2,0; rf11 0,2–2,0; CED11 0–10.000.000 (limites ampliados conforme solicitado).
- DOE 45, refinamento 48, validação 15; seed 86028157.
- Limite desta cópia: **12 trials × 1 rank = 12 processos LIGGGHTS** (12 CPUs WSL disponíveis); não altera DOE, parâmetros nem protocolo físico.

## Requisitos
- Windows + WSL2 com distro `Ubuntu-24.04`.
- Em WSL: Python 3 e LIGGGHTS executável em `/usr/local/bin/liggghts`.
- Espaço livre para dumps/logs da campanha.

## Executar e validar
1. Extraia em pasta local.
2. Execute `run_visible.cmd`; o painel fica no terminal e cada trial mantém logs próprios. A campanha cria `calibration_runs_<data_hora>` novo com 12 processos LIGGGHTS, um rank por trial.
3. Para validação opcional no WSL:
   ```bash
   python3 -m unittest discover -s tests -v
   python3 -u -m optimizer.calibration_controller --preflight --parallel-trials 12 --run-root /tmp/dem_preflight
      python3 -u -m optimizer.calibration_controller --smoke --smoke-trials 1 --parallel-trials 1 --run-root /tmp/dem_smoke
   ```
   O smoke executa uma única amostra real, não a otimização inteira.
4. Não apague pastas de run. Para parar, Ctrl+C uma vez e preserve logs/dumps.

## Artefatos de cada trial
- Após medições válidas, cada `trial_####` recebe `medicao_descarga.svg` + `.png` com os cortes em X=0, +0,1 m e −0,1 m na escala física 1:1 e os dois lados de cada corte; o CSV guarda a média de cada corte, os seis lados e as bandas efetivas. Centros e meia-largura física são multiplicados pela `nbrratio` do deck trial-local. A geometria desta cópia vai de X=−0,125 m a +0,125 m: as bandas nominais de ±0,05 m nos cortes ±0,1 m são truncadas pela parede para [−0,125,−0,05] e [+0,05,+0,125] m; imagem/CSV sinalizam o truncamento. A nova campanha permanece condicionada à revisão desta medição.
- O arquivo legado `medicao_repouso_oito_direcoes.svg` + `.png` apresenta quatro cortes diametrais (0°,45°,90°,135°), os oito lados individuais e as médias por par/global. A seleção usa a interseção das esferas das partículas com cada plano; os perfis reutilizam as bandas de 2× raio mediano, correção pelo raio individual, âncora basal e ajuste Huber da descarga. Mostra nuvem inteira projetada em cinza, centros selecionados em verde, bases em amarelo e pontos de alto resíduo em laranja incluídos com peso menor; eixo central fixo em (0,0) e referência horizontal Y=0 na base. O CSV guarda as oito faces individuais e as quatro médias por par.
- O deck de deslizamento grava somente `post/deslizamento_final.vtk`, o snapshot de partículas usado na medição. Foram removidos os dumps STL periódicos da base/caixa/parede; os limites da caixa são derivados do STL-fonte do próprio trial e do fator de escala declarado no deck. O deck de repouso não foi alterado por essa limpeza.
- Os SVG são evidência vetorial derivada dos dumps, não captura de tela da interface 3D.

## Correção e limites conhecidos
Para o deslizamento, a seleção usa uma amostra por banda vertical de 2× o raio mediano e corrige a superfície pelo raio individual. O ajuste é ancorado na base e inclui todas as amostras selecionadas da base ao topo; resíduos acima de 1,75× o raio mediano apenas reduzem o peso robusto — não removem pontos. R² é mostrado como qualidade, não como veto. Ainda são exigidos ao menos cinco bandas e tendência global para fora; um perfil que se volta para dentro não recebe ângulo. Isso mantém diagnósticos para geometrias sem face livre mensurável.

As verificações de medição incluem reprocessamento somente de leitura de dumps arquivados; esses dumps e relatórios locais não fazem parte deste repositório. Isso não é nova simulação nem validação de repouso físico. As figuras de medição são geradas pelo controlador durante a campanha; consulte os testes e o código da versão publicada para o formato ativo.

O deck de deslizamento agora usa duração fixa máxima de 15 s, gerando o snapshot final nessa marca. A parada antecipada quando a energia cinética total ficar abaixo de 1 J após 10 s não foi ativada: o binário LIGGGHTS-PUBLIC 3.8.0 não lista `fix halt`, e a alternativa condicional testada não passou no probe. 15 s também não garante por si só repouso; será preciso conferir energia/logs e movimento nos dumps reais.

Para comparar PCs, use mesmo LIGGGHTS build, decks/meshes e seed. O DOE inicial é reproduzível; ordem/quantidade das reposições pode variar se houver rejeições ou tempos de conclusão diferentes.

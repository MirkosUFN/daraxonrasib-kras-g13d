# Controle KRAS selvagem (9BG9)

Tri-complexo **daraxonrasib · CypA · KRAS WT (Gly13)** montado a partir do
cristal 9BG9 (1.58 Å) e submetido a **dinâmica molecular de 2 ns com protocolo
idêntico** ao do mutante G13D (9BG5), para isolar o efeito da mutação sobre a
estabilidade do complexo.

## Conteúdo
- `9BG9_gromacs_system.tar.gz` — sistema pronto para GROMACS (topologia, .itp
  dos ligantes, coordenadas solvatadas/neutralizadas, MDPs, index.ndx).
- `9BG9_md_trajectory.tar.gz` — .tpr + trajetória limpa (PBC tratada) do
  complexo, frame inicial e CSVs de análise.
- `complex_md_frame0.pdb` — primeiro frame do complexo (visualização).
- `fig_md_9BG9.png` — 4 painéis: RMSD backbone, RMSF por resíduo (KRAS/CypA),
  distância mínima fármaco–parceiro, contatos <4 Å ao longo do tempo.
- `fig_comparativo_G13D_vs_WT.png` — comparação direta de estabilidade.
- `md_analysis_summary_9BG9.csv` — métricas do WT.
- `comparativo_G13D_vs_WT.csv` — tabela lado a lado G13D vs WT.

## Resultado
O tri-complexo WT é mais estável (RMSD backbone 1.76 Å vs 4.21 Å; RMSF médio
0.95 Å vs 3.82 Å), e o fármaco engaja **KRAS e CypA simultaneamente em 100 %
dos frames em ambos os sistemas** — a cola molecular RAS(ON)·fármaco·CypA se
forma independentemente do status G13.

*Trajetória curta (2 ns, CPU) — estimativa de modelagem, não substitui
validação experimental nem amostragem estendida.*

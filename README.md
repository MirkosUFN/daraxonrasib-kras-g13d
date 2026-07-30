# Daraxonrasib (RMC-6236) × KRAS G13D — Estudo *in silico*

Docking comparativo (AutoDock Vina) e dinâmica molecular (GROMACS) da ação da
**daraxonrasib** sobre **KRAS G13D**, considerando o mecanismo real do fármaco.

## Mecanismo (essencial para interpretar tudo aqui)

Daraxonrasib é um **inibidor tri-complexo não-covalente RAS(ON) multi-seletivo**
("molecular glue"). Ela **não se liga ao KRAS isolado**: liga-se primeiro à
**ciclofilina A (CypA / PPIA, UniProt P62937)**, e o binário **CypA·fármaco**
engaja as alças switch-I/II do RAS no **estado ativo (GTP/GppNHp)**. Portanto o
docking e a MD só são mecanicamente válidos **com CypA presente**; KRAS isolado
serve como **controle negativo**.

- Ligante: código PDB **A1AHB**, C44H58N8O5S, MW ≈ 811, CAS 2765081-21-6,
  InChIKey `FVICRBSEYSHKFY-JYQNNKODSA-N`.
- Alvo: KRAS (UniProt **P01116**), mutação **G13D**.

## Seleção de PDBs

| PDB  | Papel | Conteúdo | Resolução |
|------|-------|----------|-----------|
| **9BG5** | Âncora | KRAS **G13D** + CypA + GppNHp (GNP) + Mg + fármaco (A1AHB); 2 tri-complexos na ASU | 1.67 Å |
| **9BG9** | Controle WT | KRAS WT (Gly13) + CypA + GppNHp + fármaco | 1.58 Å |
| **4TQA** | Controle negativo | KRAS G13D isolado, GDP, sem CypA | 1.13 Å |
| **8BLR** | Controle conformacional | KRAS4b G13D, switch-I aberto, GDP, sem CypA | 1.40 Å |
| 9BGA / 9BG6 | Contexto | G12C / G12V (opcional) | — |

## Resultados do docking (AutoDock Vina, exhaustiveness 32, seed 42)

Redocking do 9BG5 validado: **RMSD 0.57 Å** vs. cristal.

| Sistema | Mutação | Estado | CypA | Afinidade top (kcal/mol) |
|---------|---------|--------|------|--------------------------|
| **9BG5** | G13D | ativo (GppNHp) | Sim (tri) | **−10.59** |
| **9BG9** | WT | ativo (GppNHp) | Sim (tri) | **−10.03** |
| 4TQA_tri | G13D | inativo (GDP) | CypA superposta | −7.44 |
| 8BLR_tri | G13D | inativo (GDP), switch-I aberto | CypA superposta | −5.69 |
| 4TQA_kras | G13D | inativo (GDP) | Não (isolado) | −5.67 |
| 8BLR_kras | G13D | inativo (GDP), switch-I aberto | Não (isolado) | −5.10 |

**Leitura mecanística:** a afinidade só é alta quando há CypA formando o
tri-complexo no estado ativo (9BG5/9BG9 ≈ −10 kcal/mol). Sem CypA, ou no estado
inativo/switch aberto, a afinidade cai ~4–5 kcal/mol — coerente com o
"molecular glue" não reconhecer KRAS isolado.

## Análise de interações (ProLIF)

Poses top de 9BG5/9BG9 reproduzem o farmacóforo do co-cristal: contatos com o
switch-I do KRAS (Tyr32/Pro34/Thr35/Ile36 na numeração do modelo), switch-II
(Ala59, Gln61, Tyr64 com π-stacking, Met67) **e** 11–12 resíduos de CypA
(Arg55, Ile57, Phe60, Gln63, Trp121…). Os controles **KRAS-isolado registram
ZERO contatos com CypA** — confirmação do mecanismo.

## Dinâmica molecular (GROMACS 2024.1)

Sistema montado a partir do tri-complexo âncora do 9BG5 (KRAS-A G13D + CypA-D +
fármaco + GppNHp + Mg):

- **Campo de força:** amber99sb-ildn + TIP3P (ff14SB não disponível no
  GROMACS 2024.1; amber99sb-ildn é compatível com GAFF2).
- **Ligantes:** fármaco (carga 0, 116 átomos) e GNP/GppNHp (carga −4, 45 átomos)
  parametrizados via **acpype / GAFF2 (AM1-BCC)**; Mg²⁺ como íon do campo de força.
- **Caixa:** dodecaédrica, solvatada (TIP3P), 0.15 M NaCl neutralizado →
  **58 290 átomos**, sistema neutro.
- **Protocolo:** minimização (steepest descent) → NVT 100 ps → NPT 100 ps →
  produção. MDPs em `05_gromacs/mdp/`.

O sistema pronto para rodar está empacotado em
`05_gromacs/9BG5_gromacs_system.tar.gz`.

## Controle KRAS selvagem (9BG9) — validação mecanística

Para testar se a estabilidade do tri-complexo depende da mutação, montei o
**mesmo tri-complexo com KRAS selvagem (WT, Gly13)** a partir do cristal 9BG9,
usando **protocolo idêntico** ao do mutante (mesmo campo de força, mesma
parametrização de ligantes, mesma caixa/íons, mesmos MDPs, produção de 2 ns).
A comparação controlada mostra:

| Métrica | G13D (9BG5) | WT (9BG9) | unidade |
|---|---|---|---|
| RMSD backbone (média 2ª metade) | 4.21 | **1.76** | Å |
| RMSF Cα (média) | 3.82 | **0.95** | Å |
| RMSF Cα (máx) | 7.05 | **3.09** | Å |
| Fármaco engaja KRAS | 100 | 100 | % dos frames |
| Fármaco engaja CypA | 100 | 100 | % dos frames |

O tri-complexo **WT é substancialmente mais estável** sob dinâmica, enquanto o
fármaco mantém contato simultâneo com KRAS e CypA em **100 % dos frames em
ambos os sistemas** — confirmando que a arquitetura de cola molecular
(RAS(ON)·fármaco·CypA) se forma independentemente do status G13. Arquivos em
`07_controle_9BG9/` (sistema, trajetória, figuras e tabela comparativa).

## Estrutura do repositório

```
01_selecao_pdbs/  estruturas cristalográficas + justificativa da seleção
02_ligante/       daraxonrasib: SDF, PDBQT (rígido/flexível), 2D, ligantes do co-cristal
03_receptores/    receptores preparados (PDBQT) + complexos + caixas de docking
04_docking/       scores (CSV), logs Vina, poses (PDBQT)
05_gromacs/       sistema MD pronto (tar), MDPs, complexo de partida
06_analise/       fingerprint de interações (CSV)
07_controle_9BG9/ controle KRAS WT: sistema+trajetória MD, figuras, tabela comparativa
figuras/          figuras comparativas de afinidade/contatos e mapa de resíduos
```

## Ferramentas

AutoDock Vina · Meeko · RDKit · ProLIF · MDAnalysis · OpenMM/pdbfixer ·
acpype/antechamber (GAFF2) · GROMACS 2024.1 · gemmi · Open Babel.

---
*Estudo computacional exploratório. As afinidades de docking e a trajetória de
MD são estimativas de modelagem, não substituem validação experimental.*

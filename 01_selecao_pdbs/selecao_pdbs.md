# Seleção de estruturas PDB — estudo da daraxonrasib na KRAS G13D

## 1. O mecanismo determina a seleção

A **daraxonrasib (RMC-6236)** é um inibidor **tri-complexo não-covalente, RAS(ON)
multi-seletivo** ("molecular glue"). Ela **não se liga ao KRAS isolado**: primeiro
associa-se à **ciclofilina A (CypA / PPIA)** e o binário CypA·fármaco engaja a
região **switch-I / switch-II do RAS no estado ativo (ligado a GTP)**. O sítio de
ligação é uma **superfície composta** que só existe quando CypA e RAS estão juntos.

Consequência para o desenho: docking/MD da droga só é mecanisticamente válido
**na presença de CypA + KRAS no estado ativo**. KRAS isolado (apo/GDP) serve como
**controle**, não como receptor produtivo.

## 2. Verificação do conteúdo real dos PDBs

| PDB  | Sistema | Nucleotídeo | CypA | Daraxonrasib (A1AHB) | Res 13 | Resol. |
|------|---------|:-----------:|:----:|:--------------------:|:------:|:------:|
| **9BG5** | KRAS G13D + CypA + droga (2 tri-complexos na ASU) | **GppNHp (GNP, ativo)** | ✅ | ✅ (2×) | **ASP (G13D)** | 1.67 Å |
| **9BG9** | KRAS **WT** + CypA + droga | GppNHp (ativo) | ✅ | ✅ | GLY (WT) | 1.58 Å |
| 4TQA | KRAS G13D isolado | **GDP (inativo)** | ❌ | ❌ | ASP (G13D) | 1.13 Å |
| 8BLR | KRAS4b G13D (switch-I aberto) | **GDP (inativo)** | ❌ | ❌ | ASP (G13D) | 1.40 Å |
| 9BGA | KRAS G12C + CypA + droga | GppNHp | ✅ | ✅ | — | 1.41 Å |
| 9BG6 | KRAS G12V + CypA + droga | GppNHp | ✅ | ✅ | — | 1.66 Å |

A daraxonrasib no 9BG5 é o ligante de código **A1AHB** (C44H58N8O5S, MW ≈ 811 Da,
CAS 2765081-21-6).

## 3. Mapa de contatos do co-cristal 9BG5 (ground-truth)

Distância de corte 4.5 Å, por unidade tri-complexa (KRAS-A + CypA-D + droga;
KRAS-B + CypA-C + droga):

- **CypA (14 resíduos — sítio PPIase):** Arg55, Ile57, Phe60, Met61, Gln63, Gly72,
  Ala101, Asn102, Gln111, Phe113, Trp121, Leu122, His126, Arg148.
- **KRAS switch-I (4):** Tyr32, Pro34, Thr35, Ile36.
- **KRAS switch-II (4):** Ala59, Gln61, Tyr64, Met67.

Ou seja, a droga faz de ponte entre o sítio PPIase da CypA e as regiões switch do
KRAS — a assinatura do inibidor tri-complexo.

## 4. Papel de cada estrutura no estudo

- **9BG5 — sistema-âncora / controle positivo.** Co-cristal G13D. Fonte da pose
  ground-truth para o redocking de validação e ponto de partida da MD.
- **9BG9 — controle da mutação (WT).** Mesmo farmacóforo, KRAS sem mutação →
  isola o efeito específico da G13D nos contatos.
- **4TQA — controle negativo / estado inativo.** KRAS G13D apo, GDP, sem CypA.
  Demonstra que sem CypA não há bolso produtivo; também representa o estado OFF.
- **8BLR — controle conformacional.** G13D com switch-I em conformação aberta,
  estado GDP; útil para contrastar acessibilidade do switch.
- **9BGA / 9BG6 — contexto de hotspot (opcional).** G12C / G12V para situar a G13D
  frente aos hotspots do códon 12.

Estratégia adicional de docking: reconstruir um "tri-complexo-modelo" superpondo a
CypA do 9BG5 sobre 4TQA/8BLR, permitindo comparar G13D no estado GTP (9BG5) vs GDP
(4TQA/8BLR) sob o mesmo farmacóforo composto.

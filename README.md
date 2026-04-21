# WaveletBootstrap

> Quasi U-statistics for wavelet variance estimation; bootstrap inference; applications: standard error, CIs, changepoint detection, clustering of stationary time series.

This repository contains the full source code and text for the PhD thesis:
**"Quase U-estatísticas: Representação, Classificação e Detecção de Mudanças"**
*(Quasi U-statistics: Representation, Classification and Change-point Detection)*

*   **Author:** Hilário Fernandes de Araujo Júnior
*   **Advisor:** Aluísio de Souza Pinheiro
*   **Institution:** IMECC – Unicamp, 2024
*   **Thesis PDF:** [View via Repositório Unicamp](https://hdl.handle.net/20.500.12733/41334)

## 📂 Repository Structure

The project is divided into two primary environments:

```
WaveletBootstrap/
├── src/           # R scripts used to generate simulations and results
└── Tese/          # Full LaTeX source for the final thesis document
```

## 🛠️ How to Run (`src/`)

The implementation heavily utilizes core time series and wavelet functionalities.
**Required R packages:** `ltsa`, `fracdiff`, `waveslim`, `multitaper`, `fossil`, `TSclust`, `cluster`, `dplyr`.

### ⚠️ IMPORTANT: Path Setup
Every single R script requires a variable called `BASE_PATH` to be correctly pointed to this cloned repository's root directory on your local machine.

```r
# Example:
BASE_PATH <- "C:/Users/YourUser/Projects/WaveletBootstrap"
```
Place this definition at the top of any script before you execute it, ensuring all imports and `source()` calls resolve successfully.

### Script Dependency Map

Scripts in the `src/` directory are numerically serialized. Some define reusable capabilities, while others apply them towards specific simulation tests or real-world datasets:
*   `1_*.R`, `2_*.R`, `7_*.R`, `10_*.R`: Foundational **library scripts**. Run these to load models, bootstrap resamplers, and quasi-U methods.
*   The remaining scripts are heavily parameterized simulations (`*_sim.R`, `*_resampling.R`) or studies plotting true inputs (`*_aplicacoes.R`). Note that these study components rely directly on the function scopes defined in the core libraries.

## 📊 Data Availability

**Certain datasets used in this thesis are NOT INCLUDED in this repository.** This is due to proprietary licensures. However, all public equivalent sources and structural descriptors are properly credited within the chapters of the text itself. In `_aplicacoes.R` files, comments denote where custom files were inherently pulled.

## 📜 Citation

If you use code from this repository or cite the thesis directly, please use the following BibTeX format:

```bibtex
@phdthesis{araujo2024quase,
  title  = {Quase U-estatísticas: Representação, Classificação e Detecção de Mudanças},
  author = {de Araujo Júnior, Hilário Fernandes},
  year   = {2024},
  school = {Instituto de Matemática, Estatística e Computação Científica, Universidade Estadual de Campinas}
}
```

## ⚖️ License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

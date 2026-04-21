# WaveletBootstrap

> Quasi U-statistics for wavelet variance estimation; bootstrap inference; applications: standard error, CIs, changepoint detection, clustering of stationary time series.

This repository contains the full source code and research methodology for the PhD thesis:
**"Quase U-estatísticas: Representação, Classificação e Detecção de Mudanças"**
*(Quasi U-statistics: Representation, Classification and Change-point Detection)*

*   **Author:** Hilário Fernandes de Araujo Júnior
*   **Advisor:** Aluísio de Souza Pinheiro
*   **Institution:** IMECC – Unicamp, 2024
*   **Thesis PDF:** [View via Repositório Unicamp](https://hdl.handle.net/20.500.12733/41334)

---

## 📂 Repository Structure

The project is organized into the following core directories:

- `src/`: Standardized R scripts for simulations, library functions, and real-data applications.
- `Dados/`: (Optional) Raw data directory. *Note: Some proprietary datasets are excluded for licensing reasons.*

## 🛠️ Getting Started

### Prerequisites
The implementation relies on several R packages for time series analysis and wavelets:
```r
install.packages(c("ltsa", "fracdiff", "waveslim", "multitaper", "fossil", "TSclust", "cluster", "dplyr", "imputeTS", "grDevices"))
```

### ⚠️ Critical Configuration: `BASE_PATH`
To ensure file I/O works seamlessly across different machines, every script in `src/` uses a global `BASE_PATH` variable. 

**Before running any script, set this variable at the top of the file to point to your local clone of this repository:**

```r
# Example:
BASE_PATH <- "C:/Users/YourUser/Projects/WaveletBootstrap"
```

---

## 📖 Script Guide (`src/`)

The R scripts are numerically serialized according to the thesis chapters and logical dependencies:

| Category | Scripts | Description |
| :--- | :--- | :--- |
| **Libraries** | `1_`, `2_`, `7_`, `10_` | Foundational functions (Simulations, Bootstrap methods, Quasi-U machinery). |
| **Validation** | `3_`, `4_`, `5_` | Studies on bootstrap ordering, standard errors, and confidence intervals. |
| **Analysis** | `6_`, `8_`, `9_` | Characteristic scales, resampling validation, and model comparisons. |
| **Advanced** | `11_`, `12_` | Changepoint detection and time series clustering algorithms. |
| **Utilities** | `0_` | General plotting and data visualization tools. |

---

## 📊 Data Availability

Certain datasets (*SPX Returns*) are not included due to license restrictions. Please refer to the thesis text for detailed data sources and methodology for reconstruction. Publicly available datasets are referenced within the respective `*_aplicacoes.R` scripts.

## 📜 Citation

If you use this code or cite the thesis research, please use:

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

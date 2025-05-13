# Dispersive Effective Medium for Anti-Plane Dynamics of Piezo-Magneto-Elastic Multiphase Composites

This repository contains MATLAB scripts and COMSOL simulation files used to generate Figures 3–12 and Table 3 in the accepted manuscript:

**“Dispersive Effective Medium for Anti-Plane Dynamics of Piezo-Magneto-Elastic Multiphase Composites with Generalized Periodicity”**  
Authors: *Mriganka Shekhar Chaki, Julián Bravo-Castillero, and David Guinovart*  
Published in: *Proceedings of the Royal Society A*

---

## 📚 Citation

If you use this code or data in your own work, please cite our paper:

> Chaki, M.S., Bravo-Castillero, J., & Guinovart, D. (2025). Dispersive Effective Medium for Anti-Plane Dynamics of Piezo-Magneto-Elastic Multiphase Composites with Generalized Periodicity. *Proceedings of the Royal Society A*. [DOI Link to be added]

---

## 🗂️ Folder and File Mapping for Figures and Tables

| Content         | Source Folder/File Path                                          | Notes |
|----------------|------------------------------------------------------------------|-------|
| **Fig. 3**      | `Vibration_First_order_FDM/main_file.m` (ε = 0)                   | |
| **Fig. 4**      | `Vibration_First_order_FDM/main_file.m` (ε = 0.05)                | |
| **Fig. 5**      | `Vibration_First_order_FDM/main_file.m` (compare ε = 0 and 0.05)  | |
| **Figs. 6(a–c)**| `Vibration_First_order_FDM/main_file.m` (vary H, L, Vf)           | |
| **Figs. 7–8**   | `COMSOL_piezoelectric/40lamina_wavy.mph`                          | |
| **Fig. 9**      | `Validation_COMSOL/First_order_Piezoelectric_FDM/main_file.m`     | |
| **Table 3**     | FEM (COMSOL):                                                    | |
|                | `COMSOL_piezoelectric/20lamina_wavy.mph`                          | |
|                | `COMSOL_piezoelectric/40lamina_wavy.mph`                          | |
|                | `COMSOL_piezoelectric/80lamina_wavy.mph`                          | |
|                | Analytical: `Vibration_First_order_FDM/main_file.m` (ε = 0.1, 0.05, 0.025) | |
| **Fig. 10(a)**  | `Wave_Dispersion/dispersion_2D_slowness/main_file.m`             | |
| **Figs. 10(b–c)**| `Wave_Dispersion/dispersion_2D_slowness_x2/main_file.m`         | |
| **Fig. 11(a)**  | `Wave_Dispersion/dispersion_2D_position/main_file.m`             | |
| **Fig. 11(b)**  | `Wave_Dispersion/dispersion_2D_position_surf/main_file.m`        | |
| **Figs. 12(a–d)**| `Wave_Dispersion/dispersion_parametric_effect/main_file.m` (vary ε, H, L, Vf) | |

---

## 🛠️ Software Requirements

- **MATLAB** (tested on R2022a and above recommended)
- **COMSOL Multiphysics** (version 6.2)
  - Modules: *Structural Mechanics Module*

---

## ▶️ Running the Code

Each folder contains its own `main_file.m`. Follow these general steps:

1. Open MATLAB.
2. Navigate to the corresponding folder.
3. Open and run `main_file.m`.
4. Modify parameters in the script to reproduce the exact figures or table entries, as noted above.

For COMSOL `.mph` files:
1. Open COMSOL Multiphysics 6.2.
2. Load the respective `.mph` file.
3. Run the simulation and export results as required.

---

## 📬 Contact

For questions, please contact:  
📧 [guino001@umn.edu](mailto:guino001@umn.edu)

---

## 📄 License

This project is shared for academic and research purposes. Please contact the authors for reuse outside the scope of the publication.


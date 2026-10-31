Sources for the various materials data are set out below, by file, property, material and RRR (if relevant).
The "..._fixedat20k.csv" files are copies of the main files, with properties fixed at the 20K values, for the 2D diffusion test.
Full references are shown at the bottom of the file.

#### material_properties.csv

These values cover temperatures from 0-300K [0,1,2,...,300], and are formatted with 3 header rows for (property, material, RRR). Outside the range of the source data, values have been extrapolated. Copper resitivity for RRR other than 150 has been derived from thermal conductivity using Wiedemann-Franz and constant Lorenz number - this is a limitation and the default RRR of 150 uses directly measured resistivity.

| Property | Material            | Source       | Processing                                                                                       |
|----------|---------------------|--------------|--------------------------------------------------------------------------------------------------|
| k        | Copper (RRR 50-500) | [1]          | NIST fit evaluated on 1K grid                                                                    |
| k        | Silver (RRR 30)     | [2]          | Digitised                                                                                        |
| k        | Hastelloy C-276     |              | Estimated values - replacement with measured data in progress                                    |
| k        | YBCO (in-plane)     | [3]          | Above 70K: digitised. Below 70K: estimated values - replacement with   measured data in progress |
| k        | Kapton              | [1]          | Digitised                                                                                        |
| ρ        | Copper (RRR 150)    | [4] Fig. 8.2 | Digitised                                                                                        |
| ρ        | Copper (other RRR)  | [1]          | Derived from k via Wiedemann Franz, , L₀ = 2.44×10⁻⁸ WΩ/K²                                       |
| ρ        | Silver (RRR 30)     | [2]          | Digitised                                                                                        |
| ρ        | Hastelloy C-276     |              | Estimated values - replacement with measured data in progress                                    |
| ρ        | YBCO                | [3]          | Digitised                                                                                        |
| c_p      | Copper              | [1]          | NIST fit evaluated on 1K grid                                                                    |
| c_p      | Silver              | [2]          | Digitised                                                                                        |
| c_p      | Hastelloy C-276     |              | Estimated values - replacement with measured data in progress                                    |
| c_p      | YBCO                | [5]          | Converted to J / (kg K) using assumed function unit mass of 1.106×10-24 kg                       |
| c_p      | Kapton              | [1]          | Digitised                                                                                        |

#### material_properties_magres.csv

Log (B x RRR) (in the range [0.1, 5.0], B in tesla):
- Copper: [4] pp. 8-23
- Silver: [2]

#### material_properties_caxis.csv

These values cover temperatures from 0-300K [0,1,2,...,300], and are formatted with 3 header rows for (property, material, RRR).
They show the gradient and intercept of straight lines (for each temperature) fitted to data for (depth x thermal resistivity), with thicknesses of 250, 500 and 1000 nm.

Depth x thermal resistivity slope / intercept (d_x_rho_th_slope / d_x_rho_th_intercept):
- YBCO (out-of-plane): [6] Table 3

#### Sources
[1] NIST Cryogenics Technologies Group, *Material Properties* database: OFHC Copper and Polyimide (Kapton) pages.
    https://trc.nist.gov/cryogenics/materials/materialproperties.htm (accessed October 2026).

[2] D. R. Smith and F. R. Fickett, *Cryogenic Properties of Silver*, National Institute of Standards and Technology, Boulder, CO. (1994)
    https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1363.pdf

[3] M. Matusiak, K. Rogacki and B. W. Veal, "Enhancement of the Hall-Lorenz number in optimally doped YBa₂Cu₃O₇₋δ",
    *EPL* 88, 47005 (2009). https://doi.org/10.1209/0295-5075/88/47005 (preprint: arXiv:0911.0851)

[4] N. J. Simon, E. S. Drexler and R. P. Reed, *Properties of Copper and Copper Alloys at Cryogenic Temperatures*,
    NIST Monograph 177, National Institute of Standards and Technology, 1992.

[5] Specific heat of YBCO, digitised from a figure on ResearchGate ("Variation of specific heat with temperature for
    YBCO and LSCO"): https://www.researchgate.net/figure/ariation-of-specific-heat-with-temperature-for-YBCO-and-LSCO_fig1_316831649
    Original publication not yet identified - to be replaced with a primary source.

[6] Y. Murakami, H. Goto, Y. Taguchi, Y. Nagasaka et al., "Measurement of Out-of-Plane Thermal Conductivity of
    Epitaxial YBa₂Cu₃O₇₋δ Thin Films in the Temperature Range from 10 K to 300 K by Photothermal Reflectance",
    *International Journal of Thermophysics* 38, 160 (2017). https://doi.org/10.1007/s10765-017-2294-7
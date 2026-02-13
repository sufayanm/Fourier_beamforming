# Delay-and-Sum Consistent Fourier beamforming

[![DOI](https://zenodo.org/badge/1119511884.svg)](https://doi.org/10.5281/zenodo.18633756)

This repository contains a MATLAB implementation of the image reconstruction methods used in:

S. Mulani, M. S. Ziksari, A. Austeng, and S. P. Näsholm, **“Delay-and-Sum Consistent Wavenumber Algorithm,”** IEEE Transactions on Ultrasonics, 2026 (under review).

---

## Background

The Wavenumber Algorithm (WA), introduced by Hunter et al. [1], reconstructs images from full-matrix capture (FMC) datasets using Fourier-domain migration. Compared to conventional Delay-and-Sum (DAS) beamforming, WA reduces the computational cost by approximately a factor equal to the number of transducer elements.
However, as shown analytically in our paper, WA introduces an inherent amplitude weighting when expressed in the spatial domain. Using the stationary phase method, WA can be written asymptotically as DAS multiplied by a spatial filter.

As a result, WA images typically exhibit:

- Reduced brightness at larger imaging depths  
- Distorted amplitude at higher steering angles  
- Resolution changes due to spectral weighting introduced by the filter  


The **DAS-Consistent Wavenumber Algorithm (DCWA)** removes the amplitude distortion introduced by WA while preserving its computational efficiency.

DCWA modifies the WA formulation by:

1. Multiplying the Fourier-domain data by a correction factor  
2. Scaling the reconstructed image by the axial coordinate \( z \)

Through asymptotic analysis, DCWA is shown to converge to conventional DAS in the far-field limit of the transducer element. Thus, DCWA provides a practical, low-cost Fourier-domain alternative to DAS for FMC imaging without sacrificing amplitude fidelity.

---

## Repository Structure

The file `main_file.m` runs all implemented beamforming algorithms and displays the reconstructed images. The `Algorithms/` folder contains the image reconstruction methods (WA, DAS, and DCWA), while the `Functions/` folder includes supporting functions.

The repository also includes a phantom dataset acquired using a Verasonics ultrasound system and a CIRS tissue-mimicking phantom.

---

## Citation

If you use this code in your research or build upon this work, please cite the following paper:

S. Mulani, M. S. Ziksari, A. Austeng, and S. P. Näsholm, “Delay-and-Sum Consistent Wavenumber Algorithm,” IEEE Transactions on Ultrasonics, 2026 (under review).


### BibTeX

```bibtex
@article{Mulani2025DCWA,
  author  = {Mulani, Sufayan and Sotoodeh Ziksari, Mahsa and Austeng, Andreas and Näsholm, Sven Peter},
  title   = {Delay-and-Sum Consistent Wavenumber Algorithm},
  journal = {IEEE Transactions on Ultrasonics},
  year    = {2026},
  note = {Under review}
}
```

---
## Reference
[1] A. J. Hunter, B. W. Drinkwater, and P. D. Wilcox, “The wavenumber algorithm for full-matrix imaging using an ultrasonic array,” IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, vol. 55, no. 11, pp. 2450–2462, 2008.




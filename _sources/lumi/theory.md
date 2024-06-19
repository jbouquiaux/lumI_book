---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Theory
This chapter presents the theoretical formalism behind the computation the luminescence properties of point defects in solids. It starts by a section dedicated solely to the theory of luminescence, showing how the Fermi's golden rule is approximated to obtain a tractable way of computing the luminescence lineshape, following an effective phonon mode model, or a multi-phonon mode model. Then, the computational methodology used to obtain the parameters entering the lineshape equation is presented. Finally, subtilities regarding the increase of the phonon supercell size are presented. 

## From the Fermi's Golden rule to the luminescence lineshape
Following Fermi's golden rule, the absolute luminescence intensity $I(\hbar\omega)$ (number of photons per unit time per unit energy) associated to one emitting center with two states e and g, is expressed as a function of the photon energy $\hbar\omega$ as

$$
    I(\hbar \omega)=\frac{n_D \omega^3}{3 \pi \hbar \epsilon_0 c^3} |\boldsymbol{\mu_{eg}}|^2 \delta\left(E_g-E_e-\hbar \omega\right),
$$ (first_lum_intensity)

where $n_D$ is the material refractive index, $\boldsymbol{\mu_{eg}}=\langle\Psi_{e}|\boldsymbol{\mu}|\Psi_{g}\rangle$ is the total dipole matrix element. It is then assumed that the electronic part of the dipole matrix element between the excited and ground state depends weakly on the nuclear coordinates (Franck-Condon approximation). We also suppose that the nuclear motions can be expressed as a superposition of 3N harmonic normal modes of vibration $\nu$, represented by normal coordinates $Q_{\nu}$ and frequency $\omega_{\nu}$, with N being the number of atoms in the cell considered. The total vibrational state $\chi_{\boldsymbol{n}}$ is the product of 3N independent harmonic oscillator eigenfunctions $\chi_{\boldsymbol{n_{\nu}}}$ with $n_{\nu}$ the vibrational state of the $\nu$-th harmonic oscillator. We write the normalized luminescence intensity as

$$
L(\hbar\omega)=C\omega^3A(\hbar\omega),
$$ (norm_lum_intensity)

with C a normalization constant. In the case of the absolute intensity $I(\hbar\omega)=C\omega^3A(\hbar\omega)$, $C=\frac{n_D}{3 \pi \hbar \epsilon_0 c^3}|\boldsymbol{\mu}_{eg}^{el}|^2$ with $\boldsymbol{\mu}_{eg}^{el}$ the purely electronic dipole moment (inversely proportional to the radiative lifetime of the emitting center). The emission lineshape function $ A(\hbar \omega)$, which dictates the shape of the spectrum, writes

$$
    A(\hbar \omega)= \sum_{\boldsymbol{n,m}}p_{m_{\nu}}(T)|\langle\chi_{e,\boldsymbol{m}}|\chi_{g,\boldsymbol{n}}\rangle|^2 \delta\left(E_{ZPL}-E_{g,\boldsymbol{n}}+E_{e,\boldsymbol{m}}-\hbar \omega\right).
$$ (second_lum_intensity)

In equation {eq}`second_lum_intensity`, $\boldsymbol{n}$ denotes the set of 3N vibrational states {$n_1$,$n_2$,...,$n_{3N}$}, $E_{ZPL}$ is the so-called zero-phonon line energy (purely electronic energy difference), $E_{g/e,\boldsymbol{n/m}}=\sum{_\nu}n_{\nu}\hbar\omega_{g/e,\nu}$ is the vibrational energy of state $\chi_{g/e,\boldsymbol{n/m}}$, and $p_{_{\nu}}(T)$ is a Bose-Einsten thermal occupation factor. 

The vibratinal modes of the excited and ground states can be in principle different, making equation {eq}`second_lum_intensity` computationally very expensive because of the highly multidimensional integrals to be evaluated (Franck-condon overlaps). Some approximations should be made. We assume first that the vibrational modes in the excited electronic state are identical to those in the ground state. We assume also first that we work at T=0K, to alleviate the equations ($p_{m_{\nu}}(T=0)$ = 0 if $m_{\nu}$ a vibrational state of the electronic excited state is different from 0). Considering temperature dependent luminescence spectrum will be presented later. 

### Effective phonon mode model

The simplest model to start with is built when considering that the 3N phonon modes of the system can be reduced to a single effective phonon mode (how to compute this effective phonon is discussed in section ...). Suppose that a single effective mode couples to the electronic transition, the Franck-Condon overlaps will reduce to $|\langle\chi_{e,0}|\chi_{g,n}\rangle|^2$. In the harmonic approximation, the vibrational eigenfunctions $\chi$ are expressed with Hermite polynomials and their overlaps can be computed analytically : 

$$
|\langle\chi_{e,0}|\chi_{g,n}\rangle|^2=e^{-S}\frac{S^n}{n!},
$$ (HuangRhys_eq)
with $S$ the Huang-Rhys factor, a dimensionless constant giving the mean number of effective phonons of frequency $\Omega$ involved in the electronic transition,

$$
S=\frac{\frac12\Omega^2\Delta Q^2}{\hbar\omega}=\frac{\Omega\Delta Q^2}{2\hbar},
$$ (HuangRhys_factor)

and $\Delta Q$ the offset between the harmonic potential energy surface of the ground and excited state. In practice, such offset is computed as a total mass weighted atomic displacements induced by the electronic transition : 

$$
(\Delta Q)^2=\sum_{\kappa i}m_\kappa(R_{\text{e};\kappa i}-R_{\text{g};\kappa i})^2,
$$ (Delta_Q)

where i labels the cartesian axes, $\kappa$ the atoms, $m_{\kappa}$ atomic masses, $R_{e;\kappa i}$ and $R_{g;\kappa i}$ are respectively atomic positions in the excited and the ground states. 

This allows to simplify equation {eq}`second_lum_intensity` to 

$$
A(\hbar\omega)=\sum_n e^{-S}\frac{S^n}{n!}\delta(E_{\mathrm{ZPL}}-n\hbar\Omega-\hbar\omega).
$$ (lineshape_1D)

A schematic figure of such model is presented on figure {numref}`displaced_HA` for a single mode $\nu$ (which might be the effective phonon mode of this section...)

```{figure} images/displaced_HA.png
---
height: 300px
name: displaced_HA
---
Schematic representation of the origin of photoluminescence (PL) spectra. On the left, ground state (GS) and excited state (ES) energy curves are projected along phonon normal coordinate $Q_{\nu}$ (here only for one mode) and are approximated by harmonic functions with the same frequency $\omega_{\nu}$. Vibrational energy levels and corresponding eigenfunctions are shown as horizontal lines and colored areas. GS and ES are displaced by $\Delta Q_{\nu}$ , the mass-weighted displacement between the minimum of the GS and the ES curves. On the right, the PL spectrum is formed. The zero-phonon line (ZPL) comes from the transition between the first vibrational level of the ES to the first vibrational level of the GS. Other transitions give the phonon sideband (PSB). The intensity of each peak is computed with the overlap between corresponding eigenfunctions
```

In order to visualise the effect of harmonic oscillator frequency or the offset $\Delta Q$ on the luminescence lineshape, you can slide the cursors below. For the frequency, we assume harmonic oscillators with frequency $\omega=\sqrt{\frac{k}{\mu}}$ with $k$ a force/spring constant and $\mu$ the oscillator mass.  


```{code-cell}
from lineshape_slider import get_plotly_lineshape
%matplotlib inline
```
```{code-cell}
get_plotly_lineshape(slider="frequency",num_step=10)
```

```{code-cell}
get_plotly_lineshape(slider="Delta_Q",num_step=10)
```

### Multi phonon mode model and generating function

We now consider the 3N vibrational modes. Quantities of the above section are now computed "per phonon mode". The $\Delta Q_{\nu}$ (mass-weighted atomic displacement projected along phonon mode $\nu$)

$$
\Delta Q_\nu=\sum_{\kappa\alpha}\sqrt{m_\kappa}\Delta R_{\kappa\alpha}e_{\nu,\kappa\alpha}
$$ (Delta_Q_multi)

are used to compute the partial Huang-Rhys factor $S_\nu=\frac{\omega_\nu\Delta Q_\nu^2}{2\hbar}$. $e_{\nu,\kappa\alpha}$ are the phonon eigenvectors. When dealing with 3N phonons, a direct evaluation of equation {eq}`second_lum_intensity` is impractical. Instead, the so-called generating function approach is used. The lineshape function $A(\hbar\omega)$ is evaluated as the Fourier transform of the generating function $G(t)$

$$
A(\hbar\omega,T)=\int_{-\infty}^{+\infty}G(t,T)e^{i\omega t-\frac\gamma\hbar|t|-i\frac{E^\mathrm{ZPL}}\hbar t}dt,
$$ (lineshape_generating)
with

$$
G(t)=e^{S(t)-S(0)},
$$ (generating_fct)

where $S(t)=\sum_\nu S_\nu e^{i\omega_\nu t}$ is the Fourier transform of the Huang-Rhys spectral function. $\gamma$ is a constant homogeneous Lorentzian broadening associated to each vibronic transition.

While working in the time domain, the effect of temperature (transitions involving initial vibrational states $n_{\nu} \ne 0$, weighted by Bose-Einstein occupation) can be included straightforwadly by rewriting the generating function

$$
G(t,T)=e^{S(t)-S(0)+C(t,T)+C(-t,T)-2C(0,T)}
$$ (generating_fct_T)

with $C(t,T)=\sum_\nu\overline{n}_\nu(T)S_\nu e^{i\omega_\nu t}$ the Fourier transform of the temperature weighted Huang-Rhys spectral function and $\overline{n}_\nu(T)$ the average occupation number of $\nu$-th phonon mode:

$$
\overline{n}_\nu(T)=\frac{1}{e^{\frac{\hbar\omega_\nu}{k_BT}}-1}.
$$ (bose_einstein)

### Semi-classical approach 
Formula of the fwhm

$$
W(0)=S_{\mathrm{em}}\hbar\Omega_{\mathrm{g}}\sqrt{8\ln2}/\sqrt{S_{\mathrm{abs}}}.\\
W(T)=W(0)\sqrt{\coth(\hbar\Omega_\mathrm{e}/2k_BT)},
$$ (fwhm_semi_classical)

## Computational methodology 
ToWrite

### Delta SCF constrained occupation method 
ToWrite

### Consequences of using plane-wave periodic cells DFT. 
ToWrite

### Forces vs Displacements, increase of supercell size.  
For a given phonon mode $\nu$, the corresponding partial Huang-Rhys factor $S_{\nu}=\frac{\omega_{\nu}\Delta Q_{\nu}^2}{2\hbar}$ is computed thanks to the mass-weighted displacement between ground and excited states projected on this phonon mode:

$$
	\Delta Q_\nu=\sum_{\kappa\alpha}\sqrt{M_{\kappa}}\Delta R_{\kappa\alpha}{e_{\nu,\kappa\alpha}},
$$ (Delta_Q_nu_emb)

where $\kappa$ labels the atoms and $\alpha$ the cartesian direction. Under the harmonic approximation,

$$
	M_{\kappa}\omega_{\nu}^2e_{\nu,\kappa\alpha}=\sum_{\kappa'}C_{\kappa\alpha,\kappa'\alpha'}e_{\nu,\kappa\alpha},
$$ (IFCs)

with the interatomic force constants (IFCs) $C_{\kappa\alpha,\kappa'\alpha'}=\frac{\Delta F_{\kappa'\alpha'}}{\Delta R_{\kappa\alpha}}$, Eq. {eq}`Delta_Q_nu_emb` rewrites:

$$
	\Delta Q_\nu=\frac{1}{\omega_\nu^2}\sum_{\kappa\alpha}\frac{\Delta F_{\kappa\alpha} e_{\nu,\kappa\alpha}}{\sqrt{M_{\kappa}}}.
$$ (Delta_Q_nu_emb_forces)

The use of Eq.{eq}`Delta_Q_nu_emb_forces` allows increasing the supercell size and hence minimize finite-size effects. 

Indeed, the forces decay faster to zero compared to the displacements with respect to the distance from the defect as illustrated in see {numref}`Forces_vs_dis`, panel **a-b**. 

We assume that the forces computed within the small red supercell (here 288 atoms), see {numref}`Forces_vs_dis`, panel **c** , are already converged and without finite-size effect because of this rapid decay, and that the forces outside this supercell are essentially zero. 

This means that we can use the forces computed in this red supercell, and deduce the displacements at much larger distances in the blue supercell via the IFCs computed within this large supercell. 

This subtle point is actually directly included when using Eq.{eq}`Delta_Q_nu_emb_forces` in the blue supercell. This also means that the IFCs/phonons modes should be computed in the same blue supercell. This is discussed in the next section. 

```{figure} images/Forces_vs_dis.png
---
height: 300px 
name: Forces_vs_dis
---
(a) Norm of the displacement induced by an electronic transition as a function of the distance from the defect atom (Eu). (b) Norm of the forces in the ground state at the excited state equilibrium atomic positions. The decay of forces (b) is much faster than the decay of displacements (a). (c) Cartoon of a small red supercell containing substitutional defect where the forces are the ones computed with DFT. Outside this red supercell, the forces are set to zero because of their short-range decay while the displacements are non-zero and are computed with the IFCs computed in this large blue supercell. 

```

### Approximations on the phonons computation :
- size sc deltaSCF = size sc phonons bulk (use of forces or displacements)
- size sc deltaSCF = size sc phonons with defect (use of forces or displacements)
- size sc deltaSCF < size sc phonons bulk (use of forces, no embedding)
- size sc deltaSCF < size sc phonons defect (use of forces, embedding)

### Convergence of the Huang-Rhys spectral function
ToWrite


### Suggested references 
ToWrite





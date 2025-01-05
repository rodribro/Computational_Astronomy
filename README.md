# <b> Computational_Astronomy </b>

Code for classes and 3 projects.

### <b> Project 1: Polytropic Indices </b> 

* <b> Goal: </b> Find the polytrope (of index n) that best models a real star, using the Lane-Emden equations, assuming that energy production occurs through hydrogen fusion in the core, and predict a star's density, pressure, luminosity and emmisivity throughout its radius relatively to the Sun.

* <b> Stars used:</b> Epsilon Eridani and Theta Persei A

* <b> Tools: </b> Python (NumPy, SciPy and MatPlotLib)

* <b> Folders and relevant files: </b>
    * classes: 00_integrals_ode up to 04_monte_carlo_*
    * Project1

<p> </p>
<div style="display: flex; justify-content: space-between;">
  <img src="Project1/plots_project1/lane_em_index.png" alt="Lane-Emden solutions" width="18%">
  <img src="Project1/plots_project1/em_tot.png" alt="Total Emissivity by star" width="18%">
  <img src="Project1/plots_project1/temperatures.png" alt="Temperature by star" width="18%">
  <img src="Project1/plots_project1/pressures.png" alt="Pressure by star" width="18%">
  <img src="Project1/plots_project1/lum_mass_plot.png" alt="Luminosity by star" width="18%">
</div>


### <b> Project 2: Stellar Paremeters </b> 

* <b> Goal: </b>  Predict stellar parameters based on an observed spectrum through the comparison of multiple simulated spectra to find the best possible fit.

* <b> Spectra used: </b> 2 provided observed spectrums of 2 stars and hundreds of generated spectra based on [Pollux's database](https://pollux.oreme.org/).

* <b> Tools: </b> Python (NumPy, SciPy, Pandas and MatPlotLib)

* <b> Folders and relevant files: </b>
    * classes: 05_spectral_analysis up to 08_stellar_params
    * Project2





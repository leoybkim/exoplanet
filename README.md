# Exploratory Data Analysis on NASA Exoplanet Archive

## Overview

Hi there!

Before diving in, I'd like to clarify that while I'm passionate about learning and expanding my knowledge of astronomy, I'm not an expert in this field.

Although I have made an effort to ensure the accuracy of my analysis, there may be errors or omissions, and my interpretations of the data may not be entirely accurate.

If you have any questions or feedback, please feel free to reach out to me.

I hope you find my notebook informative and enjoyable!

## Quick Start

```bash
# Create virtual environement
python -m venv venv/exoplanet

# Activate virtual environment
source venv/exoplanet/bin/activate

# Install dependencies
pip install -r requirements.txt

# Create kernel for jupyter lab
python -m ipykernel install --user --name exoplanet --display-name "exoplanet (venv)"

# Start jupyter lab
jupyter-lab
```


## NASA Exoplanet Archive

Dataset source: https://exoplanetarchive.ipac.caltech.edu/docs/counts_detail.html

Dataset type: Planetary Systems Composite Parameters

Dataset: [exoplanet_dataset.csv](exoplanet_dataset.csv)

Dataset exported on **2023-05-13**

## Dataset Column Definitions

| Column Name       | Table Label                                                 | Description                                                   | 
|:------------------|:------------------------------------------------------------|:--------------------------------------------------------------|
|pl_name|            Planet Name|                                                  Planet name most commonly used in the literature|
|hostname|           Host Name|                                                    Stellar name most commonly used in the literature|
|sy_snum|            Number of Stars|                                              Number of stars in the planetary system|                        
|sy_pnum|            Number of Planets|                                            Number of planets in the planetary system|
|sy_mnum|            Number of Moons|                                              Number of moons in the planetary system|
|cb_flag|            Circumbinary Flag|                                            Flag indicating whether the planet orbits a binary system (1=yes, 0=no)|
|discoverymethod|    Discovery Method|                                             Method by which the planet was first identified|
|disc_year|          Discovery Year|                                               Year the planet was discovered|
|disc_locale|        Discovery Locale|                                             Location of observation of planet discovery (Ground or Space)|
|disc_facility|      Discovery Facility|                                           Name of facility of planet discovery observations|
|disc_telescope|     Discovery Telescope|                                          Name of telescope of planet discovery observations|
|disc_instrument|    Discovery Instrument|                                         Name of instrument of planet discovery observations|
|rv_flag|            Detected by Radial Velocity Variations|                       Detections by Radial Velocity Variations Flag indicating if the planet host star exhibits radial velocity variations due to the planet (1=yes, 0=no)|                  
|pul_flag|           Detected by Pulsar Timing Variations|                         Boolean flag indicating if the planet host star exhibits pulsar timing variations due to the planet (1=yes, 0=no)|         
|ptv_flag|           Detected by Pulsation Timing Variations|                      Boolean flag indicating if the planet host star exhibits pulsation timing variations due to the planet (1=yes, 0=no)|            
|tran_flag|          Detected by Transits|                                         Flag indicating if the planet transits its host star (1=yes, 0=no)|
|ast_flag|           Detected by Astrometric Variations|                           Flag indicating if the planet host star exhibits astrometrical variations due to the planet (1=yes, 0=no)| 
|obm_flag|           Detected by Orbital Brightness Modulations|                   Flag indicating whether the planet exhibits orbital modulations on the phase curve (1=yes, 0=no)|      
|micro_flag|         Detected by Microlensing|                                     Boolean flag indicating if the planetary system acted as a lens during an observed microlensing event (1=yes, 0=no)|
|etv_flag|           Detected by Eclipse Timing Variations|                        Boolean flag indicating if the planet induces transit timing variations on the orbit of another another planet in the system (1=yes, 0=no)|          
|ima_flag|           Detected by Imaging|                                          Flag indicating if the planet has been observed via imaging techniques (1=yes, 0=no)|
|dkin_flag|          Detected by Disk Kinematics|                                  Boolean flag indicating if the presence of the planet was inferred due to its kinematic influence on the protoplanetary disk of its host star (1=yes, 0=no)|
|pl_orbper|          Orbital Period [days]|                                        Time the planet takes to make a complete orbit around the host star or system|
|pl_orbsmax|         Orbit Semi-Major Axis [au])|                                  The longest radius of an elliptic orbit, or, for exoplanets detected via gravitational microlensing or direct imaging, the projected separation in the plane of the sky|  
|pl_rade|            Planet Radius [Earth Radius]|                                 Length of a line segment from the center of the planet to its surface, measured in units of radius of the Earth| 
|pl_radj|            Planet Radius [Jupiter Radius]|                               Length of a line segment from the center of the planet to its surface, measured in units of radius of Jupiter|
|pl_bmasse|          Planet Mass or Mass\*sin(i) [Earth Mass]|                     Best planet mass estimate available, in order of preference: Mass, M\*sin(i)/sin(i), or M\*sin(i), depending on availability, and measured in Earth masses|   
|pl_bmassj|          Planet Mass or Mass\*sin(i) [Jupiter Mass]|                   Best planet mass estimate available, in order of preference: Mass, M\*sin(i)/sin(i), or M\*sin(i), depending on availability, and measured in Jupiter masses|   
|pl_bmassprov|       Planet Mass or Mass\*sin(i) Provenance|                       Provenance of the measurement of the best mass. Options are: Mass, M\*sin(i)/sin(i), and M\*sin(i)|
|pl_dens|            Planet Density [g/cm\*\*3]|                                   Amount of mass per unit of volume of the planet|  
|pl_orbeccen|        Eccentricity|                                                 Amount by which the orbit of the planet deviates from a perfect circle| 
|pl_insol|           Insolation Flux [Earth Flux]|                                 Insolation flux is another way to give the equilibrium temperature. It's given in units relative to those measured for the Earth from the Sun.|        
|pl_eqt|             Equilibrium Temperature [K]|                                  The equilibrium temperature of the planet as modeled by a black body heated only by its host star, or for directly imaged planets, the effective temperature of the planet required to match the measured luminosity if the planet were a black body|    
|pl_orbincl|         Inclination [deg]|                                            Angle of the plane of the orbit relative to the plane perpendicular to the line-of-sight from Earth to the object|
|pl_tranmid|         Transit Midpoint [days]|                                      The time given by the average of the time the planet begins to cross the stellar limb and the time the planet finishes crossing the stellar limb.|
|ttv_flag|           Data show Transit Timing Variations|                          Flag indicating if the planet orbit exhibits transit timing variations from another planet in the system (1=yes, 0=no). Note: Non-transiting planets discovered via the transit timing variations of another planet in the system will not have their TTV flag set, since they do not themselves demonstrate TTVs.|
|pl_imppar|          Impact Parameter|                                             The sky-projected distance between the center of the stellar disc and the center of the planet disc at conjunction, normalized by the stellar radius|
|pl_trandep|         Transit Depth [%]|                                            The size of the relative flux decrement caused by the orbiting body transiting in front of the star| 
|pl_trandur|         Transit Duration [hours]|                                     The length of time from the moment the planet begins to cross the stellar limb to the moment the planet finishes crossing the stellar limb|
|pl_ratdor|          Ratio of Semi-Major Axis to Stellar Radius|                   The distance between the planet and the star at mid-transit divided by the stellar radius. For the case of zero orbital eccentricity, the distance at mid-transit is the semi-major axis of the planetary orbit.|                                     
|pl_ratror|          Ratio of Planet to Stellar Radius|                            The planet radius divided by the stellar radius|                    
|pl_occdep|          Occultation Depth [%]|                                        Depth of occultation of secondary eclipse|
|pl_orbtper|         Epoch of Periastron [days]|                                   The time of the planet's periastron passage|  
|pl_orblper|         Argument of Periastron [deg]|                                 The angular separation between the orbit's ascending node and periastron. Note: there are a varying conventions in the exoplanet literature regarding argument of periastron (or periapsis). For example, some publications refer the planet's orbit, others to the host star's reflex orbit, which differs by 180 deg. The values in the Exoplanet Archive are not corrected to a standardized system, but are as-reported for each publication.|  
|pl_rvamp|           Radial Velocity Amplitude [m/s]|                              Half the peak-to-peak amplitude of variability in the stellar radial velocity|            
|pl_projobliq|       Projected Obliquity [deg]|                                    The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet, projected into the plane of the sky. Depending on the choice of coordinate system, projected obliquity is represented in the literature as either lambda (λ) or beta (β), where λ is defined as the negative of β (i.e., λ = -β). Since λ is reported more often than β, all values of β have been converted to λ.|
|pl_trueobliq|       True Obliquity [deg]|                                         The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet|
|st_spectype|        Spectral Type|                                                Classification of the star based on their spectral characteristics following the Morgan-Keenan system| 
|st_teff|            Stellar Effective Temperature [K]|                            Temperature of the star as modeled by a black body emitting the same total amount of electromagnetic radiation|               
|st_rad|             Stellar Radius [Solar Radius]|                                Length of a line segment from the center of the star to its surface, measured in units of radius of the Sun|  
|st_mass|            Stellar Mass [Solar mass]|                                    Amount of matter contained in the star, measured in units of masses of the Sun|
|st_met|             Stellar Metallicity [dex]|                                    Measurement of the metal content of the photosphere of the star as compared to the hydrogen content| 
|st_metratio|        Stellar Metallicity Ratio|                                    Ratio for the Metallicity Value ([Fe/H] denotes iron abundance, [M/H] refers to a general metal content)|
|st_lum|             Stellar Luminosity [log(Solar)]|                              Amount of energy emitted by a star per unit time, measured in units of solar luminosities|  
|st_logg|            Stellar Surface Gravity [log10(cm/s\*\*2)]|                   Gravitational acceleration experienced at the stellar surface|                      
|st_age|             Stellar Age [Gyr]|                                            The age of the host star| 
|st_dens|            Stellar Density [g/cm\*\*3]|                                  Amount of mass per unit of volume of the star|
|st_vsin|            Stellar Rotational Velocity [km/s]|                           Rotational velocity at the equator of the star multiplied by the sine of the inclination|                
|st_rotp|            Stellar Rotational Period [days]|                             The time required for the planet host star to complete one rotation, assuming it is a solid body|          
|st_radv|            Systemic Radial Velocity [km/s]|                              Velocity of the star in the direction of the line of sight|          
|ra|                 RA [deg]|                                                     Right Ascension of the planetary system in decimal degrees|
|dec|                Dec [deg]|                                                    Declination of the planetary system in decimal degrees|
|glat|               Galactic Latitude [deg]|                                      Galactic latitude of the planetary system in units of decimal degrees|
|glon|               Galactic Longitude [deg]|                                     Galactic longitude of the planetary system in units of decimal degrees|
|elat|               Ecliptic Latitude [deg]|                                      Ecliptic latitude of the planetary system in units of decimal degrees| 
|elon|               Ecliptic Longitude [deg]|                                     Ecliptic longitude of the planetary system in units of decimal degrees|
|sy_pm|              Total Proper Motion [mas/yr]|                                 Angular change in position over time as seen from the center of mass of the Solar System|
|sy_pmra|            Proper Motion (RA) [mas/yr]|                                  Angular change in right ascension over time as seen from the center of mass of the Solar System|
|sy_pmdec|           Proper Motion (Dec) [mas/yr]|                                 Angular change in declination over time as seen from the center of mass of the Solar System|
|sy_dist|            Distance [pc]|                                                Distance to the planetary system in units of parsecs| 
|sy_plx|             Parallax [mas]|                                               Difference in the angular position of a star as measured at two opposite positions within the Earth's orbit| 
|sy_bmag|            B (Johnson) Magnitude|                                        Brightness of the host star as measured using the B (Johnson) band in units of magnitudes|
|sy_vmag|            V (Johnson) Magnitude|                                        Brightness of the host star as measured using the V (Johnson) band in units of magnitudes
|sy_jmag|            J (2MASS) Magnitude|                                          Brightness of the host star as measured using the J (2MASS) band in units of magnitudes|
|sy_hmag|            H (2MASS) Magnitude|                                          Brightness of the host star as measured using the H (2MASS) band in units of magnitudes|
|sy_kmag|            Ks (2MASS) Magnitude|                                         Brightness of the host star as measured using the K (2MASS) band in units of magnitudes|
|sy_umag|            u (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) u band, in units of magnitudes|
|sy_gmag|            g (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) g band, in units of magnitudes|
|sy_rmag|            r (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) r band, in units of magnitudes|
|sy_imag|            i (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) i band, in units of magnitudes|
|sy_zmag|            z (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) z band, in units of magnitudes|
|sy_w1mag|           W1 (WISE) Magnitude|                                          Brightness of the host star as measured using the 3.4um (WISE) band in units of magnitudes.|
|sy_w2mag|           W2 (WISE) Magnitude|                                          Brightness of the host star as measured using the 4.6um (WISE) band in units of magnitudes.|
|sy_w3mag|           W3 (WISE) Magnitude|                                          Brightness of the host star as measured using the 12.um (WISE) band in units of magnitudes|
|sy_w4mag|           W4 (WISE) Magnitude|                                          Brightness of the host star as measured using the 22.um (WISE) band in units of magnitudes|
|sy_gaiamag|         Gaia Magnitude|                                               Brightness of the host star as measuring using the Gaia band in units of magnitudes. Objects matched to Gaia using the Hipparcos or 2MASS IDs provided in Gaia DR2.|
|sy_icmag|           I (Cousins) Magnitude|                                        Brightness of the host star as measured using the I (Cousins) band in units of magnitudes.|
|sy_tmag|            TESS Magnitude|                                               Brightness of the host star as measured using the TESS bandpass, in units of magnitudes|
|sy_kepmag|          Kepler Magnitude|                                             Brightness of the host star as measured using the Kepler bandpass, in units of magnitudes|


## NASA Exoplanet Catalog

Data source: [NASA's exoplanet catalog](https://exoplanets.nasa.gov/discovery/exoplanet-catalog/)

Webscraper script: [web_scraper.py](/web_scraping/web_scraper.py)

Dataset: [nasa_exoplanet_catalog.csv](/web_scrapgin/nasa_exoplanet_catalog.csv)

Data web-scraped on **2023-05-13**




<ins>References:</ins>
- https://exoplanets.nasa.gov/
- https://en.wikipedia.org/wiki/Orbital_elements
- https://en.wikipedia.org/wiki/Exoplanet_orbital_and_physical_parameters
- https://en.wikipedia.org/wiki/Stellar_classification
- https://en.wikipedia.org/wiki/Correlation_coefficient
- https://www.nasa.gov/kepler/overview/abouttransits
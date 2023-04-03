# NASA Exoplanet Archive

Dataset source: https://exoplanetarchive.ipac.caltech.edu/docs/counts_detail.html

Dataset type: Planetary Systems Composite Parameters

Dataset exported on **2023-04-01**

As of export date, there are currently **5322** confirmed exoplanets from 3989 different planetary systems!


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
|pl_orbpererr1|      Orbital Period Upper Unc. [days]|                             |     
|pl_orbpererr2|      Orbital Period Lower Unc. [days]|                             |     
|pl_orbperlim|       Orbital Period Limit Flag|                                    |
|pl_orbsmax|         Orbit Semi-Major Axis [au])|                                  The longest radius of an elliptic orbit, or, for exoplanets detected via gravitational microlensing or direct imaging, the projected separation in the plane of the sky| 
|pl_orbsmaxerr1|     Orbit Semi-Major Axis Upper Unc. [au]|                        |          
|pl_orbsmaxerr2|     Orbit Semi-Major Axis Lower Unc. [au]|                        |          
|pl_orbsmaxlim|      Orbit Semi-Major Axis Limit Flag|                             |     
|pl_rade|            Planet Radius [Earth Radius]|                                 Length of a line segment from the center of the planet to its surface, measured in units of radius of the Earth| 
|pl_radeerr1|        Planet Radius Upper Unc. [Earth Radius]|                      |            
|pl_radeerr2|        Planet Radius Lower Unc. [Earth Radius]|                      |            
|pl_radelim|         Planet Radius Limit Flag|                                     | 
|pl_radj|            Planet Radius [Jupiter Radius]|                               Length of a line segment from the center of the planet to its surface, measured in units of radius of Jupiter|
|pl_radjerr1|        Planet Radius Upper Unc. [Jupiter Radius]|                    |
|pl_radjerr2|        Planet Radius Lower Unc. [Jupiter Radius]|                    |
|pl_radjlim|         Planet Radius Limit Flag|                                     |
|pl_bmasse|          Planet Mass or Mass\*sin(i) [Earth Mass]|                     Best planet mass estimate available, in order of preference: Mass, M\*sin(i)/sin(i), or M\*sin(i), depending on availability, and measured in Earth masses|   
|pl_bmasseerr1|      Planet Mass or Mass\*sin(i) [Earth Mass] Upper Unc.|          |
|pl_bmasseerr2|      Planet Mass or Mass\*sin(i) [Earth Mass] Lower Unc.|          |
|pl_bmasselim|       Planet Mass or Mass\*sin(i) [Earth Mass] Limit Flag|          |
|pl_bmassj|          Planet Mass or Mass\*sin(i) [Jupiter Mass]|                   Best planet mass estimate available, in order of preference: Mass, M\*sin(i)/sin(i), or M\*sin(i), depending on availability, and measured in Jupiter masses|   
|pl_bmassjerr1|      Planet Mass or Mass\*sin(i) [Jupiter Mass] Upper Unc.|        |
|pl_bmassjerr2|      Planet Mass or Mass\*sin(i) [Jupiter Mass] Lower Unc.|        |
|pl_bmassjlim|       Planet Mass or Mass\*sin(i) [Jupiter Mass] Limit Flag|        |
|pl_bmassprov|       Planet Mass or Mass\*sin(i) Provenance|                       Provenance of the measurement of the best mass. Options are: Mass, M\*sin(i)/sin(i), and M\*sin(i)|
|pl_dens|            Planet Density [g/cm\*\*3]|                                   Amount of mass per unit of volume of the planet|  
|pl_denserr1|        Planet Density Upper Unc. [g/cm\*\*3]|                        |          
|pl_denserr2|        Planet Density Lower Unc. [g/cm\*\*3]|                        |         
|pl_denslim|         Planet Density Limit Flag|                                    |
|pl_orbeccen|        Eccentricity|                                                 Amount by which the orbit of the planet deviates from a perfect circle| 
|pl_orbeccenerr1|    Eccentricity Upper Unc.|                                      |
|pl_orbeccenerr2|    Eccentricity Lower Unc.|                                      |
|pl_orbeccenlim|     Eccentricity Limit Flag|                                      |
|pl_insol|           Insolation Flux [Earth Flux]|                                 Insolation flux is another way to give the equilibrium temperature. It's given in units relative to those measured for the Earth from the Sun.|        
|pl_insolerr1|       Insolation Flux Upper Unc. [Earth Flux]|                      |           
|pl_insolerr2|       Insolation Flux Lower Unc. [Earth Flux]|                      |           
|pl_insollim|        Insolation Flux Limit Flag|                                   |
|pl_eqt|             Equilibrium Temperature [K]|                                  The equilibrium temperature of the planet as modeled by a black body heated only by its host star, or for directly imaged planets, the effective temperature of the planet required to match the measured luminosity if the planet were a black body|  
|pl_eqterr1|         Equilibrium Temperature Upper Unc. [K]|                       |             
|pl_eqterr2|         Equilibrium Temperature Lower Unc. [K]|                       |          
|pl_eqtlim|          Equilibrium Temperature Limit Flag|                           |      
|pl_orbincl|         Inclination [deg]|                                            Angle of the plane of the orbit relative to the plane perpendicular to the line-of-sight from Earth to the object|
|pl_orbinclerr1|     Inclination Upper Unc. [deg]|                                 |                  
|pl_orbinclerr2|     Inclination Lower Unc. [deg]|                                 |
|pl_orbincllim|      Inclination Limit Flag|                                       |
|pl_tranmid|         Transit Midpoint [days]|                                      The time given by the average of the time the planet begins to cross the stellar limb and the time the planet finishes crossing the stellar limb.|
|pl_tranmiderr1|     Transit Midpoint Upper Unc. [days]|                           | 
|pl_tranmiderr2|     Transit Midpoint Lower Unc. [days]|                           |
|pl_tranmidlim|      Transit Midpoint Limit Flag|                                  |
|ttv_flag|           Data show Transit Timing Variations|                          Flag indicating if the planet orbit exhibits transit timing variations from another planet in the system (1=yes, 0=no). Note: Non-transiting planets discovered via the transit timing variations of another planet in the system will not have their TTV flag set, since they do not themselves demonstrate TTVs.|
|pl_imppar|          Impact Parameter|                                             The sky-projected distance between the center of the stellar disc and the center of the planet disc at conjunction, normalized by the stellar radius|
|pl_impparerr1|      Impact Parameter Upper Unc.|                                  |
|pl_impparerr2|      Impact Parameter Lower Unc.|                                  |
|pl_impparlim|       Impact Parameter Limit Flag|                                  |
|pl_trandep|         Transit Depth [%]|                                            The size of the relative flux decrement caused by the orbiting body transiting in front of the star| 
|pl_trandeperr1|     Transit Depth Upper Unc. [%]|                                 | 
|pl_trandeperr2|     Transit Depth Lower Unc. [%]|                                 | 
|pl_trandeplim|      Transit Depth Limit Flag|                                     |
|pl_trandur|         Transit Duration [hours]|                                     The length of time from the moment the planet begins to cross the stellar limb to the moment the planet finishes crossing the stellar limb|
|pl_trandurerr1|     Transit Duration Upper Unc. [hours]|                          |
|pl_trandurerr2|     Transit Duration Lower Unc. [hours]|                          |        
|pl_trandurlim|      Transit Duration Limit Flag|                                  |
|pl_ratdor|          Ratio of Semi-Major Axis to Stellar Radius|                   The distance between the planet and the star at mid-transit divided by the stellar radius. For the case of zero orbital eccentricity, the distance at mid-transit is the semi-major axis of the planetary orbit.|               
|pl_ratdorerr1|      Ratio of Semi-Major Axis to Stellar Radius Upper Unc.|        |                          
|pl_ratdorerr2|      Ratio of Semi-Major Axis to Stellar Radius Lower Unc.|        |                          
|pl_ratdorlim|       Ratio of Semi-Major Axis to Stellar Radius Limit Flag|        |                          
|pl_ratror|          Ratio of Planet to Stellar Radius|                            The planet radius divided by the stellar radius|      
|pl_ratrorerr1|      Ratio of Planet to Stellar Radius Upper Unc.|                 |
|pl_ratrorerr2|      Ratio of Planet to Stellar Radius Lower Unc.|                 |                 
|pl_ratrorlim|       Ratio of Planet to Stellar Radius Limit Flag|                 |                 
|pl_occdep|          Occultation Depth [%]|                                        Depth of occultation of secondary eclipse|
|pl_occdeperr1|      Occultation Depth Upper Unc. [%]|                             |     
|pl_occdeperr2|      Occultation Depth Lower Unc. [%]|                             |     
|pl_occdeplim|       Occultation Depth Limit Flag|                                 | 
|pl_orbtper|         Epoch of Periastron [days]|                                   The time of the planet's periastron passage| 
|pl_orbtpererr1|     Epoch of Periastron Upper Unc. [days]|                        |          
|pl_orbtpererr2|     Epoch of Periastron Lower Unc. [days]|                        |          
|pl_orbtperlim|      Epoch of Periastron Limit Flag|                               |   
|pl_orblper|         Argument of Periastron [deg]|                                 The angular separation between the orbit's ascending node and periastron. Note: there are a varying conventions in the exoplanet literature regarding argument of periastron (or periapsis). For example, some publications refer the planet's orbit, others to the host star's reflex orbit, which differs by 180 deg. The values in the Exoplanet Archive are not corrected to a standardized system, but are as-reported for each publication.| 
|pl_orblpererr1|     Argument of Periastron Upper Unc. [deg]|                      |            
|pl_orblpererr2|     Argument of Periastron Lower Unc. [deg]|                      |            
|pl_orblperlim|      Argument of Periastron Limit Flag|                            |      
|pl_rvamp|           Radial Velocity Amplitude [m/s]|                              Half the peak-to-peak amplitude of variability in the stellar radial velocity|    
|pl_rvamperr1|       Radial Velocity Amplitude Upper Unc. [m/s]|                   |               
|pl_rvamperr2|       Radial Velocity Amplitude Lower Unc. [m/s]|                   |               
|pl_rvamplim|        Radial Velocity Amplitude Limit Flag|                         |         
|pl_projobliq|       Projected Obliquity [deg]|                                    The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet, projected into the plane of the sky. Depending on the choice of coordinate system, projected obliquity is represented in the literature as either lambda (λ) or beta (β), where λ is defined as the negative of β (i.e., λ = -β). Since λ is reported more often than β, all values of β have been converted to λ.|
|pl_projobliqerr1|   Projected Obliquity Upper Unc. [deg]|                         |         
|pl_projobliqerr2|   Projected Obliquity Lower Unc. [deg]|                         |         
|pl_projobliqlim|    Projected Obliquity Limit Flag|                               |   
|pl_trueobliq|       True Obliquity [deg]|                                         The angle between the angular momentum vector of the rotation of the host star and the angular momentum vector of the orbit of the planet|
|pl_trueobliqerr1|   True Obliquity Upper Unc. [deg]|                              |    
|pl_trueobliqerr2|   True Obliquity Lower Unc. [deg]|                              |    
|pl_trueobliqlim|    True Obliquity Limit Flag|                                    | 
|st_spectype|        Spectral Type|                                                Classification of the star based on their spectral characteristics following the Morgan-Keenan system| 
|st_teff|            Stellar Effective Temperature [K]|                            Temperature of the star as modeled by a black body emitting the same total amount of electromagnetic radiation|      
|st_tefferr1|        Stellar Effective Temperature Upper Unc. [K]|                 |                 
|st_tefferr2|        Stellar Effective Temperature Lower Unc. [K]|                 |                 
|st_tefflim|         Stellar Effective Temperature Limit Flag|                     |             
|st_rad|             Stellar Radius [Solar Radius]|                                Length of a line segment from the center of the star to its surface, measured in units of radius of the Sun|  
|st_raderr1|         Stellar Radius Upper Unc. [Solar Radius]|                     |             
|st_raderr2|         Stellar Radius Lower Unc. [Solar Radius]|                     |             
|st_radlim|          Stellar Radius Limit Flag|                                    | 
|st_mass|            Stellar Mass [Solar mass]|                                    Amount of matter contained in the star, measured in units of masses of the Sun|
|st_masserr1|        Stellar Mass Upper Unc. [Solar mass]|                         |         
|st_masserr2|        Stellar Mass Lower Unc. [Solar mass]|                         |         
|st_masslim|         Stellar Mass Limit Flag|                                      |
|st_met|             Stellar Metallicity [dex]|                                    Measurement of the metal content of the photosphere of the star as compared to the hydrogen content|
|st_meterr1|         Stellar Metallicity Upper Unc. [dex]|                         |         
|st_meterr2|         Stellar Metallicity Lower Unc. [dex]|                         |         
|st_metlim|          Stellar Metallicity Limit Flag|                               |   
|st_metratio|        Stellar Metallicity Ratio|                                    Ratio for the Metallicity Value ([Fe/H] denotes iron abundance, [M/H] refers to a general metal content)|
|st_lum|             Stellar Luminosity [log(Solar)]|                              Amount of energy emitted by a star per unit time, measured in units of solar luminosities|  
|st_lumerr1|         Stellar Luminosity Upper Unc. [log(Solar)]|                   |               
|st_lumerr2|         Stellar Luminosity Lower Unc. [log(Solar)]|                   |               
|st_lumlim|          Stellar Luminosity Limit Flag|                                |  
|st_logg|            Stellar Surface Gravity [log10(cm/s\*\*2)]|                   Gravitational acceleration experienced at the stellar surface|                 
|st_loggerr1|        Stellar Surface Gravity Upper Unc. [log10(cm/s\*\*2)]|        |                          
|st_loggerr2|        Stellar Surface Gravity Lower Unc. [log10(cm/s\*\*2)]|        |                          
|st_logglim|         Stellar Surface Gravity Limit Flag|                           |       
|st_age|             Stellar Age [Gyr]|                                            The age of the host star| 
|st_ageerr1|         Stellar Age Upper Unc. [Gyr]|                                 | 
|st_ageerr2|         Stellar Age Lower Unc. [Gyr]|                                 | 
|st_agelim|          Stellar Age Limit Flag|                                       | 
|st_dens|            Stellar Density [g/cm\*\*3]|                                  Amount of mass per unit of volume of the star|
|st_denserr1|        Stellar Density Upper Unc. [g/cm\*\*3]|                       |           
|st_denserr2|        Stellar Density Lower Unc. [g/cm\*\*3]|                       |           
|st_denslim|         Stellar Density Limit Flag|                                   |
|st_vsin|            Stellar Rotational Velocity [km/s]|                           Rotational velocity at the equator of the star multiplied by the sine of the inclination|       
|st_vsinerr1|        Stellar Rotational Velocity [km/s] Upper Unc.|                |                  
|st_vsinerr2|        Stellar Rotational Velocity [km/s] Lower Unc.|                |                  
|st_vsinlim|         Stellar Rotational Velocity Limit Flag|                       |           
|st_rotp|            Stellar Rotational Period [days]|                             The time required for the planet host star to complete one rotation, assuming it is a solid body|     
|st_rotperr1|        Stellar Rotational Period [days] Upper Unc.|                  |                
|st_rotperr2|        Stellar Rotational Period [days] Lower Unc.|                  |                
|st_rotplim|         Stellar Rotational Period Limit Flag|                         |         
|st_radv|            Systemic Radial Velocity [km/s]|                              Velocity of the star in the direction of the line of sight|    
|st_radverr1|        Systemic Radial Velocity Upper Unc. [km/s]|                   |                
|st_radverr2|        Systemic Radial Velocity Lower Unc. [km/s]|                   |               
|st_radvlim|         Systemic Radial Velocity Limit Flag|                          |        
|ra|                 RA [deg]|                                                     Right Ascension of the planetary system in decimal degrees|
|dec|                Dec [deg]|                                                    Declination of the planetary system in decimal degrees|
|glat|               Galactic Latitude [deg]|                                      Galactic latitude of the planetary system in units of decimal degrees|
|glon|               Galactic Longitude [deg]|                                     Galactic longitude of the planetary system in units of decimal degrees|
|elat|               Ecliptic Latitude [deg]|                                      Ecliptic latitude of the planetary system in units of decimal degrees| 
|elon|               Ecliptic Longitude [deg]|                                     Ecliptic longitude of the planetary system in units of decimal degrees|
|sy_pm|              Total Proper Motion [mas/yr]|                                 Angular change in position over time as seen from the center of mass of the Solar System|
|sy_pmerr1|          Total Proper Motion Upper Unc [mas/yr]|                       |
|sy_pmerr2|          Total Proper Motion Lower Unc [mas/yr]|                       |
|sy_pmra|            Proper Motion (RA) [mas/yr]|                                  Angular change in right ascension over time as seen from the center of mass of the Solar System|
|sy_pmraerr1|        Proper Motion (RA) [mas/yr] Upper Unc|                        |
|sy_pmraerr2|        Proper Motion (RA) [mas/yr] Lower Unc|                        |
|sy_pmdec|           Proper Motion (Dec) [mas/yr]|                                 Angular change in declination over time as seen from the center of mass of the Solar System|
|sy_pmdecerr1|       Proper Motion (Dec) [mas/yr] Upper Unc|                       |
|sy_pmdecerr2|       Proper Motion (Dec) [mas/yr] Lower Unc|                       |
|sy_dist|            Distance [pc]|                                                Distance to the planetary system in units of parsecs| 
|sy_disterr1|        Distance [pc] Upper Unc|                                      |
|sy_disterr2|        Distance [pc] Lower Unc|                                      |
|sy_plx|             Parallax [mas]|                                               Difference in the angular position of a star as measured at two opposite positions within the Earth's orbit| 
|sy_plxerr1|         Parallax [mas] Upper Unc|                                     |
|sy_plxerr2|         Parallax [mas] Lower Unc|                                     |
|sy_bmag|            B (Johnson) Magnitude|                                        Brightness of the host star as measured using the B (Johnson) band in units of magnitudes|
|sy_bmagerr1|        B (Johnson) Magnitude Upper Unc|                              |
|sy_bmagerr2|        B (Johnson) Magnitude Lower Unc|                              |
|sy_vmag|            V (Johnson) Magnitude|                                        Brightness of the host star as measured using the V (Johnson) band in units of magnitudes
|sy_vmagerr1|        V (Johnson) Magnitude Upper Unc|                              |
|sy_vmagerr2|        V (Johnson) Magnitude Lower Unc|                              |
|sy_jmag|            J (2MASS) Magnitude|                                          Brightness of the host star as measured using the J (2MASS) band in units of magnitudes|
|sy_jmagerr1|        J (2MASS) Magnitude Upper Unc|                                |
|sy_jmagerr2|        J (2MASS) Magnitude Lower Unc|                                |
|sy_hmag|            H (2MASS) Magnitude|                                          Brightness of the host star as measured using the H (2MASS) band in units of magnitudes|
|sy_hmagerr1|        H (2MASS) Magnitude Upper Unc|                                |
|sy_hmagerr2|        H (2MASS) Magnitude Lower Unc|                                |
|sy_kmag|            Ks (2MASS) Magnitude|                                         Brightness of the host star as measured using the K (2MASS) band in units of magnitudes|
|sy_kmagerr1|        Ks (2MASS) Magnitude Upper Unc|                               |
|sy_kmagerr2|        Ks (2MASS) Magnitude Lower Unc|                               |
|sy_umag|            u (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) u band, in units of magnitudes|
|sy_umagerr1|        u (Sloan) Magnitude Upper Unc|                                |
|sy_umagerr2|        u (Sloan) Magnitude Lower Unc|                                |
|sy_gmag|            g (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) g band, in units of magnitudes|
|sy_gmagerr1|        g (Sloan) Magnitude Upper Unc|                                |
|sy_gmagerr2|        g (Sloan) Magnitude Lower Unc|                                |
|sy_rmag|            r (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) r band, in units of magnitudes|
|sy_rmagerr1|        r (Sloan) Magnitude Upper Unc|                                |
|sy_rmagerr2|        r (Sloan) Magnitude Lower Unc|                                |
|sy_imag|            i (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) i band, in units of magnitudes|
|sy_imagerr1|        i (Sloan) Magnitude Upper Unc|                                |
|sy_imagerr2|        i (Sloan) Magnitude Lower Unc|                                |
|sy_zmag|            z (Sloan) Magnitude|                                          Brightness of the host star as measured using the Sloan Digital Sky Survey (SDSS) z band, in units of magnitudes|
|sy_zmagerr1|        z (Sloan) Magnitude Upper Unc|                                |
|sy_zmagerr2|        z (Sloan) Magnitude Lower Unc|                                |
|sy_w1mag|           W1 (WISE) Magnitude|                                          Brightness of the host star as measured using the 3.4um (WISE) band in units of magnitudes.|
|sy_w1magerr1|       W1 (WISE) Magnitude Upper Unc|                                |
|sy_w1magerr2|       W1 (WISE) Magnitude Lower Unc|                                |
|sy_w2mag|           W2 (WISE) Magnitude|                                          Brightness of the host star as measured using the 4.6um (WISE) band in units of magnitudes.|
|sy_w2magerr1|       W2 (WISE) Magnitude Upper Unc|                                |
|sy_w2magerr2|       W2 (WISE) Magnitude Lower Unc|                                |
|sy_w3mag|           W3 (WISE) Magnitude|                                          Brightness of the host star as measured using the 12.um (WISE) band in units of magnitudes|
|sy_w3magerr1|       W3 (WISE) Magnitude Upper Unc|                                |
|sy_w3magerr2|       W3 (WISE) Magnitude Lower Unc|                                |
|sy_w4mag|           W4 (WISE) Magnitude|                                          Brightness of the host star as measured using the 22.um (WISE) band in units of magnitudes|
|sy_w4magerr1|       W4 (WISE) Magnitude Upper Unc|                                |
|sy_w4magerr2|       W4 (WISE) Magnitude Lower Unc|                                |
|sy_gaiamag|         Gaia Magnitude|                                               Brightness of the host star as measuring using the Gaia band in units of magnitudes. Objects matched to Gaia using the Hipparcos or 2MASS IDs provided in Gaia DR2.|
|sy_gaiamagerr1|     Gaia Magnitude Upper Unc|                                     |
|sy_gaiamagerr2|     Gaia Magnitude Lower Unc|                                     |
|sy_icmag|           I (Cousins) Magnitude|                                        Brightness of the host star as measured using the I (Cousins) band in units of magnitudes.|
|sy_icmagerr1|       I (Cousins) Magnitude Upper Unc|                              |
|sy_icmagerr2|       I (Cousins) Magnitude Lower Unc|                              | 
|sy_tmag|            TESS Magnitude|                                               Brightness of the host star as measured using the TESS bandpass, in units of magnitudes|
|sy_tmagerr1|        TESS Magnitude Upper Unc|                                     |
|sy_tmagerr2|        TESS Magnitude Lower Unc|                                     |
|sy_kepmag|          Kepler Magnitude|                                             Brightness of the host star as measured using the Kepler bandpass, in units of magnitudes|
|sy_kepmagerr1|      Kepler Magnitude Upper Unc|                                   |
|sy_kepmagerr2|      Kepler Magnitude Lower Unc|                                   |




<ins>References:</ins>
- https://exoplanets.nasa.gov/
- https://en.wikipedia.org/wiki/Orbital_elements
- https://en.wikipedia.org/wiki/Exoplanet_orbital_and_physical_parameters
- https://en.wikipedia.org/wiki/Stellar_classification
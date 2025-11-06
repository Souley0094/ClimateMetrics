# Climate metrics!
#
# This package aim to make avaible somes climates metrics such as GWP (global warming potential) and AGTP (absolute global temperature change)
# of CO2 N2O and CH4, and albedo change based on Sieber et al. 2020 DOI:10.1111/gcbb.12682.


if(!require(dbplyr)){install.packages("dplyr");  library(dplyr)}
if(!require(dbplyr)){install.packages("ggplot2");  library(ggplot2)}
if(!require(dbplyr)){install.packages("cowplot");  library(cowplot)}

#' Global mean temperature Change of albedo and GHG
#'
#' @param RF a value of a vector of values of albedo or GHG radiaive forcings
#' @param duration  An integer guiding the  number of years for duration of emission scenario
#' @param effect A string representing the effect e.g "albedo", "CO2" for carbon dioxyde etc...
#'
#' @returns A vector of values corresponding to the AGTP for a time horizon defined
#' @export
#'
#' @examples AGTP(RF = -3, duration = 30, effect = "albedo")
#' @examples AGTP(RF = -3, duration = 30, effect = "N2O")

AGTP <- function(RF,duration, effect){

  Ag <-  510064472*10**6 # Earth surface area (m2)

  A <-  10000 # Functional unit 1 ha = 10,000 m2

  Area <- Ag/A # convert into global RF W/ha

  TH <- 100 # time horizon, i.e. evaluation period

  H <- 0:TH # years for IRF


  if(effect == "albedo"){
    AGTP_alb <- matrix(0, nrow = TH + 1, ncol = 2)

    # Boucher & Reddy constants for temperature response

    c <- c(0.631, 0.429) # cj, components of climate sensitivity

    d <- c(8.4, 409.5) # dj, short and long response timescale

    for (j in c(0, 1)) {
      AGTP_alb[, j + 1] <- c[j + 1] * ((1 - exp(-H / d[j + 1])) -  (1 - exp(-(H - 1) / d[j + 1])) * ifelse(H - 1 >= 0, 1, 0)) # Heaviside step function
    }

    # Sum across the columns to combine j components

    AGTP_alb <- rowSums(AGTP_alb)

    # Ensure RF is a vector

    temp_alb <- function(RF) {
      # Pad RF with zeros to length TH longer
      RF <- c(RF, rep(0, TH))

      # Convolution (np.convolve with mode="full")
      dT <- convolve(RF, rev(AGTP_alb), type = "open")

      # Crop to TH+1 values (R is 1-based indexing)
      return(dT[1:(TH + 1)])
    }

    dT <- temp_alb(rep(RF,30))

    return(dT)

  }else if(effect == "CO2"){

    # Radiative efficiency per kg gas, including indirect effects
    RE_CO2 <- 1.75171717997627e-15

    # Metric values for GHG decay functions
    # Metric values for CO2 (AR5 based on Joos et al., 2013)

    a0 <- 0.2173
    a1 <- 0.2240
    a2 <- 0.2824
    a3 <- 0.2763 # Decay coefficients

    tao1 <- 394.4
    tao2 <- 36.54
    tao3 <- 4.304 # Decay time scales

    a <- c(a0, a1, a2, a3)  # Combine decay coefficients into a vector

    tao <- c(0, tao1, tao2, tao3)  # Combine decay time scales into a vector

    c <- c(0.631, 0.429) # cj, components of climate sensitivity

    d <- c(8.4, 409.5) # dj, short and long response timescale


    # Initialize AGTP_CO2 as a matrix of zeros

    AGTP_CO2 <- matrix(0, nrow = TH + 1, ncol = 8)

    # Loop over j
    for (j in c(0, 1)) {

      AGTP_CO2[, j + 1] <- a[1] * c[j + 1] * (1 - exp(-H / d[j + 1]))

      # Loop over i
      for (i in c(1, 2, 3)) {

        AGTP_CO2[, j * 3 + 1 + i] <- a[i + 1] * tao[i + 1] * c[j + 1] / (tao[i + 1] - d[j + 1]) * (exp(-H / tao[i + 1]) - exp(-H / d[j + 1]))
      }
    }

    # Sum across j and i components, then multiply with RE_CO2
    AGTP_CO2 <- RE_CO2 * rowSums(AGTP_CO2)

    # Convolution of emission vector (annual) and AGTP = analytical solution of integral
    temp_CO2 <- function(E) {
      # Pad E with zeros to match time horizon
      E <- c(E, rep(0, TH))

      # Convolution (equivalent to np.convolve(mode='full'))
      dT <- convolve(E, rev(AGTP_CO2), type = "open")

      # Crop to TH+1 elements
      return(dT[1:(TH + 1)])
    }

    dT_CO2 <- temp_CO2(rep(RF,duration))

    return(dT_CO2)

  }else if(effect == "N2O"){

    # Radiative efficiency per kg gas, including indirect effects

    RE_N2O <- 3.56791115325733e-13

    # Metric values for N2O

    taoN2O <- 121

    c <- c(0.631, 0.429) # cj, components of climate sensitivity

    d <- c(8.4, 409.5) # dj, short and long response timescale


    # Initialize AGTP_CO2 as a matrix of zeros

    AGTP_N2O <- matrix(0, nrow = length(H), ncol = 2)

    for (j in 0:1) {
      # R is 1-based
      AGTP_N2O[, j+1] <- taoN2O * c[j+1] / (taoN2O - d[j+1]) *
        (exp(-H / taoN2O) - exp(-H / d[j+1]))
    }

    # Sum across j and i components, then multiply with RE_CO2
    AGTP_N2O <- RE_N2O * rowSums(AGTP_N2O)

    temp_N2O <- function(E) {
      # Pad E with zeros to match time horizon
      E <- c(E, rep(0, TH))

      # Convolution (equivalent to np.convolve(mode='full'))
      dT <- convolve(E, rev(AGTP_N2O), type = "open")

      # Crop to TH+1 elements
      return(dT[1:(TH + 1)])
    }

    dT_N2O <- temp_N2O(rep(RF,duration))

    return(dT_N2O)

  }else if(effect == "CH4"){

    # Radiative efficiency per kg gas, including indirect effects

    RE_CH4 <- 2.10762811889549e-13

    # Metric values for N2O

    taoCH4 <- 12.4

    c <- c(0.631, 0.429) # cj, components of climate sensitivity

    d <- c(8.4, 409.5) # dj, short and long response timescale

    # Initialize AGTP_CO2 as a matrix of zeros

    AGTP_CH4 <- matrix(0, nrow = length(H), ncol = 2)

    for (j in 0:1) {
      # R is 1-based
      AGTP_CH4[, j+1] <- taoCH4 * c[j+1] / (taoCH4 - d[j+1]) *
        (exp(-H / taoCH4) - exp(-H / d[j+1]))
    }


    # Sum across j and i components, then multiply with RE_CO2
    AGTP_CH4 <- RE_CH4 * rowSums(AGTP_CH4)

    temp_CH4 <- function(E) {
      # Pad E with zeros to match time horizon
      E <- c(E, rep(0, TH))

      # Convolution (equivalent to np.convolve(mode='full'))
      dT <- convolve(E, rev(AGTP_CH4), type = "open")

      # Crop to TH+1 elements
      return(dT[1:(TH + 1)])
    }

    dT_CH4 <- temp_CH4(rep(RF,duration))

    return(dT_CH4)

  }else{
    print("Error! put the right effect or gas")
  }


}


#Global warming potential of albedo based on Bright et al. 2016
#GWP = RF(i)/AGWP_CO2 with RF daily or Monthly radiative forcing and AGWP_CO2
# the absolute global warming potential of CO2 gas expressed in kg CO2 at time horizon (TH)



#' Global warming potential of albedo from Bright et al. 2016
#'
#' @param RF a value of a vector of values of albedo or GHG radiaive forcings
#' @param TH  An integer representing the time horizon, i.e. evaluation period
#'
#' @returns A values in Kg CO2 eq per ha integrating the GWP of albedo of a given TH e.g. 20, 100
#' @export
#'
#' @examples GWP_albedo(-3, 20)
#' @examples GWP_labedo(-3,100)

GWP_albedo <- function(RF, TH){

  Ag <-  510064472*10**6 # Earth surface area (m2)

  A <-  10000 # Functional unit 1 ha = 10,000 m2

  Area <- Ag/A # convert into global RF W/ha

  AGWP_CO2 <- function(TH) {
    # Atmospheric concentration
    c0_CO2 <- 391       # background CO2 (ppm)
    c1_CO2 <- c0_CO2 + 1

    # Constants
    m_air <- 28.97       # kg/mol
    m_CO2 <- 44.0098     # kg/mol
    m_tot <- 5.1352e18   # total mass of atmosphere (kg)

    # Convert ppm to kg CO2
    ppmv_to_kg_CO2 <- 1 / ((m_air / m_CO2) * (1e6 / m_tot))

    alpha_CO2 <- 5.35

    # Radiative forcing
    RE_CO2_ppm <- alpha_CO2 * log(c1_CO2 / c0_CO2)
    RE_CO2_kg  <- RE_CO2_ppm / ppmv_to_kg_CO2

    # CO2 decay function
    a0 <- 0.2173; a1 <- 0.2240; a2 <- 0.2824; a3 <- 0.2763
    tao1 <- 394.4; tao2 <- 36.54; tao3 <- 4.304

    # AGWP cumulative forcing
    AGWP <- RE_CO2_kg * (
      a0*TH +
        a1*tao1*(1 - exp(-TH/tao1)) +
        a2*tao2*(1 - exp(-TH/tao2)) +
        a3*tao3*(1 - exp(-TH/tao3))
    )

    #return(AGWP)
  }

  RF <- RF/Area
  GWP <- RF/AGWP_CO2(TH)

  return(GWP)
}



#' Global warming potential metric for GHG
#'
#' @param x A numeric value or vector of numeric values of a GHG emission expressed into CO2 e.g 44/12 to expressed mass of C to mass of CO2 and 44/28 to convert from Mass N to mass of N2O
#' @param IPCC_factor A value corresponding to IPCC conversion factor of a GHG e.g 1 for CO2, 28 for CH4 and 273 for N2O at 100 years time horizon (TH) refer to AR6 (IPCC; 2023)
#'
#' @returns Global warming potential values expressed in kg CO2 eq for a given TH
#' @export
#'
#' @examples GWP(-280.3,273) # GWP of N2O for 100 years TH
#' @examples GWP(-120.12,28) # GWP of CH4 for 100 years TH

GWP <- function(x, IPCC_factor){

  GWP <- x*IPCC_factor

  return(GWP)

}

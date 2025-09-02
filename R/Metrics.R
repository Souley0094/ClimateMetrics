# Climate metrics!
#
# This package aim to make avaible somes climates metrics such as GWP (global warming potential) and AGTP (absolute global temperature change)
# of CO2 N2O and CH4, and albedo change


if(!require(dbplyr)){install.packages("dplyr");  library(dplyr)}
if(!require(dbplyr)){install.packages("ggplot2");  library(ggplot2)}
if(!require(dbplyr)){install.packages("cowplot");  library(cowplot)}

AGTP <- function(RF,duration, effect){
  TH <- 100 # time horizon, i.e. evaluation period
  if(effect == "albedo"){
    AGTP_alb <- matrix(0, nrow = TH + 1, ncol = 2)

    # Boucher & Reddy constants for temperature response

    c <- c(0.631, 0.429) # cj, components of climate sensitivity

    d <- c(8.4, 409.5) # dj, short and long response timescale

    H <- seq(0, TH, 1) # years for IRF

    for (j in c(0, 1)) {
      AGTP_alb[, j + 1] <- c[j + 1] * ((1 - exp(-H / d[j + 1])) -  (1 - exp(-(H - 1) / d[j + 1])) * ifelse(H - 1 >= 0, 1, 0)) # Heaviside step function
    }

    # Sum across the columns to combine j components

    AGTP_alb <- rowSums(AGTP_alb)

    # Ensure RF is a vector

    RF <- as.numeric(RF)

    # Pad RF with zeros to match the length of H

    RF <- c(RF, rep(0, TH))

    # Convolution sum of RF and AGTP_alb

    dT <- convolve(RF, rev(AGTP_alb), type = "open")

    # Crop the array to return only the first TH + 1 elements

    dT <- dT[1:(TH + 1)]
    return(dT)

  }else if(effect == "C"){

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
    return(AGTP_CO2)

  }else if(effects == "N2O"){

    # Radiative efficiency per kg gas, including indirect effects

    RE_N2O <- 3.56791115325733e-13

    # Metric values for N2O

    taoN2O <- 121


    # Initialize AGTP_CO2 as a matrix of zeros

    AGTP_N2O <- matrix(0, nrow = TH + 1, ncol = 8)

    # Loop over j
    for (j in c(0, 1)) {

      AGTP_N2O[, j + 1] <- a[1] * c[j + 1] * (1 - exp(-H / d[j + 1]))

      # Loop over i
      for (i in c(1, 2, 3)) {

        AGTP_N2O[, j * 3 + 1 + i] <- a[i + 1] * tao[i + 1] * c[j + 1] / (tao[i + 1] - d[j + 1]) * (exp(-H / tao[i + 1]) - exp(-H / d[j + 1]))
      }
    }

    # Sum across j and i components, then multiply with RE_CO2
    AGTP_N2O <- RE_N2O * rowSums(AGTP_N2O)

    return(AGTP_N2O)

  }else{

    # Radiative efficiency per kg gas, including indirect effects

    RE_CH4 <- 2.10762811889549e-13

    # Metric values for N2O

    taoN2O <- 12.4


    # Initialize AGTP_CO2 as a matrix of zeros

    AGTP_CH4 <- matrix(0, nrow = TH + 1, ncol = 8)

    # Loop over j
    for (j in c(0, 1)) {

      AGTP_CH4[, j + 1] <- a[1] * c[j + 1] * (1 - exp(-H / d[j + 1]))

      # Loop over i
      for (i in c(1, 2, 3)) {

        AGTP_CH4[, j * 3 + 1 + i] <- a[i + 1] * tao[i + 1] * c[j + 1] / (tao[i + 1] - d[j + 1]) * (exp(-H / tao[i + 1]) - exp(-H / d[j + 1]))
      }
    }

    # Sum across j and i components, then multiply with RE_CO2
    AGTP_CH4 <- RE_CH4 * rowSums(AGTP_CH4)

    return(AGTP_CH4)

  }



}


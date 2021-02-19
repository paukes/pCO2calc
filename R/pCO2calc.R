#' Calculation of pCO2 within a water sample using Henry's Law
#'
#' This takes the CO2 using the head-space equilibirum technique and calculates water sample pCO2.
#' @param df Dataframe with all the values
#' @param temp_field_C Measured temperature of the water in the field (C)
#' @param p_field_kPa Measured barometric pressure in the field (kPa)
#' @param temp_lab_C Measured temperature in the lab during analysis (C)
#' @param p_lab_kPa Measured barometric pressure in the lab (kPa)
#' @param salt_added_g Measured weight of salt added in vial (g; can be 0 if nothing added)
#' @param vial_g Weight of an empty vial (g)
#' @param vial_full_g Weight of a vial full of water (g)
#' @param vial_HS_g Weight of the vial after headspace created (g)
#' @param he_inj_mL Volume of helium injected while sample drawn out (mL)
#' @param measCO2_ppm Measured CO2 from gas chromatograph (ppm)
#' @param atmCO2_ppm Current atmospheric concentration of CO2 at time of sampling (ppm)
#' @keywords CO2, headspace equilibrium, Henry's Law,
#' @export
#' @examples
#' lab_datasheet <- data.frame(Sample = c('Fun Lake', 'Not Fun Lake'), temp_field_C = c(25, 4),p_field_kPa = c(100, 93.5),temp_lab_C = c(22, 22),p_lab_kPa = c(101,101),KCl_added_g = c(0.05, 0.05),exe_g = c(19.8, 19.7),exe_full_g = c(26, 26),exe_headspace_g = c(24.5, 24.3),he_inj_mL = c(10, 10),measCO2_ppm = c(850, 405))
#' pCO2calc(lab_datasheet, 'temp_field_C', 'p_field_kPa', 'temp_lab_C', 'p_lab_kPa', 'KCl_added_g','exe_g', 'exe_full_g', 'exe_headspace_g', 'he_inj_mL', 'measCO2_ppm', 410)



pCO2calc <- function(df, temp_field_C, p_field_kPa, temp_lab_C, p_lab_kPa, salt_added_g,
                     vial_g, vial_full_g, vial_HS_g, he_inj_mL, measCO2_ppm, atmCO2_ppm){

  ## Check to see if using values or specific columns for each parameter ##
  temp_field_C <- if(is.numeric(temp_field_C) == TRUE) {temp_field_C} else {df[[temp_field_C]]}
  p_field_kPa <- if(is.numeric(p_field_kPa) == TRUE) {p_field_kPa} else {df[[p_field_kPa]]}
  temp_lab_C <- if(is.numeric(temp_lab_C) == TRUE) {temp_lab_C} else {df[[temp_lab_C]]}
  p_lab_kPa <- if(is.numeric(p_lab_kPa) == TRUE) {p_lab_kPa} else {df[[p_lab_kPa]]}
  salt_added_g <- if(is.numeric(salt_added_g) == TRUE) {salt_added_g} else {df[[salt_added_g]]}
  vial_g <- if(is.numeric(vial_g) == TRUE) {vial_g} else {df[[vial_g]]}
  vial_full_g <- if(is.numeric(vial_full_g) == TRUE) {vial_full_g} else {df[[vial_full_g]]}
  vial_HS_g <- if(is.numeric(vial_HS_g) == TRUE) {vial_HS_g} else {df[[vial_HS_g]]}
  he_inj_mL <- if(is.numeric(he_inj_mL) == TRUE) {he_inj_mL} else {df[[he_inj_mL]]}
  measCO2_ppm <- if(is.numeric(measCO2_ppm) == TRUE) {measCO2_ppm} else {df[[measCO2_ppm]]}
  atmCO2_ppm <- if(is.numeric(atmCO2_ppm) == TRUE) {atmCO2_ppm} else {df[[atmCO2_ppm]]}


  ## Convert units ##
  temp_field_K <- temp_field_C + 273.15;
  temp_lab_K <- temp_lab_C + 273.15;
  p_field_atm <- p_field_kPa * 0.009869;
  p_lab_atm <- p_lab_kPa * 0.009869;
  vol_water_L <- (vial_full_g - vial_g)/1000; #assuming density of 1.00 (g ~ mL)
  vol_HS_L <- (vial_full_g - vial_HS_g)/1000;
  vol_waterHS_L <- (vol_water_L - vol_HS_L);
  he_inj_L <- he_inj_mL/1000

  ## Define Constants ##
  R_gas <- 0.08205737
  #Weiss constants for Kh
  a1 <- -58.0931
  a2 <- 90.5069
  a3 <- 22.294
  b1 <- 0.027766
  b2 <- -0.025888
  b3 <- 0.0050578
  #Determine the salinity
  sal_permil <- salt_added_g / vol_water_L

  ## Calculate Henry Constants ##
  #Kh of He at lab temperature (using van't Hoff equation)
  Kh_he <- 0.00038 * exp(92 * ((1/temp_lab_K) - (1/298.15)))
  #Kh of CO2 at lab temperature and salinity (using Weiss 1974)
  #(remember that log() in R is ln() )
  Kh_co2 <- exp(a1 + (a2*(100/temp_lab_K)) + (a3 * log(temp_lab_K/100)) + (sal_permil * (b1 + (b2 * (temp_lab_K/100)) + (b3 * ((temp_lab_K/100)^2)))))
  #Kh of CO2 at field temperature (with 0 salinity in freshwaters)
  Kh_co2_field <- exp(a1 + (a2*(100/temp_field_K)) + (a3 * log(temp_field_K/100)))


  #Calculate the new pressure in the headspace in the vial
  p_He_atm <- ((p_lab_atm * he_inj_L)/(R_gas * temp_lab_K)) / ((vol_waterHS_L * Kh_he) + (vol_HS_L / (R_gas * temp_lab_K)))

  #Calculate the pCO2 in headspace
  pCO2_HS <- (measCO2_ppm / (10^6)) * p_He_atm

  #Calculate [CO2] in water with headspace
  CO2_vial_mol.L <- Kh_co2 * pCO2_HS

  #Calculate the total [CO2] in the vial from the sum of moles of CO2 in headspace and moles of CO2 in water, divided by the amount of water in the vial after headspace
  n_HS <- (pCO2_HS * vol_HS_L)/(R_gas * temp_lab_K)
  n_diss <- CO2_vial_mol.L * vol_waterHS_L
  co2_orig <- (n_HS + n_diss) / vol_waterHS_L

  #Calculate the %-saturation relative to the atmosphere
  CO2_perc_sat <- (co2_orig / (Kh_co2_field * (atmCO2_ppm * (10^-6) * p_field_atm))) * 100

  #Get units of pCO2
  pCO2_uM <- co2_orig * 10^6
  pCO2_uatm <- (co2_orig / Kh_co2_field) * 10^6

  #Export these results to the original dataframe
  df$pCO2_uM <- pCO2_uM
  df$pCO2_uatm <- pCO2_uatm
  df$pCO2_perc_sat <- CO2_perc_sat

  return(df)

}

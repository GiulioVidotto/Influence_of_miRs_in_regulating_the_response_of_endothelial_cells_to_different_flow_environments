# Function: get_Kd
#
# This function calculates the equilibrium dissociation constant (Kd) from the Gibbs free energy change (delta_G) 
# and the temperature (T_constant) using the following thermodynamic formula:
#
#  Kd = exp(-delta_G / (R * T))
#
# where:
#   - delta_G: The Gibbs free energy change in kcal/mol.
#   - T_constant: The temperature in Celsius (which is internally converted to Kelvin).
#   - R: The gas constant (1.987 x 10^-3 kcal/mol·K).
#
# Input:
#   - delta_G: The Gibbs free energy change (in kcal/mol).
#   - T_constant: The temperature in degrees Celsius.
#
# Output:
#   - Kd: The calculated equilibrium dissociation constant.

get_Kd <- function(delta_G, T_costant) {
  if (!is.numeric(delta_G) || !is.numeric(T_costant)) {
    stop("Both delta_G and T_constant must be numeric values.")
  }
  if (T_costant < -273.15) {
    stop("Temperature cannot be below absolute zero (in Celsius).")
  }
  # The measure unit of delta_G is kcal/mol 
  # Convert the temperature from Celsius to Kelvin
  T_costant <- T_costant + 273.15 # K
  # R: Gas constant
  R_costant <- 1.987e-3 # Kcal / mol·K 
  # Calculate the equilibrium constant
  Kd <- exp(-1 * (delta_G / (R_costant * T_costant)))
  return(Kd)
}
module Constants
# Module with useful global constants

# Light speed in m/s
const c = 3e8;

# Complex water refractive index (for wavelength 3,2 cm and 10 celsius degrees temperature)
const m = 7.8 + im*2.44;

# Parameter Km for water
const Km = (m^2 - 1) / (m^2 + 2);

# Square of Km abs
const Km_sq_abs = abs(Km)^2;

# Reciever noise PSD (4kT)
const N0 = 1.656e-20;

end

add_subdirectory("1-03-poisson")

add_subdirectory("2-06-bc-newton")

add_subdirectory("3-07-navier-stokes")

add_subdirectory("5-04-complex-adapt")

add_subdirectory("6-05-hcurl-adapt")

add_subdirectory("7-3-05-newton-heat-rk")

add_subdirectory("A-electrostatic-simple")

add_subdirectory("B-electrostatic-adapt")

add_subdirectory("C-subdomains")

add_subdirectory("D-subdomains-multi")

add_subdirectory("E-singlemesh")

# Commented out due to eigenproblems solver being reworked.
#add_subdirectory("8-01-eigenvalue")

# Still in the old forms.
IF(WITH_TRILINOS)
  add_subdirectory("9-02-trilinos-nonlinear")
ENDIF(WITH_TRILINOS)

add_subdirectory("linear-advection-dg-adapt")

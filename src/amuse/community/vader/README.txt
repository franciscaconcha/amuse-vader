VADER is a code simulating the evolution of viscous thin accretion
disks. It is developed by Mark Krumholz and John Forbes [1]. 

Note that the version in AMUSE is an old version of the code;
the main branch has been updated to Python 3, while AMUSE works
with Python 2. The code is available on github, with the version
used here having the identifier [2], submitted on march 5th 2017.

VADER is a very flexible code. A large number of viscous
parameters, such as the boundary conditions and mass source
function, can be implemented via specialized c code. However,
that does mean that VADER must be re-compiled if the user desires
to implement their own code. This has the following steps:

1:  Write your own functions. These must be in a file with the
	name userFunc_{PROB}.c, where {PROB} is the name of your
	problem/set of functions. It is convenient to just copy
	userFunc_none.c, and replace any empty functions with the
	desired ones. These functions will only be used if the
	corresponding function parameters (i.e. alpha_function) are
	set to true. 

2:  User-defined parameters can be passed to the userFunc through
	the array pointed to by params. In VADER, their data type can
	be freely chosen and defined, but in the interface every
	parameter must be a float/double. Before assigning, the number
	of user parameters must be defined (through the parameter
	number_of_user_parameters), and then each parameter is
	assigne	by calling set_parameter(i, X), where i is the
	parameter index and X the parameter itself. These values must
	be assigned without units. 

3:  Compile the code using the following:
	- make clean
	- make PROB={PROB}
	Leaving PROB={PROB} out compiles VADER with userFunc_none. 


Running the simulation takes a couple of steps. These include:

1:  Initialize the code. The interface passes all quantities
	in cgs units; if working in physical units, no converter is
	needed. 
	Vader has its own diagnostic outputs, which can be enabled via
	the verbosity parameter (0 through 3). Note that for this to
	show on terminal, the initialization must include
	" redirection='none' ".

2:  Call the initialize_code function.

3:  Initialize a grid, either with a Keplerian or flat rotation
	curve. Free and tabulated rotation curves are currently
	not available.

4:  Assign parameters. Default values can be found in DOCU.txt.

5:  Assign column density and pressure distributions. Note that
	the pressure is actually vertically-integrated pressure, and
	has units of [pressure]*[length]. Internal energy need only
	be assigned if equation_of_state_function is set to True.

6:  Evolve model as desired.


Details on the various parameters and (grid) properties can be 
found in DOCU.txt.

For any questions on VADER, consult the manual included with the
git repository. For any questions on the interface, contact the
coder at wilhelm@strw.leidenuniv.nl.





[1] Krumholz, M. R. and Forbes, J. C., VADER: A Flexible, Robust,
Open-Source Code for Simulating Viscous Thin Accretion Disks, 
Astronomy and Computing, Vol. 11 (2015)

[2] ef719e8b068ab101b3fdb35bead6a317a0c1ae78

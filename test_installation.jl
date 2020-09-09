using DifferentialEquations
using ReactionMechanismSimulator

file_path = "mech/superminimal.rms"; # Change it to your personal path
#load mechanism dictionary
phaseDict = readinput(file_path);
#mechanism dictionaries index:  phaseDict[phasename]["Species" or "Reactions"]
spcs = phaseDict["phase"]["Species"];
rxns = phaseDict["phase"]["Reactions"];

#######################    INPUT     #######################
V0 = 8.314 * 1000 / 1e5;  # [m^3] # Comparison with chemkin
P0 = 1.0e+5; #[Pa]
# T0 = 1000; #[K]
H2 = 0.67;  #[mol]
O2 = 0.33;  #[mol]
N0 = H2 + O2;
############################################################

# Define the phase (how species thermodynamic and kinetic properties calculated)
ig = IdealGas(spcs,rxns,name="phase");
# Define initial condition of the reactor
initialconds = Dict(["V"=>V0,"P"=>P0, "H2"=>H2,"O2"=>O2]);
# Define the domain (encodes how system thermodynamic properties calculated)
domain,y0 = ConstantPDomain(phase=ig,initialconds=initialconds);

#######################    INPUT     #######################
t_final = 2.1; #[s]
solver = DifferentialEquations.CVODE_BDF();
abstol = 1e-20;
reltol = 1e-12;
############################################################

react = Reactor(domain,y0,(0.0,t_final)); #Create the reactor object

#solve the ode associated with the reactor
sol = solve(react.ode,solver, abstol=reltol,reltol=reltol);
bsol = Simulation(sol,domain);

sol(2.0000)[domain.indexes[end]]

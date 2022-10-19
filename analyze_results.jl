using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

# path to save files
global path = "$(@__DIR__)/plots/immature_fully_healing/$(FILE)/"

# # clean old files 
# if isdir(path)
#     rm(path, recursive = true)
# end

mkpath(path)

# global out_path = "$(@__DIR__)/data/$(FILE)/"
global out_path = "$(@__DIR__)/data/immature_fully_healing/$(FILE)/"

# Global variables
yr2sec = 365*24*60*60
# comment this part if there is nothing in event_time temporarily

# Read data
event_time = readdlm(string(out_path, "event_time.out"), header=false)
tStart = event_time[:,1]
println("Start time of all seismic events(s):",tStart) 
tEnd = event_time[:,2]
println("Duration of all seismic events(s):",tEnd-tStart)


hypo = event_time[:,3]
d_hypo = event_time[:,4]    # unit: m 
println("Cumulative slips when earthquakes happen:",d_hypo) 
println("Depth of all seismic events:",hypo)

event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
indx = Int(length(event_stress[1,:])/2)

taubefore = event_stress[:,1:indx]
tauafter = event_stress[:,indx+1:end]
stressdrops = taubefore .- tauafter

# coseismic slip on fault for all different events(row)
delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
println("Dimension of cosesimic slip:",size(delfafter))
sliprate = readdlm(string(out_path, "sliprate.out"), header=false)   # every 10 timesteps
println("Dimension of sliprate:",size(sliprate))

println("Total number of all seismic events:",size(delfafter,1))
println("Total number of all on-fault GLL nodes:",size(delfafter,2))     

# displacement on fault line for different time 
delfsec = readdlm(string(out_path, "delfsec.out"))
# print(size(delfsec))
delfyr = readdlm(string(out_path, "delfyr.out"))
# print(size(delfyr))


# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)

Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
println("Dimension of FltX:",size(FltX))
cca = params[4,:]
ccb = params[5,:]
a_b = cca .- ccb
Lc = params[6,:]

# Index of fault from 0 to 18 km
flt18k = findall(FltX .<= 18)[1]

time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
t = time_vel[:,1]             # all real timsteps
Vfmax = time_vel[:,2]
Vsurface = time_vel[:,3]
alphaa = time_vel[:,4]         # initial background rigidity ratio

# stress = readdlm(string(out_path, "stress.out"), header=false)
# start_index = get_index(stress', taubefore')

# Event_details

rho1 = 2670
vs1 = 3462
rho2 = 2670
vs2 = sqrt(alphaa[1])*vs1
mu = rho2*vs2^2    # to calculate seismic moment
println("Shear modulus of damage zone:",mu)

Mw, del_sigma, fault_slip, rupture_len =
        moment_magnitude_new(mu, FltX, delfafter', stressdrops');

println("Moment magnitudes of all seismic events:", Mw)
println("Average stress drops of all seismic events(MPa):", del_sigma)
println("Average fault slips of all seismic events(m):", fault_slip)
println("Rupture lengths along depth of all seismic events(km):", rupture_len./1e3)



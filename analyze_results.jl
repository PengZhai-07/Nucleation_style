using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

# Global variables
yr2sec = 365*24*60*60

# Read data
time_vel = readdlm(string(out_path, "time_velocity.out"), header=false, Float64)
t = time_vel[:,1]         # all real timsteps
Vfmax = time_vel[:,2]
alphaa = time_vel[:,3]         # initial background rigidity ratio
#b_value = time_vel[:,5]


# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)
Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
#println("Dimension of FltX:",size(FltX))
cca = params[4,:]
ccb = params[5,:]
a_b = params[6,:]
Lc = params[7,:]

event_time = readdlm(string(out_path, "event_time.out"), header=false)
tStart = event_time[:,1]
println("Start time of all seismic events(s):",tStart) 
tEnd = event_time[:,2]
println("Duration of all seismic events(s):",tEnd-tStart)

hypo = event_time[:,3]
d_hypo = event_time[:,4]    # unit: m 
println("Cumulative slips when earthquakes happen:",d_hypo) 
println("Depth of all seismic events:",hypo)

sliprate = readdlm(string(out_path, "sliprate.out"), header=false, use_mmap=true)  
println("Dimension of sliprate:",size(sliprate))

# weakeningrate = exp.(readdlm(string(out_path, "weakeningrate.out"), header=false,use_mmap=true))
# println("Dimension of weakeningrate:",size(weakeningrate))

# coseismic slip on fault for all different events(row)
delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
println("Dimension of cosesimic slip:",size(delfafter))
N_events = size(delfafter,1)   # here the number of event should depend on the event_time.out file
println("Total number of all seismic events:", N_events)
println("Total number of all on-fault GLL nodes:",size(delfafter,2))    

# displacement on fault line for different time 
# delfsec = readdlm(string(out_path, "delfsec.out"))   # every 0.1 second, but if timestep is large enough, it should be every timestep!!
delfyr = readdlm(string(out_path, "delfyr.out"))
# print(size(delfyr))

delfsec_et = readdlm(string(out_path, "delfsec_each_timestep.out"), header=false, use_mmap=true)    # every 10 timesteps in coseismic phase
# # index_ds_start, index_ds_end = get_index_delfsec(N_events, delfsec_et)      

println(size(delfsec_et))
temp = readdlm(string(out_path, "delfsec_each_timestep_endline.out"), header=false)
index_ds_end::Vector{Int} = temp[:,1]
println(index_ds_end)
index_ds_start::Vector{Int64} = zeros(length(index_ds_end))
index_ds_start[2:end] = (index_ds_end .+ 1)[1:end-1]
index_ds_start[1] = 1
println(index_ds_start)
# here the number of effective events should depend on the event_time.out file 
# and we ignore some small events
 
event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
indx = Int(length(event_stress[1,:])/2)

taubefore = event_stress[:,1:indx]
tauafter = event_stress[:,indx+1:end]
stressdrops = taubefore .- tauafter      # unit: MPa 

stress = readdlm(string(out_path, "stress.out"), header=false)   # unit: MPa  timesteps/10, shear stress on fault line points
# get the start and end time of every seismic event
index_start, index_end = get_index(t, tStart, tEnd)         # 
println(index_start)
println(index_end)

#Event_details
G = 3e10
# comment this part if there is nothing in event_time temporarily
rho1::Float64 = 2670
vs1::Float64 = sqrt(G/rho1)
rho2 = rho1
vs2 = sqrt(alphaa[1])*vs1   # initial shear wave velocity
mu = rho2*vs2^2    # to calculate seismic moment
println("Shear modulus of damage zone:",mu)

Mw, del_sigma, fault_slip, rupture_len, scaled_energy, radiation_eff =
        moment_magnitude_new(mu/1e6, FltX, delfafter', stressdrops', delfsec_et,index_ds_start, index_ds_end, stress,index_start, index_end);   # Time*L

println("Moment magnitudes of all seismic events:", Mw)
println("Average stress drops of all seismic events(MPa):", del_sigma)
println("Average fault slips of all seismic events(m):", fault_slip)
println("Rupture lengths along depth of all seismic events(km):", rupture_len./1e3)

println("Saled energy of all seismic events:", scaled_energy)
println("Radiation efficiency of all seismic events:", radiation_eff)

open(string(path,"Scaled energy info.out"), "w") do io
        for i in eachindex(Mw)
                write(io, join(hcat(Mw[i], del_sigma[i], fault_slip[i], rupture_len[i]/1e3, 
                scaled_energy[i], radiation_eff[i]), " "), "\n") 
        end
end
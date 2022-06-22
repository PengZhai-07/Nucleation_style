using DelimitedFiles

# FILE = "1000_4_1"          # 2.43h: 8735.684543 seconds (1.31 G allocations: 253.616 GiB, 0.09% gc time, 0.15% compilation time)
FILE = "1000_6_1"          # 6.5h: 23304.773837 seconds (3.43 G allocations: 675.554 GiB, 0.10% gc time, 0.04% compilation time)
# FILE = "48000.0_1000_6_0.9"      # 6.5h: 22120.908916 seconds (2.75 G allocations: 608.486 GiB, 0.08% gc time, 0.29% compilation time)
alpha = 1.0

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")
 
# path to save files
global path = "$(@__DIR__)/plots/$(FILE)/"
mkpath(path)

global out_path = "$(@__DIR__)/data/$(FILE)/"

# Global variables
yr2sec = 365*24*60*60

# Read data
event_time = readdlm(string(out_path, "event_time.out"), header=false)
tStart = event_time[:,1]
#println(tStart)
tEnd = event_time[:,2]

hypo = event_time[:,3]
d_hypo = event_time[:,4]    # unit: m 
print(d_hypo) 
print(hypo)

event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
indx = Int(length(event_stress[1,:])/2)
# println(length(event_stress[1,:]))  # 962
# println(indx)          # 481
taubefore = event_stress[:,1:indx]
tauafter = event_stress[:,indx+1:end]

# coseismic slip on fault for all different events(row)
delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
# print(size(delfafter))
#  slip = readdlm(string(out_path, "slip.out"), header=false)
sliprate = readdlm(string(out_path, "sliprate.out"), header=false)
# println(size(delfafter,1))
# println(size(delfafter,2))

# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)

Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
cca = params[4,:]
ccb = params[5,:]
a_b = cca .- ccb
Lc = params[6,:]

# Index of fault from 0 to 18 km
flt18k = findall(FltX .<= 18)[1]

time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
t = time_vel[:,1]
Vfmax = time_vel[:,2]
Vsurface = time_vel[:,3]
alphaa = time_vel[:,4]


rho1 = 2670
vs1 = 3464
rho2 = 2500
vs2 = alpha*vs1
mu = rho2*vs2^2

# displacement on fault line for different time 
delfsec = readdlm(string(out_path, "delfsec.out"))
# print(size(delfsec))
delfyr = readdlm(string(out_path, "delfyr.out"))
# print(size(delfyr))
stress = readdlm(string(out_path, "stress.out"), header=false)


start_index = get_index(stress', taubefore')
stressdrops = taubefore .- tauafter

Mw, del_sigma, fault_slip, rupture_len =
        moment_magnitude_new(mu, FltX, delfafter', stressdrops', t);

using HyperFSI

mesh_file = joinpath(@__DIR__, "..", "Probe_in_flow.inp")
out_boundaries = ["Line1", "Line2", "Line3", "Line4", "Line5", "Line6"]

Δt = 1e-5
n_step = 10000
n_fq = 200

geo = Post2D(mesh_file, out_boundaries)
pos = geo.pos
area = geo.area
h = maximum(area)^(1/2)
r_h = h * 4
vol = h * area

#temperature_dependent_thermophysical properties
function f_kp_cc(kc, kp, T)
    kc_actul = 4.732+2.837e-2*T-2.288e-5*T^2+6.368e-9*T^3
    kp_actul = kc_actul/kc * kp
    return kp_actul
end

function f_cv_cc(cv, T)
    return -38.0528 + 0.00041618*T +741.254/T - 0.707584*√T + 19.0915*log(T)
end

body = Body(BBTTMaterial(), pos, vol)
    material!(body, horizon=r_h, E=1.9e11, rho=7930.0, epsilon_c=0.2, kc=16.3, 
                    aph=17.3e-6, cv=500.0, rft=293.0, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick = h)  #Fe
    point_set!(p->(0.011 <= p[2] <= 0.051 && p[1] <= 0.008), body, :TC)
    material!(body, :TC; horizon=r_h, E=9e10, rho=2650, epsilon_c=0.12, kc=11,
                    aph=1.4e-6, cv=60, rft=298.0, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick = h, kp_T=f_kp_cc, cv_T=f_cv_cc)
    point_set!(p->((p[2] < 0.011 || p[2] > 0.051) && p[1] <= 0.011), body, :Cu)                       
    material!(body, :Cu; horizon=r_h, E=1.1e11, rho=8960.0, epsilon_c=0.2, kc=285.0, 
                    aph=16.5e-6, cv=385.0, rft=293.0, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick = h)                      
    point_set!(p->((0.011 < p[2] < 0.015 || 0.047 < p[2] < 0.051) && 0.008 <= p[1] <= 0.027), body, :AlO)  
    material!(body, :AlO; horizon=r_h, E=3.0e11, rho=3700.0, epsilon_c=0.18, kc=20.0, 
                    aph=7.0e-6, cv=750.0, rft=293.0, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick = h) 
    point_set!(p->((0.0285 < p[2] < 0.0335) && 0.008 <= p[1] <= 0.027), body, :GFRP)  
    material!(body, :GFRP; horizon=r_h, E=0.1e11, rho=1500.0, epsilon_c=0.12, kc=0.3, 
                    aph=7.0e-6, cv=800.0, rft=293.0, h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick = h)                        
    point_set!(p->((0.003 <=p[2] < 0.008 || 0.054 <= p[2] <= 0.059) && 0.011 <= p[1] <= 0.032), body, :cool)  

# Give the simple boundary conditions
hsource = zeros(1, size(pos,2))
for i in geo.bc_nodes
    if pos[1, i] < 0.028
        hsource[1, i] = 1.0e10
    end
end 

hsource_databc!(body, hsource, :all_points)
temperature_bc!(t -> 298.0, body, :cool)
temperature_ic!(body, :all_points, 600.0)

vv = Thermstep(steps = n_step,stepsize = Δt, safety_factor = 0.9)
job = Job(body, vv; freq=n_fq, path=joinpath(@__DIR__, "PD_results"), 
          fields=(:temperature, :hsource))
A = FSI_submit(job, "T")






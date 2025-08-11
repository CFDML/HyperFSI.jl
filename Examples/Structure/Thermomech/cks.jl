using HyperFSI

l = 0.1
w = 0.1
a = 0.05
b = 0.025
h = 0.0005
Δx = l/200
δ = 4.015Δx 
dt = 1.0e-8
nth = 9000

pos, vol =  uniform_box(l, w+3Δx, h, Δx;center=(0.0, -1.5Δx, 0.0))
body = Body(BBTMMaterial(), pos, vol)
material!(body; horizon=δ, E=1.90e11, rho=8000, epsilon_c=0.0103, kc=16.2, aph=17.6e-6, cv=477.0, rft=285.0, 
            h=25.0, hσ=5.73e-8, hϵ=0.9, tem∞=0.0, thick=h) 

point_set!(p -> (-w/2 < p[2] < -w/2+b && p[1]<= -l/2+4Δx), body, :vx)
point_set!(p -> (p[2] < -w/2), body, :fy)
point_set!((p -> p[1] < -l/2+a && -w/2+b < p[2] < -w/2+b+3Δx) , body, :ck_up)
point_set!((p -> p[1] < -l/2+a && -w/2+b-3Δx < p[2] < -w/2+b) , body, :ck_dw)
point_set!((p -> p[1] > l/2-3Δx) , body, :no_fail_zone)

velocity_bc!(t -> t < 1e-6 ? t/1e-6*16.5 : 16.5, body, :vx, 1)
velocity_bc!(t -> 0.0, body, :fy, 2)
precrack!(body, :ck_up, :ck_dw) 
no_failure!(body, :no_fail_zone)

vv = Thermomechstep(steps=nth, stepsize = dt)
job = Job(body, vv; freq=10, path="results/wk", fields=(:displacement, :velocity, :damage, :temperature))
A =  FSI_submit(job, "TM")



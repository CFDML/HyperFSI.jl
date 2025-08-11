using HyperFSI

    pos = [ -0.25   0.25  -0.25  0.25
             -0.25  -0.25   0.25  0.25
               0.0    0.0    0.0   0.0]
    vol = [0.125, 0.125, 0.125, 0.125]
    height = 0.5
    body = Peridynamics.Body(BBTMaterial(), pos, vol)
    Peridynamics.material!(body; horizon=2, E=2.0e11, rho=1.0, epsilon_c=1e-2, kc=1.0, aph=1.0e-6, cv=1.0, rft=273, 
            h=25.0, hσ=1e-8, hϵ=1.0, tem∞=0.0, thick=height) 

    Peridynamics.point_set!(p -> p[1]< 0.0, body, :radiation)
    Peridynamics.point_set!(p -> p[1]> 0.0, body, :convention)
    temperature_ic!(body, :all_points, 1000.0)

    vv = Thermstep(steps=1, stepsize = 1e-6) # Time solver setup for thermal diffusion simulation in structure
    job = Peridynamics.Job(body, vv)
    ret = FSI_submit(job, "T"; quiet=true)


function logo_init_logs(options::Peridynamics.AbstractJobOptions)
    Peridynamics.mpi_isroot() || return nothing
    Peridynamics.print_log(Hyperfsi_banner())
    Peridynamics.print_log(Peridynamics.get_run_info())
    Peridynamics.set_progress_bars!()
    options.export_allowed || return nothing
    mkpath(options.vtk)
    logo_init_logfile(options)
    return nothing
end

function logo_init_logfile(options::Peridynamics.AbstractJobOptions)
    open(options.logfile, "w+") do io
        write(io, get_logfile_head_fsi())
        write(io,Hyperfsi_banner(color=false))
        write(io, Peridynamics.get_run_info())
    end
    return nothing
end

function get_logfile_head_fsi()
    msg = "LOGFILE CREATED ON "
    msg *= Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
    msg *= "\n"
    msg *= "\n"
    return msg
end

function Hyperfsi_banner(; color::Bool=true, indentation::Int=10)
    indent = indentation > 0 ? " "^indentation : ""
    if color
        hyper_color = "\e[38;2;255;69;0m"  
        fsi_color = "\e[38;2;65;105;225m"      
        bold = "\e[1m"
        reset = "\e[0m"
    else
        hyper_color = fsi_color = bold = reset = ""
    end

    msg  = indent * "$(hyper_color)    __  __                           $(reset)$(fsi_color)______ _____  ____$(reset)\n"
    msg *= indent * "$(hyper_color)   / / / /__  __ ____   ___   _____ $(reset)$(fsi_color)/ ____// ___/ /  _/$(reset)\n"
    msg *= indent * "$(hyper_color)  / /_/ // / / // __ \\ / _ \\ / ___/$(reset)$(fsi_color)/ /_    \\__ \\  / / $(reset)\n"
    msg *= indent * "$(hyper_color) / __  // /_/ // /_/ //  __// /   $(reset)$(fsi_color)/ __/   ___/ /_/ / $(reset)\n"
    msg *= indent * "$(hyper_color)/_/ /_/ \\__, // .___/ \\___//_/   $(reset)$(fsi_color)/_/     /____//___/ $(reset)\n"
    msg *= indent * "$(hyper_color)       /____//_/                   $(reset)\n"
    
    msg *= indent * "\n"
    msg *= indent * "$(hyper_color)Hyper$(reset)$(fsi_color)FSI$(reset): Extensible framework for coupled fluid-structure-thermal multiphysics\n"
    msg *= indent * "Copyright © 2025 Shiwei Hu\n\n"
    
    return msg
end

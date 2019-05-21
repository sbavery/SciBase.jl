using Pkg

if lowercase(get(ENV, "CI", "false")) == "true"

    let basepython = get(ENV, "PYTHON", "python3")
        envpath = joinpath(@__DIR__, "env")
        run(`virtualenv --python=$basepython $envpath`)

        if Sys.iswindows()
            python = joinpath(envpath, "Scripts", "python.exe")
        else
            python = joinpath(envpath, "bin", "python3")
        end

        run(`$python -m pip install numpy`)
        run(`$python -m pip install scipy`)
        run(`$python -m pip install matplotlib`)

        ENV["PYTHON"] = python
        Pkg.build("PyCall")
    end
end
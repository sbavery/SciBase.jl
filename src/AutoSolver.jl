module AutoSolver

#files = [
#		]

#for f in files
#    include(f * ".jl")
#end

function modName(mod::Module)
    if mod == Main
        mod_name = string(mod)
    else
        mod_name = string(mod)[6:end]
    end

    return mod_name
end

# check if a symbol is a nested module
issubmodule(m::Module, s::Symbol) = isa((@eval $(m).$(s)),Module) && !occursin(string(s),modName(m))
    
# get a list of all submodules of module m
function submodules(m::Module, mod_arr::Array{Module})
    modules = filter(x->issubmodule(m,x),names(m,all=true))
    push!(mod_arr,m)
    for mod in modules
        submodules((@eval $(m).$(mod)), mod_arr)
    end
end

function findArgs(func::Function)
    methods_str = methods(func)
    arg_ind_1 = findlast("(",string(methods_str))[1]+1
    arg_ind_2 = findlast(")",string(methods_str))[1]-1

    arg_s = string(methods_str)[arg_ind_1:arg_ind_2]
    args = []
    str = ""
    isArg = true

    for char in arg_s
        str = str
        if isArg == false && (char == ',' || char == ';')
            isArg = true
        elseif isArg == true
            if char == ',' || char == ';'
                push!(args,Symbol(strip(str)))
                str = ""
            elseif char == '=' || char == ':'
                isArg = false
                push!(args,Symbol(strip(str)))
                str = ""
            else
                str = string(str,char)
            end
        end
    end
    if str != ""
        push!(args,Symbol(strip(str)))
    end

    return args
end

function moduleDict(mod::Module, funcs::Dict)
    mod_name = modName(mod)
    
    for func in names(mod, all=true)
        if (func != :eval && func != :include && 
            func != Symbol(mod_name) && func != :Base &&
            func != :Core)
            if !occursin("#",string(func))
                funcType = string(typeof(@eval $(mod).$(func)))
                if occursin("typeof",funcType)
                    func_args = findArgs(@eval $(mod).$(func))
                    funcs[func] = func_args
                end
            end
        end
    end
    return funcs
end

function evalVariableInputs(mod::Module, func::Symbol, args::Array{Any})
    global output = 0
    args = tuple(args...)
    output = (@eval $(mod).$(func)($(args)...))
    
    return output
end

function varExplode(var)
    var_arr = []
    var_name = string(var)
    var_string = ""

    for c in var_name
        var_string = string(var_string,c)
        push!(var_arr,Symbol(var_string))
    end
    
    return var_arr
end

function varExplodeDict(var_dict::Dict)
    inp_dict = Dict()
    out_dict = Dict()
    for pair in var_dict
        out_dict[pair[1]] = varExplode(pair[1])
        for var in pair[2]
            inp_dict[var] = varExplode(var)
        end
    end
        
    return out_dict, inp_dict
end

function varMapArray(d::Dict,var::Symbol)
    map_array = []
    
    for pair in d
        for var_n in pair[2]
            if var_n == var
                push!(map_array, pair[1])
            end
        end
    end
    return map_array
end

function crossVarMap(d1::Dict, d2::Dict)
    d12 = Dict()
    d21 = Dict()
    
    for pair in d1
        map_array = []
        for var in pair[2]
            map_array = vcat(map_array, varMapArray(d2,var))
        end
        d12[pair[1]] = vcat(pair[1], map_array)
    end
    
    for pair in d2
        map_array = []
        for var in pair[2]
            map_array = vcat(map_array, varMapArray(d1,var))
        end
        d21[pair[1]] = vcat(pair[1], map_array)
    end
    return d12,d21
end

function argCheck(mod::Module, arg::Symbol, arg_arr::Array{Any}, inp_arr::Array{Any}, inp_map::Dict)
    valid_inputs = 0
    mods = [Main,mod]
    
    for arg_inp in inp_map[arg]
        for m in mods
            if isdefined(m,arg_inp) && !isa((@eval $(m).$(arg_inp)),Function)
                if @eval $(m).$(arg_inp) != nothing
                    push!(inp_arr, @eval $(m).$(arg_inp))
                    push!(arg_arr, arg_inp)
                    valid_inputs = 1
                    break
                end
            else
                valid_inputs = 0
            end
        end
        if valid_inputs == 1; break; end
    end
    
    return valid_inputs
end

function outputDict(funcs::Dict, outputs::Dict; mod::Module=Main, inp_map::Dict=Dict())
    valid_inputs = 0
    global inp = []
    global arg_arr = Array{Symbol}
    all_solved = false
    none_solved = false
    num_solved = 0
    global output = 0
    
    while all_solved == false && none_solved == false
        for func in funcs
            if !haskey(outputs,func[1]) && !occursin(string(func[1]),"outputDict")
                valid_inputs = 0
                inp = []
                arg_arr = []
                
                for arg in func[2]
                    valid_inputs = argCheck(AutoSolver, arg, arg_arr, inp, inp_map)
                    
                    if arg == func[2][end] && valid_inputs == 1
                        try
                            output = evalVariableInputs(mod, func[1], inp)
                            println("Func: $(func[1])$(tuple(func[2]...))")
                            println("Inputs:")
                            for i = 1:length(arg_arr)
                                println("$(arg_arr[i]) = $(inp[i])")
                            end
#                             println("Input Args: $(tuple(arg_arr...))")
#                             println("Inputs: $(tuple(inp...))")
                            println("Output: $(output)\n")
                            num_solved+=1
                            outputs[func[1]] = output
                            out = string(func[1])
                            @eval $(Symbol(out))=output
                        catch
                            println("Func: $(func[1])$(tuple(func[2]...))")
                            println("Inputs:")
                            for i = 1:length(arg_arr)
                                println("$(arg_arr[i]) = $(inp[i])")
                            end
                            println("Output: Error\n")
                        end
                    elseif valid_inputs == 0
                        break
                    end
                end
            end
        end
        if num_solved == length(funcs)
            all_solved = true
        elseif num_solved == 0
            none_solved = true
        end
        num_solved = 0
    end
end

function initArgCheck(funcs::Dict, init_args::Array{Any})
    for func in funcs
        for arg in func[2]
            if isdefined(Main,arg) && @eval Main.$(arg) != nothing
                push!(init_args,arg)
            end
        end
    end
end

function resetArgs(funcs::Dict, init_args::Array{Any})
    for func in funcs
        for arg in func[2]
            if !any(x->x==arg, init_args)
                @eval $(arg)=nothing
            end
        end
    end
end

function autoSolveOutputs(mods::Array{Module})
    for mod in mods
        in_mod = mod
        @show mod
        funcs = Dict()
        outputs = Dict()
        moduleDict(in_mod,funcs)
        if !isempty(funcs)
            println("::::MODULE $(mod)::::")
            out_d, in_d = varExplodeDict(funcs)
            inp_map, out_map = crossVarMap(in_d,out_d)
            println()
            outputDict(funcs, outputs, mod=in_mod, inp_map=inp_map)
            @show outputs
            println()
        end
    end
end

end # AutoSolver
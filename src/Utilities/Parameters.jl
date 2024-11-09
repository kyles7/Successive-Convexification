module Parameters

export load_parameters

using YAML

"""
    load_parameters(file_path::String) -> Dict

Loads simulation and system parameters from a YAML configuration file into a Dict.

# Arguments
- `file_path::String`: Path to the YAML configuration file.

# Returns
- `params::Dict`: Dictionary containing all parameters.
"""
function load_parameters(file_path::String)
    # Load the YAML file
    params = YAML.load_file(file_path)
    
    # Convert nested structures (like x0 and u_guess) to arrays
    params["x0"] = collect(params["x0"])
    params["u_guess"] = collect(params["u_guess"])
    
    # If x_target is provided, convert it to an array
    if haskey(params, "x_target")
        params["x_target"] = collect(params["x_target"])
    end
    
    return params
end

end # module Parameters

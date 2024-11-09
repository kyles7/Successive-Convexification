# Parameters.jl

module Parameters

using YAML

function load_parameters(file_path::String)
    params = YAML.load_file(file_path)
    return params
end

end # module Parameters

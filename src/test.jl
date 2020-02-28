include("prm.jl")

using Random

# TODO: Get rid of this at some point. Good for initial testing though
# for repeatability
Random.seed!(1234)

const RotationMatrix = Array{Float64,2}

function generate_obstacles(n::Int64)::Array{Number,2}
    out = reshape(collect(Number, 1:5*n),(n,5))

    for i = 1:n
        out[i,:] = [(rand(Float64,(1,4))*10) rand((0,1))::Int64][:]
    end

    return out
end

function rotation_matrix(theta_x, theta_y, theta_z)::RotationMatrix
    cx = cos(theta_x)
    sx = sin(theta_x)
    cy = cos(theta_y)
    sy = sin(theta_y)
    cz = cos(theta_z)
    sz = sin(theta_z)

    rx = [1 0 0; 0 cx -sy; 0 sx cx]
    ry = [cy 0 sy; 0 1 0; -sy 0 cy]
    rz = [cz -sz 0; sz cz 0; 0 0 1]
    return rz * ry * rx
end

function generate_states()::Tuple{HomogeneousMatrix, HomogeneousMatrix}
    Rs = rotation_matrix(rand(), rand(), rand())
    ds = [-2.0; -2.0; 4.0]
    Rg = rotation_matrix(rand(), rand(), rand())
    dg = [12.0; 12.0; 12.0]
    s = [Rs ds; 0 0 0 1]
    g = [Rg dg; 0 0 0 1]
    return (s, g)
end

O = generate_obstacles(5)

pprint_matrix(O)

PRM(generate_states()..., O)

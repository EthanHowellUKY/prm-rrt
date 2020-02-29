using LinearAlgebra

macro c_assert(boolean)
    message = string("Assertion: ", boolean, " failed")
    :($(esc(boolean)) || error($message))
end

struct Cylinder
    p::Array{Float64,1}
    r::Float64
end

struct Sphere
    p::Array{Float64,1}
    r::Float64
end

mutable struct World
    cylinders::Array{Cylinder,1}
    spheres::Array{Sphere,1}
end

const HomogeneousMatrix = Array{Float64, 2}

# World constructor
World() = World([],[])

function add_object(w::World, data::Array{Number,1})
    # We expect a 1x5 row vector
    @c_assert size(data) == (5,)

    # Convert our values to float64s if necessary
    # This just makes the API a bit more friendly with
    # the abstract Number type as an input
    p = convert(Array{Float64}, data[1:3])
    r = convert(Float64, data[4])

    if data[5] == 0 # Sphere
        s = Sphere(p,r)
        push!(w.spheres, s)
    else
        c = Cylinder(p,r)
        push!(w.cylinders, c)
    end
end

function check_state_collision(w::World, hm::HomogeneousMatrix)::Bool
    s = hmatrix_to_sphere(hm)
    check_state_collision(w, s)
end

function check_state_collision(w::World, s::Sphere)::Bool
    for c in w.cylinders
        if check_collision(s, c)
            return true
        end
    end

    for ws in w.spheres
        if check_collision(s, ws)
            return true
        end
    end

    return false
end

function check_path_collision(w::World,
                              a::Array{Float64,1},
                              b::Array{Float64,1})::Bool
    # Parameterize our line as a + tx, where x = b - a s.t.
    # when t = 1, a + tx = b
    x = b - a
    # Take L2 norm of x. This is distance of the line between
    # the two points, and should tell us roughly how many points
    # to sample along the line. Since our radius is 1, we'll check
    # the points at each radius length, or 2*r*round(||x||)
    frac = ceil(l2_norm(x)) * 2

    i = 0.0
    while i <= 1.0
        current_point = a + x * i
        s = Sphere(current_point,1)
        if check_state_collision(w,s)
            return true
        end
        i += 1 / frac
    end

    return false
end

function check_collision(a::Sphere, b::Sphere)::Bool
    d = a.p - b.p
    dp = dot(d, d)
    min_dist = (a.r + b.r) ^ 2
    if dp < min_dist
        return true
    else
        return false
    end
end

function check_collision(a::Sphere, b::Cylinder)::Bool
    # First, check if we are too high to collide
    if a.p[3] - a.r > b.p[3]
        return false
    elseif a.p[3] < b.p[3]
        # Next, if our center is still within the height
        # of the cylinder, do circle collision check
        s_xy = a.p[1:2]
        c_xy = b.p[1:2]
        d = s_xy - c_xy
        dp = dot(d, d)
        min_dist = (a.r + b.r) ^ 2
        if dp < min_dist
            return true
        else
            return false
        end
    else
        # Finally, the hardest case is when the ball
        # would hit the top edge of the cylinder
        # First, we need to find the radius of the
        # sphere at the height of the edge, then do a
        # simple circle test again
        zd = a.p[3] - b.p[3]
        @c_assert zd > 0 && zd < a.r
        r = sqrt(a.r^2 - zd^2)
        s_xy = a.p[1:2]
        c_xy = b.p[1:2]
        d = s_xy - c_xy
        dp = dot(d, d)
        min_dist = (r + b.r)^2
        if dp < min_dist
            return true
        else
            return false
        end
    end
end

function hmatrix_get_d(m::HomogeneousMatrix)::Array{Float64,1}
    return m[1:3,4]
end

function hmatrix_to_sphere(m::HomogeneousMatrix)::Sphere
    return Sphere(hmatrix_get_d(m),1)
end

function l2_norm(v::Array{Float64,1})
    sqrt(dot(v,v))
end

function pprint_matrix(m::Array{T,2}) where {T}
    rows, _ = size(m)
    for i=1:rows
        println(m[i,:])
    end
end

using LinearAlgebra
using Distributions

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

function check_state_collision(w::World, pt::Array{Float64,1})::Bool
    s = point_to_sphere(pt)
    check_state_collision(w, s)
end

function check_state_collision(pt1::Array{Float64,1}, pt2::Array{Float64,1})::Bool
    s1 = point_to_sphere(pt1)
    s2 = point_to_sphere(pt2)
    check_collision(s1, s2)
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

function check_state_collision(pt1::Array{Float64,1}, pt2::Array{Float64,1}, O::Array{Float64, 2})::Bool
    numO = length(O[:,1]);

    s1 = point_to_sphere(pt1);
    s2 = point_to_sphere(pt2);

    for ii = 1:numO
        if O[ii,5] == 0
            obj = Sphere(O[ii,1:3], O[ii,4]);
        else
            obj = Cylinder(O[ii,1:3], O[ii,4]);
        end

        s1_coll_obj = check_collision(s1, obj);
        s2_coll_obj = check_collision(s2, obj);
        s1_coll_s2 = check_collision(s1,s2);

        if s1_coll_obj || s2_coll_obj || s1_coll_s2
            retval = true;
        else
            retval = false;
            break
        end
    end

    return retval
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

function hmatrix_to_sphere(m::HomogeneousMatrix)::Sphere
    return Sphere(m[1:3,4],1)
end

function point_to_sphere(pt::Array{Float64,1})::Sphere
    return Sphere(pt,1)
end

function pprint_matrix(m::Array{T,2}) where {T}
    rows, _ = size(m)
    for i=1:rows
        println(m[i,:])
    end
end

function dist(n::Array{Float64,1}}, p::Array{Float64,1})::Float64
    return sqrt( (n[1] - p[1])^2 + (n[2] - p[2])^2 + (n[3] - p[3])^2 )
end

function steer(p_rand::Array{Float64,1}, p_n::Array{Float64,1}, val::Float64, eps::Int64)::Array{Float64,1}
    p_new::Array{Float64,1}

    if val >= eps
        for ii = 1:length(p_rand)
            push!(p_new, p_n[ii] + (p_rand[ii]-p_n[ii])/dist(p_rand, p_n));
        end
    else
        for ii = 1:length(p_rand)
            push!(p_new, p_rand[ii]);
        end
    end

    return p_new
end

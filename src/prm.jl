# Attribution: http://www.cs.columbia.edu/~allen/F15/NOTES/Probabilisticpath.pdf
# Notes used for PRM

include("utils.jl")

mutable struct Graph
    V::Array{Array{Float64,1},1}
    E::Array{Tuple{Int64,Int64},1}
end

Graph() = Graph([],[])

function generate_graph(n::Int64,
                        k::Int64,
                        s::HomogeneousMatrix,
                        g::HomogeneousMatrix,
                        w::World)::Graph
    # TODO: Look into using column major ordering for efficiency
    # Preallocate memory for efficiency since we know how many
    # nodes we're going to need
    V = zeros(n,3)

    # Generate n valid points through random sampling
    # TODO: What bounds do we have? Is z = 0.0 a boundary? Any
    # guarantees on x,y boundaries as well?
    # Maybe define limits based on start and goal state range?
    i = 1
    while i <= n
        p = [(rand() - 0.5) (rand() - 0.5) rand()] .* 40
        s = Sphere(p, 1.0)
        if !check_state_collision(w, s)
            # If we find a non-colliding point, add it
            V[i,:] = p
        end
        i+=1
    end

    # Once we have n points, we need to connect the k nearest
    # neighboring points
    for j=1:n
        closest_neighbors = zeros(k)
        idxs = Int64[]
        for ci = eachindex(closest_neighbors)
            closest_neighbors[ci] = Inf
            push!(idxs, -1)
        end
    end
end

function PRM(s, g, O)
    # First, build the world
    world = World()

    # Next, add our objects to the world
    rows, _ = size(O)
    for i = 1:rows
        add_object(world, O[i,:])
    end

    println("Start state:")
    pprint_matrix(s)
    println("\nGoal state:")
    pprint_matrix(g)

    if check_state_collision(world, s)
        error("Invalid starting state. Collision detected.")
    elseif check_state_collision(world, g)
        error("Invalid goal state. Collision detected.")
    end
end

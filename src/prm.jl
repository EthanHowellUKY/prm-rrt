# Attribution: http://www.cs.columbia.edu/~allen/F15/NOTES/Probabilisticpath.pdf
# Notes used for PRM

include("utils.jl")

using LinearAlgebra

struct Graph
    V::Array{Float64,2}
    E::Array{Int64,2}
end

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
            i += 1
        end
    end

    # Once we find our other points, add our start and goal states
    # so that paths can be generated for them
    insert!(V, 1, hmatrix_get_d(s))
    push!(V, hmatrix_get_d(g))

    # Once we have n points, we need to connect the k nearest
    # neighboring points. These edges will be E
    E = reshape(collect(1:(n+2)*10),(n+2,10))
    closest_neighbors = zeros(k)

    for j=1:n+2
        for ci = eachindex(closest_neighbors)
            closest_neighbors[ci] = Inf
            E[j,ci] = -1 # Sentinel value
        end

        current_point = V[j,:]

        # Now look through all points where j != h
        # If its closer than what we've seen so far, add it
        for h=1:n+2
            if h == j
                continue
            end

            test_point = V[h,:]
            diff = test_point - current_point
            dist = dot(diff, diff)
            for (idx, val) in enumerate(closest_neighbors)
                if dist < val
                    # If this is a close point, then do
                    # the collision check since it is expensive
                    if !check_path_collision(w,
                                             current_point,
                                             test_point)
                        insert!(closest_neighbors, idx, dist)
                        E[j,idx] = h
                        # We want to keep our list size k, so pop
                        # off the last value in the list since it's
                        # not within the k closest neighbors
                        pop!(closest_neighbors)
                    end

                    break
                end
            end
        end
    end

    # Return the generated graph with points
    Graph(V,E)
end

# Attribution: https://brilliant.org/wiki/dijkstras-short-path-finder/
function shortest_path(g::Graph)::Array{Int64,1}
    n = size(g.V)[1]
   
    weights = zeros(n+2)
    for i=2:n+2
        weights[i] = Inf
    end

    Q = collect(2:n+2)

    cur_node = 1
    cur_dist = weights[cur_node]

    while size(Q)[1] > 0
        # Do work
        cur_edges = g.E[cur_node,:]
        for edge in cur_edges
            # Check for sentinel values
            if edge < 1
                break
            end

            p1 = g.V[cur_node,:]
            p2 = g.V[edge,:]
            diff = p1 - p2
            dist = dot(diff,diff)
            new_dist = cur_dist + dist
            if new_dist < weights[edge]
                weights[edge] = new_dist
            end
        end

        # ...then find min index, pop it and continue again until
        # Q is empty
        cur_node = min_dist_index(Q,dist)
        cur_dist = weights[cur_node]
        remove!(Q,cur_node)
    end

    return weights
end

function min_dist_index(Q::Array{Int64,1}, dist::Array{Float64,1})
    min_i = 0
    min_val = Inf
    for i in Q
        if dist[i] < min_val
            min_val = dist[i]
            min_i = i
        end
    end

    return i
end

function find_path(g::Graph, d::Array{Float64,1})
    start_idx = 1 # Start from start node in graph
    goal_idx = size(g.V)[1]
    find_path(g,d,goal_idx,start_idx)
end

function find_path(g::Graph, d::Array{Float64,1}, idx, end_idx)
    if idx == end_idx
        return (g.V[goal_idx,:])
    end

    edges = g.E[idx,:]
    min_dist = d[idx]
    idxs = []
    for edge in edges
        if d[edge] < min_dist
            push!(idxs, edge)
        end
    end

    # Sort all of the idxs that are shorter
    sort!(idxs)

    for idx in idxs
        val = find_path(g,d,idx,end_idx)
        # If we find a valid decreasing distance match from
        # the end, return it concatenated with the total tuple
        # of points
        if val
            return (val..., g.V[idx,:])
        end
    end

    # If none of our paths matched, return false
    return false
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

    # Generate our valid vertex/edge graph based on collision from
    # obstacles in the world.
    graph = generate_graph(1000, 10, s, g, world)
    # Find the shortest path to the start node for graph
    dist = shortest_path(graph)
    # Finally, take shortest graph representation and plot path
    # between the final node and the start node
    path = find_path(graph,dist)

    println("Final path vertices")
    len = size(path,1)
    for i=1:len
        println("Point ", i, ": ", path[i])
    end

    # TODO: Need to convert 3d points to 4x4 homogeneous matrices
    # At this point, path is a tuple of 3d translations moving between
    # the start and end state. I believe these are currently in the main
    # frame defined by the world origin, or (0,0,0) with no rotation.
    # That means we'll need to convert these to translations from the
    # point of view of the robot, namely through setting up an initial
    # homogeneous matrix and then running through identity rotations
    # until the goal state is acheived. At that point, we'll need to
    # do one last rotation to make sure we acheive both the translation
    # and rotation of the goal state as required by the assignment
end

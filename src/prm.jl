# Attribution: http://www.cs.columbia.edu/~allen/F15/NOTES/Probabilisticpath.pdf
# Notes used for PRM

include("utils.jl")

mutable struct Graph
    V::Array{Array{Float64,1},1}
    E::Array{Tuple{Int64,Int64},1}
end

mutable struct Point
    coords::Array{Float64,1}
    cost::Float64
    parent::Float64
end

# Leave as immutable so that the final path cannot be altered
struct Path
    coords::Array{Float64,1}
    cost::Float64
    parent::Float64
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

function RRT(s::HomogeneousMatrix, g::HomogeneousMatrix, O::Array{Float64,2})::Tuple{Array{Float64,1},2}

    # s -> 4x4 HomogeneousMatrix representing initial state
    # g -> 4x4 HomogeneousMatrix representing goal state
    # O -> Nx5 Matrix representing objects in the space

    x_max = 0;
    y_max = 0;
    z_max = 0;
    eps = 20;
    pt_fnd = false;
    numNodes = 4000;

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

    # Sanity check to validate the initial state and goal state
    if check_state_collision(world, s)
        error("Invalid starting state. Collision detected.")
    elseif check_state_collision(world, g)
        error("Invalid goal state. Collision detected.")
    end

    init_coords = s[1:3,4]; # Coordinates of start state
    goal_coords = g[1:3,4]; # Coordinates of goal state

    init_state = Path(init_coords, 0, 0);
    nodes = [init_state];

    x_max = maximum(Graph.V[:,1]); ################# Note sure if this is where the domain will be defined
    y_max = maximum(Graph.V[:,2]); ################# Used as placeholder for time being
    z_max = maximum(Graph.V[:,3]); #################

    for ii = 1:numNodes

        # Generate a random point somwhere within the bounds of the environment
        p_rand = [floor(Int64,rand(Uniform(0,1),1)*x_max), floor(Int64,rand(Uniform(0,1),1)*y_max), floor(Int64,rand(Uniform(0,1),1)*z_max)];

        # Set exit flag to false then check to see if we are close to the goal state
        pt_fnd = false;
        for jj = 1:length(nodes)
            if norm(nodes[jj].coords - goal_coords) <= 25*sqrt(2)
                pt_fnd = true;
                break
            end
        end

        # Break out of the loop if the goal state has been reached
        if pt_fnd
            break
        end

        ndist = [];
        for jj = 1:length(nodes)
            n = nodes[jj];                                              # x-y-z coords of current node
            tmp_dist = dist(n.coords, p_rand);                          # Distance between current node and randomly sampled point
            push!(ndist, tmp_dist);
        end

        (val, idx) = findmin(ndist);                                    # Get the value and index associated with closest node
        p_near = nodes[idx];                                            # Get the coordinates of the node closest to the randomly sampled point

        p_new = Point(steer(p_rand, p_near.coords, val, eps), 0, 0);    # Move in the direction of the nearest point
        if !check_state_collision(p_near.coords, p_new.coords)
            p_new.cost = dist(p_new.coords, p_near.coords) + p_near.cost;

            p_nearest = [Point([0,0,0],0,0)];
            r = 60;

            for jj = 1:length(nodes)
                if !check_state_collision(nodes[jj].coords, p_new.coords) && dist(nodes[jj].coords, p_new.coords) <= r
                    push!(p_nearest, nodes[jj]);
                end
            end

            p_min = p_near;
            c_min = p_new.cost;

            for jj = 1:length(p_nearest)
                p_curr = p_nearest[jj];
                c_curr = p_curr.cost + dist(p_curr.coords, p_new.coords)
                if !check_state_collision(p_curr.coords, p_new.coords) && c_curr < c_min
                    p_min = p_curr;
                    c_min = c_curr;
                end
            end

            for jj = 1:length(nodes)
                if nodes[jj].coords == p_min.coords
                    p_new.parent = jj;
                end
            end
            new_node = Path(p_new.coords, p_new.cost, p_new.parent)
            push!(nodes, new_node);
        end
    end

    D = [];
    for jj = 1:length(nodes)
        curr_dist = dist(nodes[jj].coords, goal_coords);
        push!(D, curr_dist);
    end

    (val, idx) = findmin(D);
    p_final = nodes(idx);
    goal_state = Point(goal_coords, 0, idx);
    push!(nodes, goal_state);

    final_path = zeros(length(nodes),3);
    for jj = 1:length(nodes)
        for kk = 1:3
            final_path[jj,kk] = nodes[jj].coords[kk];
        end
    end

    return tuple(final_path)
end

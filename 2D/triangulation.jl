using DelaunayTriangulation
using Plots
pythonplot()

Lpts = [(0.0,0.0), (0.0,2.0), (1.0,2.0), (1.0,1.0), (2.0,1.0), (2.0,0.0)]

boundary_nodes = [6, 5, 4, 3, 2, 1, 6]

tri = triangulate(Lpts; boundary_nodes)

A = get_area(tri)
refine!(tri; min_angle=30.0, max_area = 0.1*A)

# Needed for FEM later
nodes = get_points(tri)
elements = collect(each_solid_triangle(tri))

p = plot(aspect_ratio=:equal, legend=false)

for (i,j,k) in elements
    xs = [nodes[i][1], nodes[j][1], nodes[k][1], nodes[i][1]]
    ys = [nodes[i][2], nodes[j][2], nodes[k][2], nodes[i][2]]

    plot!(p, xs, ys, lw=1)
end

display(p)
readline()


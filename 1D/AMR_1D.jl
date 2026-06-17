using Plots
pythonplot()

using LinearSolve

function stiff_mat(nodes, a, kappa)
	n = length(nodes)-1
	A = zeros(n+1, n+1)
	for i in 1:n
		h = nodes[i+1] - nodes[i]
		xmid = (nodes[i+1] + nodes[i])/2
		amid = a(xmid)
		A[i,i] = A[i,i] + amid/h
		A[i, i+1] = A[i, i+1] - amid/h
		A[i+1, i] = A[i+1, i] - amid/h
		A[i+1, i+1] = A[i+1, i+1] + amid/h
	end
	A[1,1] = A[1,1] + kappa[1]
	A[n+1, n+1] = A[n+1, n+1] + kappa[2]
	return A
end

function load_vec(x, f, kappa, g)
	n = length(x)-1
	b = zeros(n+1)
	for i in 1:n
		h = x[i+1] - x[i]
		b[i] = b[i] + f(x[i])*(h/2)
		b[i+1] = b[i+1] + f(x[i + 1])*(h/2)
	end
	b[1] = b[1] + kappa[1]*g[1]
	b[end] = b[end] + kappa[2]*g[2]
	return b
end

# Solve -∂^2u/∂x^2 = δ(x)
delta(x) = exp(-100*abs(x-0.5)^2)

function mesh_refine(x, f)
	n = length(x)-1
	residuals = zeros(n)
	for i in 1:n
		h = x[i + 1] - x[i]
		a = f(x[i])
		b = f(x[i+1])
		t = (a^2 + b^2) * (h/2)
		residuals[i] = h*sqrt(t)
	end
	ref_par = 0.9
	for i in 1:length(residuals)
		if residuals[i] > ref_par*maximum(residuals)
			append!(x, (x[i+1] + x[i])/2)
		end
	end
	x = sort!(x)
	return x
end

function poisson_solver()
	h = 0.1
	x = collect(0:h:1)
	for k in 1:25
		x = mesh_refine(x, delta)
	end
	a(x) = 1
	kappa = [1e12, 1e12]
	g = [0,0]
	A = stiff_mat(x,a,kappa)
	b = load_vec(x, delta, kappa, g)
	weights_prob = LinearProblem(A,b)
	weights_sol = solve(weights_prob)
	weights = weights_sol.u
	return x,  weights
end
h = 0.1
x_test, u_test = poisson_solver()

println("Number of points from refinement: ", length(x_test))
p = plot(x_test, u_test, seriestype = :path, marker = :circle, markercolor = :red, label="numerical FEM", xlabel="x", ylabel="u(x)", title="solution to 1-D Poisson equation, 25 mesh refinements")
display(p)
readline()



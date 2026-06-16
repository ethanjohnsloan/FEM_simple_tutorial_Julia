using Plots
pythonplot()

using LinearSolve

using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

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


# Example from 2.4.2, pp.37 of Larson & Bengzon
conductivity(x) = 0.1*(5 - 0.6*x)
source(x) = 0.03*(x-6)^4

function poisson_solver()
	h = 0.1 #mesh size
	x = 2:h:8 #initializes mesh
	kappa = [1e6, 0]
	g = [-1, 0]
	A = stiff_mat(x, conductivity, kappa)
	b = load_vec(x, source, kappa, g)
	weights_prob = LinearProblem(A,b)
	weights_sol = solve(weights_prob)
	weights = weights_sol.u
	return weights
end
h = 0.1
x_test = 2:h:8
u_test = poisson_solver()


p = plot(x_test, u_test, label="numerical FEM", xlabel="x", ylabel="T", title="temperature profile of a rod with source")
display(p)
readline()

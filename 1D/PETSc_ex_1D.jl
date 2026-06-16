# Same as FEM_1D.jl, but using PETSc to solve
# PETSc makes this seem a bit more complicated, but
# it is better suited for MPI implementation

using PETSc, MPI
using SparseArrays
using Plots
pythonplot()

PetscScalar = Float64
petsclib = PETSc.getlib(; PetscScalar = PetscScalar)
PETSc.initialize(petsclib)

# Different initializations, suggested by devs (laplacian.jl)
# Initialize PETSc with logging active and printed to REPL. This is useful to detect memory leaks
#PETSc.initialize(petsclib; log_view = true)    

# Initialize PETSc with logging active and written to "logfile.txt"
#PETSc.initialize(petsclib; log_view = true,  options=[":logfile.txt"])    # write log to a file

conductivity(x) = 0.1*(5 - 0.6*x)
source(x) = 0.03*(x-6)^4

# Set the total number of grid points
Nq = 61

# number of interior points
n = Nq - 2

# create the interior grid and get the grid spacing
x = range(PetscScalar(2), length = Nq, stop = 8)[2:(end - 1)]
Δx = PetscScalar(x.step)

kappa = [1e6, 0]
g = [-1, 0]

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


A =  sparse(stiff_mat(x, conductivity, kappa))

b = load_vec(x, source, kappa, g)

# Set up the PETSc solver
ksp = PETSc.KSP(petsclib, MPI.COMM_SELF, A; ksp_monitor = true);
b_petsc = PETSc.VecSeq(petsclib, b)

# Solve!
v = ksp \ b_petsc;

# Converts PETSc vectors into Julia arrays
x_julia = x[:]
v_julia = v[:]


# Plot
p = plot(x_julia, v_julia, label="numerical FEM", xlabel="x", ylabel="T", title="temperature profile of a rod with source")
display(p)
readline()
PETSc.destroy(ksp)
PETSc.finalize(petsclib)

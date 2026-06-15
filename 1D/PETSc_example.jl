using PETSc, Plots, SparseArrays
pythonplot()

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib, log_view=false)
PetscScalar = petsclib.PetscScalar
PetscInt = petsclib.PetscInt

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

# Same as FEM_1D.jl, but solved using PETSc... a bit overkill, just for demonstration
# Example from 2.4.2, pp.37 of Larson & Bengzon
conductivity(x) = 0.1*(5 - 0.6*x)
source(x) = 0.03*(x-6)^4

function poisson_solver_PETSc()
	h = 0.1 #mesh size
	x = 2:h:8 #initializes mesh
	kappa = [1e6, 0]
	g = [-1, 0]
	A_dense = stiff_mat(x, conductivity, kappa)
	b_dense = load_vec(x, source, kappa, g)
	n = length(b_dense)
	A = PETSc.MatSeqAIJ(A_dense) # PETSc sparse matrix format for single MPI process
	b = PETSc.VecSeq(zeros(n)) #sequential vector
	u = PETSc.VecSeq(zeros(n))
	ksp = PETSc.KSP() #kyrlov subspace solver
	PETSc.KSPSetOperators(ksp, A, A) #tells PETSc which matrix represents Au=b
	PETSc.KSPSetType(ksp, PETSc.KSPCG) #chooses CG, since A is symmetric and positive definite
	pc = PETSc.KSPGetPC(ksp) #gets preconditioner to improve convergence
	PETSc.PCSetType(pc, PETSc.PCJACOBI) #Choose Jacobi preconditioner
	PETSc.KSPSetTolerances(ksp, 1e-10, PETSc.PETSC_DEFAULT, PETSc.PETSC_DEFAULT, 1000) #sets tolerances for convergence
	PETSc.KSPSetUp(ksp) #builds internal data structures
	PETSc.KSPSolve(ksp, b, u) #solves
	u_sol =  Array(u)
	return x, u_sol
end

x, u = poisson_solver_PETSc()

p = plot(x_test, u_test, label="numerical FEM", xlabel="x", ylabel="T", title="temperature profile of a rod with source")
display(p)
readline()


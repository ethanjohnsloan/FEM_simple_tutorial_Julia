using Plots
pythonplot()

using LinearSolve

#assembles the mass matrix M
function mass_mat(x) #x is vector of partitions
	n = length(x) - 1
	M = zeros(n+1, n+1)
	for i in 1:n
		h = x[i+1] - x[i]
		M[i,i] = M[i,i] + h/3
		M[i,i+1] = M[i,i+1] + h/6
		M[i+1,i] = M[i+1,i] + h/6
		M[i+1,i+1] = M[i+1,i+1] + h/3
	end
	return M
end


#assembles the load vector b
function load_vec(x, f) #takes partition and function
	n = length(x)-1
	b = zeros(n+1)
	for i in 1:n
		h = x[i+1] - x[i]
		b[i] = b[i] + f(x[i])*(h/2)
		b[i+1] = b[i+1] + f(x[i + 1])*(h/2)
	end
	return b
end

# assembles weights for projection
function weights(A, b)
	prob = LinearProblem(A, b)
	sol = solve(prob)
	xi = sol.u
	return xi
end

function hat(i, x, nodes)
	n = length(nodes)
	if i>1 && x>=nodes[i-1] && x <= nodes[i]
		return (x - nodes[i-1])/(nodes[i] - nodes[i-1])
	elseif i<n && x>= nodes[i] && x<= nodes[i+1]
		return (nodes[i+1] - x)/(nodes[i+1] - nodes[i])
	else
		return 0
	end
end

# uses weights and hat basis functions to construct projection
function projection(x, nodes, weights)
	Phf = 0
	for i in eachindex(nodes)
		Phf = Phf + weights[i]*hat(i, x, nodes)
	end
	return Phf
end


# Example parameters and test, proof it works!
numnodes = 20
nodes = range(0,10,numnodes)
f_test(x) = sin(x)*0.05*x^2

x_cont_test = 0:0.01:10

M_test = mass_mat(nodes)
b_test = load_vec(nodes, f_test)

xi_test = weights(M_test, b_test)

proj_test = projection.(x_cont_test, Ref(nodes), Ref(xi_test))

p = plot(x_cont_test, [f_test.(x_cont_test) proj_test], labels=["actual function" "projection"])
diffplot = plot(x_cont_test, f_test.(x_cont_test) .-  proj_test)
display(p)
println("Press [Enter] to exit...")
readline()



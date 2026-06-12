using Plots
pythonplot()

function stiff_mat(nodes, a, kappa)
	n = length(nodes)-1
	A = zeros(n+1, n+1)
	for i in 1:n
		h = x[i+1] - x[i]
		xmid = (x[i+1] + x[i])/2
		amid = a[xmid]
		A[i,i] = A[i,i] + amid/h
		A[i, i+1] = A[i, i+1) - amid/h
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
end


#p = plot(x, hat.(x, 4, Ref(nodes)))
#display(p)
#readline()

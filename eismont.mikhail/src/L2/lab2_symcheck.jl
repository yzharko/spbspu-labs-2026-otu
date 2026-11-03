using Pkg
Pkg.add("Giac")

using Giac
using Giac.Commands: simplify

@giac_var s

W1 = 1 / (s + 1)
W2 = 2
W3 = 5 / (s^2 + 3 * s + 2)

W = (1 + W2 * W3 + W3) / (1 + W1 * W2 * W3 + W2 * W3 + 2 * W3)
println(simplify(W))

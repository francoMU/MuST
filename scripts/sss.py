from sympy import symbols, Symbol, Function, jn, yn, simplify, Derivative

r, kappa, E = symbols('r, kappa, E')

l = symbols('l', integer=True)

Q = Function('Q')(r)
P = Function('P')(r)
V = Function('V')(r)

C = r * yn(l, r) * Q - r * yn(l, r).diff(r) * P

Q_exp = (Derivative(P, r) - P / r)

d2Pdr2 = 2 * (V - E + (l * (l + 1)) / (2 * r ** 2)) * P

dCdr = simplify(C.diff(r))


print(simplify(dCdr.subs(Q, Q_exp)))

dCdr = simplify(dCdr.subs(Q, Q_exp))

print("""

""")

print(simplify(dCdr.subs(Derivative(P, (r, 2)), d2Pdr2)))

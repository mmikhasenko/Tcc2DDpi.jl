# kinematic functions
s23(v) = v.s + v.msq[1] + v.msq[2] + v.msq[3] - v.s12 - v.s13

# scalar products
p_p12(v)  = (v.s + v.s12 - v.msq[3]) / 2
p12_p2(v) = (v.s12 + v.msq[2] - v.msq[1]) / 2
p_p2(v)   = (v.s + v.msq[2] - v.s13) / 2
# 
p_p13(v)  = (v.s + v.s13 - v.msq[2]) / 2
p13_p3(v) = (v.s13 + v.msq[3] - v.msq[1]) / 2
p_p3(v)   = (v.s + v.msq[3] - v.s12) / 2
#
p1_p2(v) = (v.s12 - v.msq[1] - v.msq[2]) / 2
p1_p3(v) = (v.s13 - v.msq[1] - v.msq[3]) / 2
p2_p3(v) = (s23(v) - v.msq[2] - v.msq[3]) / 2
#
p12_p13(v) = (2 * (v.s + v.msq[1]) - s23(v) - v.s12 - v.s13) / 2
# 	
p12_p3(v) = (v.s - v.s12 - v.msq[3]) / 2
p13_p2(v) = (v.s - v.s13 - v.msq[2]) / 2


# covariant expressions
λ(x, y, z) = x^2 + y^2 + z^2 - 2x * y - 2y * z - 2z * x
A(v) = λ(v.s12, v.msq[1], v.msq[2]) / (4 * v.s12) +
       (p_p12(v) * p12_p2(v) / v.s12 - p_p2(v))^2 / v.s
B(v) = λ(v.s13, v.msq[1], v.msq[3]) / (4 * v.s13) +
       (p_p13(v) * p13_p3(v) / v.s13 - p_p3(v))^2 / v.s
C(v) = D(v) + E(v)
D(v) =
    p12_p3(v) * p12_p2(v) / v.s12 +
    p13_p3(v) * p13_p2(v) / v.s13 -
    p12_p13(v) * p12_p2(v) * p13_p3(v) / (v.s12 * v.s13) -
    p2_p3(v);
E(v) =
    p_p12(v) * p_p13(v) *
    p12_p2(v) * p13_p3(v) / (v.s * v.s12 * v.s13) +
    p_p2(v) * p_p3(v) / v.s -
    p_p12(v) * p_p3(v) * p12_p2(v) / (v.s12 * v.s) -
    p_p13(v) * p_p2(v) * p13_p3(v) / (v.s13 * v.s)
G(v) = (
    2 * p1_p2(v) * p1_p3(v) * p2_p3(v) -
    v.msq[2] * p1_p3(v)^2 -
    v.msq[3] * p1_p2(v)^2) / (2 * v.s)
#

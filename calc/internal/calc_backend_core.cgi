#!/usr/bin/env sage

import sys
import json
import sage.all as sage
import warnings
from math import ceil

# Turn off deprecation warning
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Set the content type to JSON
print("Content-Type: application/json")
print()

# Parse the form data
request_body = sys.stdin.read()
request_data = json.loads(request_body)

# Extract parameters from the JSON data
x = sage.var('x')
function = sage.SR(request_data['function'])
a = float(request_data['a'])
b = float(request_data['b'])
n = int(request_data['n'])
max_error = float(request_data['max_error'])
rule = request_data['rule']

def integrate(f, a, b, n, max_error, rule):
    def midpoint_rule(fcn, a, b, n):
        Deltax = (b - a) / n
        xs = [a + Deltax * i for i in range(n + 1)]
        ysmid = [fcn((xs[i] + xs[i + 1]) / 2) for i in range(n)]
        return Deltax * sum(ysmid)

    def trapezoid_rule(fcn, a, b, n):
        Deltax = (b - a) / n
        coeffs = [2] * (n - 1)
        coeffs = [1] + coeffs + [1]
        valsf = [fcn(x=(a + Deltax * i)) for i in range(n + 1)]
        return (Deltax / 2) * sum([coeffs[i] * valsf[i] for i in range(n + 1)])

    def simpsons_rule(fcn, a, b, n):
        Deltax = (b - a) / n
        n2 = int(n / 2)
        coeffs = [4, 2] * n2
        coeffs = [1] + coeffs[:n - 1] + [1]
        valsf = [fcn(x=(a + Deltax * i)) for i in range(n + 1)]
        return (Deltax / 3) * sum([coeffs[i] * valsf[i] for i in range(n + 1)])

    def diff_upper_bound(fcn, a, b, rule):
        diff2 = abs(sage.diff(sage.diff(fcn, x), x))
        critical_points = sage.solve(sage.diff(diff2, x) == 0, x)
        critical_points = [p for p in critical_points if a < p < b]
        if not critical_points:
            return max(fcn(a), fcn(b))
        return max([fcn(p) for p in critical_points])

    def error(f, a, b, n, rule):
        if rule == "s":
            return diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 5 / (180 * n ** 4)
        elif rule == "m":
            return diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 3 / (24 * n ** 2)
        elif rule == "t":
            return diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 3 / (12 * n ** 2)

    def n_c(f, a, b, max_error, rule):
        if rule == "s":
            return ceil(diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 5 / (180 * max_error) ** (1 / 4))
        elif rule == "m":
            return ceil(diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 5 / (24 * max_error) ** (1 / 2))
        elif rule == "t":
            return ceil(diff_upper_bound(fcn=f, a=a, b=b, rule=rule) * (b - a) ** 5 / (12 * max_error) ** (1 / 2))

    def over_under(f, a, b, rule):
        if rule == "s":
            return "error"
        f_2prime = sage.diff(sage.diff(f, x), x)
        if all(f_2prime.subs(x, xi) > 0 for xi in [xi for xi in [a, b] if xi != a and xi != b]):
            if rule == "m":
                return "underestimate"
            elif rule == "t":
                return "overestimate"
        elif all(f_2prime.subs(x, xi) < 0 for xi in [xi for xi in [a, b] if xi != a and xi != b]):
            if rule == "m":
                return "overestimate"
            elif rule == "t":
                return "underestimate"
        else:
            return "error"

    if n == 0:
        if rule == "s":
            answer = simpsons_rule(f, a, b, n_c(f, a, b, max_error, rule)).n()
        elif rule == "m":
            answer = midpoint_rule(f, a, b, n_c(f, a, b, max_error, rule)).n()
        elif rule == "t":
            answer = trapezoid_rule(f, a, b, n_c(f, a, b, max_error, rule)).n()
        return {
            "answer": float(answer),
            "n": float(n_c(f, a, b, max_error, rule)),
            "error": float(max_error),
            "over_under": over_under(f, a, b, rule)
        }
    else:
        if rule == "s":
            answer = simpsons_rule(f, a, b, n).n()
        elif rule == "m":
            answer = midpoint_rule(f, a, b, n).n()
        elif rule == "t":
            answer = trapezoid_rule(f, a, b, n).n()
        return {
            "answer": float(answer),
            "n": float(n),
            "error": float(error(f, a, b, n, rule)),
            "over_under": over_under(f, a, b, rule)
        }

# Perform the calculation
result = integrate(function, a, b, n, max_error, rule)

# Print the result as JSON
print(json.dumps(result))

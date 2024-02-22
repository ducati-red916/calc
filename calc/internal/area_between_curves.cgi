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
function_1 = sage.SR(request_data['function_1'])
function_2 = sage.SR(request_data['function_2'])
a = float(request_data['a'])
b = float(request_data['b'])



def area(f, g, a, b):
    x = sage.var('x')  # Define the variable symbolically
    
    # Solving the equation f(x) = g(x) to find intersections
    intersections = sage.solve(f == g, x)
    decimal_intersections = [intersection.rhs() for intersection in intersections]

    if a == 0 and b == 0:
        a = decimal_intersections[0]
        b = decimal_intersections[-1]
    decimal_intersections = [p for p in decimal_intersections if a < p < b]

    answer = 0
    if decimal_intersections:
        for decimal_intersection in decimal_intersections:
            if f(decimal_intersection - 0.0001) > g(decimal_intersection - 0.0001):
                answer += sage.integral((f - g).simplify_full(), x, a, decimal_intersection)
            else:
                answer += sage.integral((g - f).simplify_full(), x, a, decimal_intersection)

            if f(decimal_intersection + 0.0001) < g(decimal_intersection + 0.0001):
                answer += sage.integral((f - g).simplify_full(), x, decimal_intersection, b)
            else:
                answer += sage.integral((g - f).simplify_full(), x, decimal_intersection, b)
    else:
        if f(a) == g(a):
            if f(a + 0.0001) > g(a + 0.0001):
                answer = sage.integral((f - g).simplify_full(), x, a, b)
            else:
                answer = sage.integral((g - f).simplify_full(), x, a, b)
        else:
            if f(a) > g(a):
                answer = sage.integral((f - g).simplify_full(), x, a, b)
            else:
                answer = sage.integral((g - f).simplify_full(), x, a, b)

    return answer




# Perform the calculation
result_expression = area(function_1, function_2, a, b)

# Convert the result to a string
result_str = str(result_expression)

# Print the result as JSON
print(json.dumps({"area": result_str}))

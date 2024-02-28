#!/usr/bin/env sage

import sys
import json
import sage.all as sage
import warnings
import numpy as np
import sympy as sp
import os
import matplotlib.pyplot as plt

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
	filepath = "/calc/internal/temp/graph.png"
	if os.path.exists(filepath):
		os.remove(filepath)
	x = sp.symbols('x')
	
	# Convert to SymPy functions
	f_sym = sp.sympify(f)
	g_sym = sp.sympify(g)
	
	# Find intersections
	intersections = sp.solve(f_sym - g_sym, x)
	decimal_intersections = [float(intersection.evalf()) for intersection in intersections]
	
	if a == 0 and b == 0:
		a = decimal_intersections[0]
		b = decimal_intersections[-1]
	decimal_intersections = [p for p in decimal_intersections if a < p < b]
	
	answer = 0
	if decimal_intersections:
		for decimal_intersection in decimal_intersections:
			if f_sym.subs(x, decimal_intersection - 0.0001) > g_sym.subs(x, decimal_intersection - 0.0001):
				answer += sp.integrate((f_sym - g_sym), (x, a, decimal_intersection))
			else:
				answer += sp.integrate((g_sym - f_sym), (x, a, decimal_intersection))
			
			if f_sym.subs(x, decimal_intersection + 0.0001) < g_sym.subs(x, decimal_intersection + 0.0001):
				answer += sp.integrate((f_sym - g_sym), (x, decimal_intersection, b))
			else:
				answer += sp.integrate((g_sym - f_sym), (x, decimal_intersection, b))
	else:
		if f_sym.subs(x, a) == g_sym.subs(x, a):
			if f_sym.subs(x, a + 0.0001) > g_sym.subs(x, a + 0.0001):
				answer = sp.integrate((f_sym - g_sym), (x, a, b))
			else:
				answer = sp.integrate((g_sym - f_sym), (x, a, b))
		else:
			if f_sym.subs(x, a) > g_sym.subs(x, a):
				answer = sp.integrate((f_sym - g_sym), (x, a, b))
			else:
				answer = sp.integrate((g_sym - f_sym), (x, a, b))
	
	# Generate x values
	x_values = np.linspace(min(a, b) - 1, max(a, b) + 1, 1000)
	# Generate y values for each function
	f_lambda = sp.lambdify('x', f_sym)
	g_lambda = sp.lambdify('x', g_sym)
	f_values = f_lambda(x_values)
	g_values = g_lambda(x_values)
	
	# Plot the functions
	plt.plot(x_values, f_values, label='f(x)')
	plt.plot(x_values, g_values, label='g(x)')
	
	# Plot intersections
	for intersection in intersections:
		plt.plot(float(intersection), float(f_sym.subs('x', intersection)), 'ro')  # Plotting intersection points
	
	# Shade the area between the functions using the updated bounds
	x_shade = np.linspace(a, b, 1000)  # Using updated bounds
	f_shade = f_lambda(x_shade)
	g_shade = g_lambda(x_shade)
	plt.fill_between(x_shade, f_shade, g_shade, where=(f_shade > g_shade), color='lightgrey', alpha=0.5)
	plt.fill_between(x_shade, f_shade, g_shade, where=(f_shade < g_shade), color='lightblue', alpha=0.5)
	
	# Add labels and legend
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Visualization of f(x) and g(x)')
	plt.legend()
	
	# Show plot
	plt.grid(True)
	plt.savefig(filepath)
	
	return answer, filepath



area_result, filepath = area(function_1, function_2, a, b)

# Calculate the area between the curves
result = {
	"area": str(area_result),  # Convert to float
	"filepath": filepath
}

# Print the result as JSON
print(json.dumps(result))

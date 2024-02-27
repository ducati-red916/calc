#!/usr/bin/env sage

import sys
import json
import sage.all as sage
import warnings
import math
import numpy as np
import matplotlib.pyplot as plt
import os

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
	filepath="/calc/internal/temp/graph.png"
	if os.path.exists(filepath):
		os.remove(filepath)
	def midpoint_rule(fcn,a,b,n):
		Deltax = (b-a)*1.0/n
		xs=[a+Deltax*i for i in range(n+1)]
		ysmid=[fcn((xs[i]+xs[i+1])/2) for i in range(n)]

		return Deltax*sum(ysmid)
	
	def trapezoid_rule(fcn,a,b,n):
		Deltax = (b-a)*1.0/n
		coeffs = [2]*(n-1)
		coeffs = [1]+coeffs+[1]
		valsf = [fcn(x=(a+Deltax*i)) for i in range(n+1)]

		return (Deltax/2)*sum([coeffs[i]*valsf[i] for i in range(n+1)])
	
	def trapezoid_rule_plot(f, a, b, n, filename):
		# Define the x values for plotting
		x_values = np.linspace(a, b, 1000)
	
		# Evaluate the function at x_values
		y_values = [f(x) for x in x_values]
	
		# Calculate the partition points
		partition_points = np.linspace(a, b, n+1)
	
		# Calculate the partition width
		partition_width = (b - a) / n

		# Plotting
		plt.plot(x_values, y_values, label=f'Function: $y = {f}$')
		plt.xlabel('x')
		plt.ylabel('y')

		# Plot trapezoids, endpoints, and outlines
		for i in range(n):
			x0 = partition_points[i]
			x1 = partition_points[i+1]
			y0 = f(x0)
			y1 = f(x1)
			plt.fill([x0, x1, x1, x0], [0, 0, y1, y0], color='orange', alpha=0.3)
			plt.plot([x0, x1], [y0, y1], 'ko')	  # Black dots at endpoints
			plt.plot([x0, x0, x1, x1, x0], [0, y0, y1, 0, 0], 'r-', linewidth=0.5)	  # Red outline

		plt.title('Function Plot with Trapezoids, Endpoints, and Outlines')
		plt.grid(True)
		plt.legend()
		plt.savefig(filename)
		
	def midpoint_rule_plot(f, a, b, n, filename):
		# Define the x values for plotting
		x_values = np.linspace(a, b, 1000)
	
		# Evaluate the function at x_values
		y_values = [f(x) for x in x_values]
	
		# Calculate the partition points
		partition_points = np.linspace(a, b, n+1)
	
		# Calculate the partition width
		partition_width = (b - a) / n

		# Plotting
		plt.plot(x_values, y_values, label=f'Function: $y = {f}$')
		plt.xlabel('x')
		plt.ylabel('y')

		# Plot rectangles, midpoints, and outlines
		for i in range(n):
			x0 = partition_points[i]
			x1 = partition_points[i+1]
			midpoint = (x0 + x1) / 2  # Calculate midpoint
			y_mid = f(midpoint)	 # Evaluate function at midpoint
			plt.fill([midpoint - partition_width/2, midpoint + partition_width/2, midpoint + partition_width/2, midpoint - partition_width/2], 
					[0, 0, y_mid, y_mid], color='orange', alpha=0.3)  # Rectangle centered at midpoint
			plt.plot(midpoint, y_mid, 'ko')	 # Black dot at midpoint
			plt.plot([x0, x1], [y_mid, y_mid], 'ko')  # Black dots at endpoints
			plt.plot([midpoint - partition_width/2, midpoint - partition_width/2, midpoint + partition_width/2, midpoint + partition_width/2, midpoint - partition_width/2], 
					[0, y_mid, y_mid, 0, 0], 'r-', linewidth=0.5)  # Red outline
	
		plt.title('Function Plot with Rectangles (Midpoint Rule), Midpoints, and Outlines')
		plt.grid(True)
		plt.legend()
		plt.savefig(filename)

	
	def simpsons_rule(fcn,a,b,n):
		if n%2 != 0:
			n+=1
		Deltax = (b-a)*1.0/n
		n2=int(n/2)
		coeffs = [4,2]*n2
		coeffs = [1] +coeffs[:n-1]+[1]
		valsf = [fcn(x=(a+Deltax*i)) for i in range(n+1)]

		return (Deltax/3)*sum([coeffs[i]*valsf[i] for i in range(n+1)])

	def diff_upper_bound(fcn, a, b, rule):
		#x=var('x')
		if rule != "s":
			diff2 = abs(sage.diff(sage.diff(fcn, x), x))
		else:
			diff2 = abs(sage.diff(sage.diff(sage.diff(sage.diff(fcn, x), x), x), x))
		critical_points = sage.solve(sage.diff(diff2, x) == 0, x)
		critical_points = [p for p in critical_points if a < p < b]
		if not critical_points:
			return max(fcn(a), fcn(b))
		return max([fcn(p) for p in critical_points])

	def error(f, a, b, n, rule):
		if rule == "s":
			return(diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**5/(180*n**4))
		elif rule == "m":
			return(diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**3/(24*n**2))
		elif rule == "t":
			return(diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**3/(12*n**2))
		
	def n_c(f, a, b, max_error, rule):
		if rule == "s":
			return(math.ceil((diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**5/(180*max_error)**(1/4))))
		elif rule == "m":
			return(math.ceil((diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**5/(24*max_error)**(1/2))))
		elif rule == "t":
			return(math.ceil((diff_upper_bound(fcn=f, a=a, b=b, rule=rule)*(b-a)**5/(12*max_error)**(1/2))))
		
	
	def over_under(f, a, b, rule):
		if rule == "s":
			return("error")
		f_2prime = sage.diff(sage.diff(f, x), x)
		if all(f_2prime.subs(x, xi) > 0 for xi in [xi for xi in [a, b] if xi != a and xi != b]):
			if rule == "m":
				return("underestimate")
			elif rule == "t":
				return("overestimate")
		elif all(f_2prime.subs(x, xi) < 0 for xi in [xi for xi in [a, b] if xi != a and xi != b]):
			if rule == "m":
				return("overestimate")
			elif rule == "t":
				return("underestimate")
		else:
			return("error")

	if (n==0):
		if rule == "s":
			answer = simpsons_rule(f,a,b,n_c(f,a,b,max_error,rule)).n()
		elif rule == "m":
			answer = midpoint_rule(f,a,b,n_c(f,a,b,max_error,rule)).n()
			midpoint_rule_plot(f,a,b,n_c(f,a,b,max_error,rule),filepath)
		elif rule == "t":
			answer = trapezoid_rule(f,a,b,n_c(f,a,b,max_error,rule)).n()
			trapezoid_rule_plot(f,a,b,n_c(f,a,b,max_error,rule),filepath)
		return {
			"answer": float(answer),
			"n": float(n_c(f,a,b,max_error,rule)),
			"error": float(max_error),
			"over_under": over_under(f, a, b, rule),
			"filepath": filepath
		}

	else:
		if rule == "s":
			answer = simpsons_rule(f,a,b,n).n()
		elif rule == "m":
			answer = midpoint_rule(f,a,b,n).n()
			midpoint_rule_plot(f,a,b,n,filepath)
		elif rule == "t":
			answer = trapezoid_rule(f,a,b,n).n()
			trapezoid_rule_plot(f,a,b,n,filepath)
		return {
			"answer": float(answer),
			"n": float(n),
			"error": float(error(f, a, b, n, rule)),
			"over_under": over_under(f, a, b, rule),
			"filepath": filepath
		}

# Perform the calculation
result = integrate(function, a, b, n, max_error, rule)

# Print the result as JSON
print(json.dumps(result))

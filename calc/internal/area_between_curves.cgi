#!/usr/bin/env sage

import sys
import json
import sage.all as sage
import warnings
import numpy as np
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

# Define the x range for plotting
x_values = np.linspace(a, b, 400)

# Evaluate the functions at x_values
y_values_1 = [function_1(x_val).n() for x_val in x_values]
y_values_2 = [function_2(x_val).n() for x_val in x_values]

# Plot the functions
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values_1, label='Function 1: $y = {}$'.format(function_1))
plt.plot(x_values, y_values_2, label='Function 2: $y = {}$'.format(function_2))

# Shade the area between the curves
plt.fill_between(x_values, y_values_1, y_values_2, where=(y_values_1 > y_values_2), color='orange', alpha=0.3)
plt.fill_between(x_values, y_values_2, y_values_1, where=(y_values_2 > y_values_1), color='orange', alpha=0.3)

# Highlight the x-axis
plt.axhline(0, color='black', linewidth=0.5)

# Set labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Area Between Curves')

# Add legend
plt.legend()

# Save the plot as an image
graph_path = "/calc/internal/temp/area_between_curves_graph.png"
plt.savefig(graph_path)

# Close the plot to release memory
plt.close()

# Calculate the area between the curves
def area(f1, f2, a, b):
    return sage.integral(abs(f1 - f2), x, a, b).n()

# Perform the calculation
result = {
    "area": area(function_1, function_2, a, b),
    "graph_path": graph_path
}

# Print the result as JSON
print(json.dumps(result))

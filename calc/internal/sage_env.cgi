#!/bin/sh

# Set the DOT_SAGE environment variable
DOT_SAGE=/tmp/sage_apache
export DOT_SAGE

# Parse the query string to get the value of the "action" parameter
query_string="$QUERY_STRING"
action_value=$(echo "$query_string" | awk -F'=' '{print $2}')
sage "$action_value"
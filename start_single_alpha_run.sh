#!/bin/sh

# Inside every single folder the following sequence of commands should be used
  
input_a='alpha_lifetime.inp'
  
input_a_offset='alpha_lifetime_offset.inp'

input_a_total='alpha_lifetime_total.inp'
 

# Compute offset time
cp $input_a_offset $input_a

/usr/bin/time -o "cpu_time_offset.dat" ./gorilla_applets_main.x

# Compute total time
cp $input_a_total $input_a

/usr/bin/time -o "cpu_time_total.dat" ./gorilla_applets_main.x


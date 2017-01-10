#!/bin/bash
# This example bash script:
# - calls the X0h program from http://x-server.gmca.aps.anl.gov/
# - gets the scattering factors chi_0r and chi_0i
# - saves them into a file as a function of energy
#
# The script needs the following external programs:
# 1) lynx : text-based web browser
# 2) sed  : GNU stream editor
# 3) bc   : arbitrary-precision calculator language
#
# On Microsoft Windows this script can run under
# Cygwin (http://www.cygwin.com).
#
# Original version: Mojmir Meduna  27.08.2003
# Current version: Sergey Stepanov 31.12.2005
# Non-integer math: Jackson Williams 11.06.2014
# Altered by Eric Landahl 12.07.2016 for his own purposes
#=========================================================

#--------------Input parameters---------------------------
### Output file:
file=example.dat

### Energy range and the number of energy points
E1=8.8976		#start energy
E2=12.398		#end energy
n=2		#number of pts (please stay within a few dozen!)

### Energy step:
dE=$(echo "($E2-$E1)/($n-1)" | bc -l)

###--------------Parameters of X0h CGI interface--------------
xway=2           # 1 - wavelength, 2 - energy, 3 - line type
#wave=$1         # works with xway=2 or xway=3
# line=Cu-Ka1    # works with xway=3 only
line=            # works with xway=3 only

###Target:
coway=0          # 0 - crystal, 1 - other material, 2 - chemicalformula
###Crystal
code=Germanium   # works with coway=0 only
###Other material:
amor=            # works with coway=1 only
###Chemical formula:
chem=            # works with coway=2 only
###and density (g/cm3):
rho=             # works with coway=2 only

### Miller indices:
i1=4
i2=0
i3=0

###Database Options for dispersion corrections df1, df2:
### -1 - Automatically choose DB for f',f"
###  0 - Use X0h data (5-25 keV or 0.5-2.5 A) -- recommended for Bragg diffraction.
###  2 - Use Henke data (0.01-30 keV or 0.4-1200 A) -- recommended for soft x-rays.
###  4 - Use Brennan-Cowan data (0.03-700 keV or 0.02-400 A)
### 10 - Compare results for all of the above sources.
df1df2=-1

modeout=1	# 0 - html out, 1 - quasy-text out with keywords
detail=0	# 0 - don't print coords, 1 = print coords
###-----------------------------------------------------------



{
###--------Print header------------------------
echo "#Energy, xr0, xi0, xrh, xih, delta, QB"

###--------Loop over energy--------------------
for ((i=0; i < $n; i++))	# loop over energy points
do
  wave=$(echo "$E1+$dE*$i" | bc -l)

###--------Building address--------------------
  address='http://x-server.gmca.aps.anl.gov/cgi/x0h_form.exe?xway='$xway
  address=$address'&wave='$wave
  address=$address'&line='$line
  address=$address'&coway='$coway
  address=$address'&code='$code
  address=$address'&amor='$amor
  address=$address'&chem='$chem
  address=$address'&rho='$rho
  address=$address'&i1='$i1
  address=$address'&i2='$i2
  address=$address'&i3='$i3
  address=$address'&df1df2='$df1df2
  address=$address'&modeout='$modeout
  address=$address'&detail='$detail

###-----------Connect & Download-----------------
### Find line with keyword, erase everything before the data, and print:
# x=(`lynx -dump $address | sed -n '/xr0=/{s/.*xr0=//;p;};/xi0=/{s/.*xi0=//;p;}'`)
#  x=(`curl --silent $address | sed -n '/xr0=/{s/.*xr0=//;p;};/xi0=/{s/.*xi0=//;p;};/xrhsg=/{s/.*xrhsg=//;p;};/xihsg=/{s/.*xihsg=//;p;};/delta=/{s/.*delta=//;p;};/QB=/{s/.*QB=//;p;}'`)
 x=(`curl --silent $address | sed -n '/xr0=/{s/.*xr0=//;p;};/xi0=/{s/.*xi0=//;p;};/xrhsg=/{s/.*xrhsg=//;p;};/xihsg=/{s/.*xihsg=//;p;}'`)
 y=(`curl --silent $address | sed -n '/delta=/{s/.*delta=//;p;};/QB=/{s/.*QB=//;p;}'`)

###--------Print current point-------------------
  echo "   ${wave}, ${x[0]}, ${x[1]}, ${x[2]}, ${x[3]}, ${y[0]}, ${y[1]}"

done
} > $file

# Get rid of annoying linefeed characters
tr -d "\015" <example.dat >example2.dat
rm example.dat
mv example2.dat example.dat


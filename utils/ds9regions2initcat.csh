#!/bin/tcsh
#=======================================================================
#+
# NAME:
#   ds9reg2initcat.csh
#
# PURPOSE:
#   Read in a ds9 region file and output a LT-format Nebula catalog, for 
#   use in initializing LT. 
#
# COMMENTS:
#   Positions must be in WCS degrees. 
#   First source is taken to be the galaxy, all others are the pt srcs.
#
# INPUTS:
#   regfile           ds9 regions file
#  
# OPTIONAL INPUTS:
#   -h --help 
#   -f --filters K    No. of filters in dataset [1] 
#
# OUTPUTS:
#   catalog           regfile:r.cat           
#
# EXAMPLES:
#
#
# BUGS:
#
# REVISION HISTORY:
#   2014-05-06  started Marshall and Agnello (UCSB)
#-
#=======================================================================

# Options and arguments:

set help = 0
set regfiles = ()
set Nbands = 1

while ( $#argv > 0 )
   switch ($argv[1])
   case -h:
      set help = 1
      shift argv
      breaksw
   case --help:
      set help = 1
      shift argv
      breaksw
   case -f:
      shift argv
      set Nbands = $argv[1]
      shift argv
      breaksw
   case --filters:
      shift argv
      set Nbands = $argv[1]
      shift argv
      breaksw
   case *:
      set regfiles = ( $regfiles $argv[1] )
      shift argv
      breaksw
   endsw
end

#-----------------------------------------------------------------------

# Catch stupidities, set up variables:

if ($help || $#regfiles == 0) then
  more `which $0`
  goto FINISH
endif

echo "Processing $#regfiles ds9 region files..."

set SED = ""
foreach k ( `seq $Nbands` )
    set SED = "$SED 99"
end    

#-----------------------------------------------------------------------

# Loop over region files:

foreach regfile ( $regfiles )

    # Read in ra and dec values, in degrees. 
    # NB first row is taken to be the galaxy, all the rest are point sources.

    set x = `grep point=cross $regfile | \
             cut -d'(' -f2 | cut -d')' -f1 | sed s/','/' '/g`

    # Start catalog:
    set date = `date`
    set catalog = $regfile:r.cat
    echo "# Dummy LensTractor-format catalog containing Nebula positions" > $catalog
    echo "# Data source: $regfile" >> $catalog
    echo "# Written by $0:t on $date" >> $catalog
    if ($#x == 10) then
        set model = 'Nebula4'
    else if ($#x == 6) then
        set model = 'Nebula2'
    else 
        echo "ERROR: confused by region file, found $#x position values but was expecting 6 or 10..."
        goto FINISH
    endif 
    echo "# Model: $model" >> $catalog
    echo "# " >> $catalog

    set line = `echo "$x[1]  $x[2]  $SED 99 99 99  $x[3]  $x[4] $SED  $x[5] $x[6] $SED"`
    if ($model == 'Nebula4') then 
        set line = `echo "$line  $x[7]  $x[8] $SED  $x[9]  $x[10] $SED"`
    endif
    
    echo "$line" >> $catalog
    wc -l $catalog

end

#-----------------------------------------------------------------------

FINISH:

#=======================================================================

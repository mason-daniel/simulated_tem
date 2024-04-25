#!/bin/bash
#   short script to test output of d2bi code.
#   returns "PASS" if the output of the command
#        ./bin/d2bi -f ../Data/testhcpMonovac.xyz -a0 3 -lattice hcp -k 0,0,0,1 -g 1,0,-1,0  -nousePhase -xifile ../Polycrystal-Analysis/Data/extinctionDistances.Zr_V200_T300.dat 
#   gives a single vacancy TEM image
#   usage:
#       testD2BI.sh [data_directory]
  


if [ "$#" -ne 1 ]; then
    echo "usage: testD2BI.sh [data_directory]"
    echo "FAIL"
else
    #   there should be a line in the output reading
    #   d2bi info - minmaxavg d2bimg         0.000000    0.000166    0.000050 (diffracted beam)
    #   so multiply columns x 1000000 to give 0,132,15                                                                                           

    data_directory=$1
    test_xyz_file=test.i111_loop_W.lammps
    target_min=0.000028         #   target output  
    target_max=0.230127
    target_avg=0.081501    
    accuracy=0.0039             #    = 1/256, the limit of greyscale intensity resolution 
    xi_file=${data_directory}/extinctionDistances.W_V200_T300.dat

    mpirun -n 1 ../bin/d2bi  -f ${data_directory}/${test_xyz_file} -a0 3.1652 -lattice bcc -k 0,0,1 -g 2,0,0 -xifile ${xi_file} -ng 1 -o "" > d2biresult 

    echo "target result"
    echo "minmaxavg d2bimg " ${target_min} " " ${target_max} " " ${target_avg} 

    echo "result from d2bi"
    echo `grep "minmaxavg d2bimg" d2biresult `


    fdelta=`grep "minmaxavg d2bimg" d2biresult |    \
        awk -v target_min="${target_min}" -v target_max="${target_max}" -v target_avg="${target_avg}" -v accuracy="${accuracy}" \
        'function max(a,b){return a>b?a:b} ; function abs(v) {return v < 0 ? -v : v} ;  \
        {print int( max(  abs( ($6 - target_min)/accuracy + 0.5 )  , max(                   \
                          abs( ($7 - target_max)/accuracy + 0.5 )  ,                        \
                          abs( ($8 - target_avg)/accuracy + 0.5 ) ) ) )}' `

    echo "fdelta = " ${fdelta} 
     
    if [ "${fdelta}" -eq "0" ]; then
        echo -e "\033[1;32m PASS"
        rm d2biresult
    else
        cat d2biresult
        echo -e "\033[0;31m FAIL"
    fi
    
fi
  
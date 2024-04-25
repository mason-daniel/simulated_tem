#!/bin/bash
#   short script to test output of findVoids code.
#   returns "PASS" if the output of the command
#        ./bin/findVoids.exe -f ../Data/Simple_W_cases/monovac_W_8.xyz -a0 3.1652
#   gives a single vacancy
#   usage:
#       testVoids.sh [data_directory]
 

if [ "$#" -ne 1 ]; then
    echo "usage: testVoids.sh [data_directory]"
    echo "FAIL"
else
    #   there should be a line in the output reading
    #   area,volume,volume/omega0,nVoids,nVacs       36.01321      15.80517       0.99684             1             1
    #   so extract column 4, multiply by 10, should give = 10 
    ../bin/findVoids -f $1/monovac_W_8.xyz -a0 3.1652 > voidresult
    nvactimes10=`grep "nVoids" voidresult | awk '{print int($4*10+0.5)}'` 
    if [ "${nvactimes10}" == "10" ]; then
        echo -e "\033[1;32m PASS"
        rm voidresult
    else
        echo -e "\033[0;31m FAIL"
    fi
    
fi

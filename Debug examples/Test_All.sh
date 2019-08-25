#!/bin/bash

# Define the path to the working directory. This path is system dependent and needs to be changed on different systems
WorkDirectory='/Users/haverkort/ProgramOutput/Test Quanty'

# Test if the working directory excists, if not create it
if [ ! -e "$WorkDirectory" ]
then
  echo "Working directory did not excist:"
  echo "  Made directory: $WorkDirectory"
  mkdir "$WorkDirectory"
fi

# Within the working directory create a unique directory for the current calculation
WorkDirectory="$WorkDirectory/"`date +%Y.%m.%d.%H.%M`" _running"

if [ ! -e "$WorkDirectory" ]
then
  echo "Running test set in directory:"
  echo "$WorkDirectory"
  mkdir "$WorkDirectory"
else
  echo "Working directory already excists"
  echo "stoped test"
  exit 1
fi

#Define a logfile
LogFile="$WorkDirectory/logfile"

# All tests should pass, if not we count the number of failures
NumberOfFailures=0

# move all examples to the working directory
echo "Time " `date +"%Y %m %d %H:%M:%S"` | tee "$LogFile"
cp -r Reference/* "$WorkDirectory"
echo "Copied all files to the working directory" | tee -a "$LogFile"

# move to working directory
echo "Time " `date +"%Y %m %d %H:%M:%S"` | tee -a "$LogFile"
cd "$WorkDirectory"
echo "Changed to working directory $WorkDirecotry" | tee -a "$LogFile"

for directory in */
do
  # move to subdirectory in workdirectory
  echo "Time " `date +"%Y %m %d %H:%M:%S"` | tee -a "$LogFile"
  cd "$directory"
    echo "Changing to directory $directory" | tee -a "$LogFile"
    if [ "$directory" == "Objects/" ]
    then
      for directory2 in */
      do
        # move to subdirectory in workdirectory
        echo "Time " `date +"%Y %m %d %H:%M:%S"` | tee -a "$LogFile"
        cd "$directory2"
          echo "Changing to directory $directory2" | tee -a "$LogFile"
          for file in *.Quanty
          do
            line="--------------------------------------------------"
            line="$file ${line:${#file}}"
            echo "  Running $line " `date +"%d %H:%M:%S"` | tee -a "$LogFile"
            Quanty "$file" > "${file%???????}.out"
            exitval=$?
            if [ "$exitval" != "0" ]
            then
              echo "FAILURE with exit code $exitval" | tee -a "$LogFile"
              NumberOfFailures=$((NumberOfFailures + exitval))
            fi
          done
        cd ..
      done
    else
      for file in *.Quanty
      do
        line="--------------------------------------------------"
        line="$file ${line:${#file}}"
        echo "  Running $line " `date +"%d %H:%M:%S"` | tee -a "$LogFile"
        Quanty "$file" > "${file%???????}.out"
        exitval=$?
        if [ "$exitval" != "0" ]
        then
          echo "FAILURE with exit code $exitval" | tee -a "$LogFile"
          NumberOfFailures=$((NumberOfFailures + exitval))
        fi
      done
    fi
  cd ..
done










echo "  Time " `date +"%Y %m %d %H:%M:%S"` | tee -a "$LogFile"
# If the test pased change the name of the directory from _running to _finished_successfully
# If the test failed chang the name of the direcotry from _runnign to _fialed with the number of failures added to it
if [ $NumberOfFailures -eq 0 ]
then
  echo "Test succeded moved directory to"
  echo "${WorkDirectory%????????}_finished_successfully"
  mv "$WorkDirectory" "${WorkDirectory%????????}_finished_successfully"
else
  echo "Test failed moved directory to"
  echo "${WorkDirectory%????????}_failed $NumberOfFailures"
  mv "$WorkDirectory" "${WorkDirectory%????????}_failed $NumberOfFailures"
fi

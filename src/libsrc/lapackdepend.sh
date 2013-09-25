#!/bin/bash

name=$1

src=./SRC
srcfiles=" $(find $src -name $name)"
n=1

# --------------------------------------------------------------------# 

for ((;;))
do

  newfiles=""

  for f in "$srcfiles"
  do
    procedures=$( \
      cat $f \
      | sed '/^[^ ]/d' \
      | sed 'N;N;N;N;N; s/\n *\$//g' \
      | sed 'N;N;N;N;N; s/\n *\$//g' \
      | grep EXTERNAL \
      | sed 's/^ *EXTERNAL */ /g;' \
      | sed "s/, */\\
 /g" )
  
    procedures=$(echo "$procedures" | sort | uniq)
  
    for proc in $procedures
    do
      newfiles="${newfiles}
 $(find $src -name *.f | grep -i $proc.f)"
    done
#    newfiles=$(echo "$newfiles" | sed '/^ *$/d' | sort | uniq)

  done

  srcfiles=$(echo "$srcfiles
$newfiles" | sed '/^ *$/d' | sort | uniq)
  nfiles=$(echo "$srcfiles" | wc -l)

  if [ $nfiles -eq $n ]
  then 
      break
  fi

  n=$nfiles

done

echo "$srcfiles"

# --------------------------------------------------------------------# 

exit 0
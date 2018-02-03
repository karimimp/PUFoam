#!/bin/bash
./prepTestCase.sh
blockMesh
setFields
rm -fv log
echo "Running . . ."
PUFoam >& log

echo "Done."
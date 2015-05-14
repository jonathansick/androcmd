#!/bin/bash

getCert

cd $HOME/code/m31hst
git pull
python setup.py develop

cd $HOME/code/androcmd
git pull
python setup.py develop

cd $HOME/code/starfisher
git pull
python setup.py develop

cd $STARFISH

NAME="$1_b$2"
phat_baseline_test.py $NAME $2
tar -zcf ${NAME}.tar.gz $NAME
vcp ${NAME}.tar.gz vos:jonathansick/phat

echo DONE

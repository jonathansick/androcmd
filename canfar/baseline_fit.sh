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

cd $HOME/code/padova
git pull
python setup.py develop

cd $STARFISH

NAME="$1_b$2"
OUTPUT_NAME="${NAME}_$3"

if [ ! -f $NAME ]; then
    vcp vos:jonathansick/phat/${NAME}.tar.gz ${NAME}.tar.gz
    tar -zxf ${NAME}.tar.gz $NAME
fi

phat_baseline_test.py $NAME $2 --fit $3
tar -zcf ${OUTPUT_NAME}.tar.gz $NAME/$3
vcp $STARFISH/${OUTPUT_NAME}.tar.gz vos:jonathansick/phat

echo DONE

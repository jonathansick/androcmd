#!/bin/bash

# Arguments
# 1 - brick
# 2 - patch number(s) (in brick); comma delimited
# 3 - vos dir. e.g. phat/patches

export HOME="/home/jonathansick"
export PATH="/home/jonathansick/anaconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games"
export PHATV2DATA="/home/jonathansick/phat_data"
export DRAINEDATA="/home/jonathansick/draine_data"
export STARFISH="/home/jonathansick/code/starfish"

cd $HOME

pwd

echo "Variables"
echo $DRAINEDATA
echo $PHATV2DATA
echo $1
echo $2
echo $3

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

# Get the JSON dataset describing patches
PATCH_INFO_PATH="patches.json"
if [ ! -f $PATCH_INFO_PATH ]; then
    vcp vos:jonathansick/$3/${PATCH_INFO_PATH} ${PATCH_INFO_PATH}
fi

fit_phat_patches.py $1 --patches $2 --json ${PATCH_INFO_PATH} --vodir vos:jonathansick/$3

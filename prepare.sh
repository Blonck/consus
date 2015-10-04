#!/bin/bash

GIT_REPOS=("https://github.com/philsquared/Catch.git"
           "https://github.com/davisking/dlib.git"
           "https://github.com/RLovelett/eigen.git"
           "https://github.com/google/benchmark.git") 
GIT_BRANCHS=("master" "master" "master" "master")
ZIPS=("http://www.agner.org/optimize/vectorclass.zip" )
ZIPS_NAMES=("vectorclass")

if [ "${#GIT_REPOS[@]}" -ne "${#GIT_BRANCHS[@]}" ]; then
  "error in preparation script"
  exit 1
fi

if [ "${#ZIPS[@]}" -ne "${#ZIPS_NAMES[@]}" ]; then
  "error in preparation script"
  exit 1
fi

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
INSTALL_DIR="${DIR}/extern"
mkdir -p $INSTALL_DIR
echo "external libraries installed in $INSTALL_DIR"

cd $INSTALL_DIR

# install gitrepos
for (( i=0; i<${#GIT_REPOS[@]}; i++ ));
do
  gitrepo=${GIT_REPOS[i]}
  GITDIR=${gitrepo##*/}
  GITDIR=${GITDIR%.*}
  echo $GITDIR
  if [ -d $GITDIR ]
  then
    cd $GITDIR
    git checkout master
    git pull
    git checkout ${GIT_BRANCHS[i]}
    cd $INSTALL_DIR
  else
    git clone $gitrepo
    cd $GITDIR
    git checkout ${GIT_BRANCHS[i]}
    cd $INSTALL_DIR
  fi
done

# download and extract zip files
for (( i=0; i<${#ZIPS[@]}; i++ ));
do
  ZIPDIR=$INSTALL_DIR/${ZIPS_NAMES[i]}
  if [ -d $ZIPDIR ]
  then
    echo "$ZIPDIR already exists"
  else
    mkdir $ZIPDIR
    wget -P $ZIPDIR ${ZIPS[i]}
    cd $ZIPDIR
    unzip *.zip
  fi
done

# create symlinks to git hooks
cd $DIR

echo "enable git hooks"

cp -vf git-hooks/post-checkout .git/hooks/post-checkout
chmod u+x .git/hooks/post-checkout
cp -vf git-hooks/post-commit .git/hooks/post-commit
chmod u+x .git/hooks/post-commit
cp -vf git-hooks/post-merge .git/hooks/post-merge
chmod u+x .git/hooks/post-merge
cp -vf git-hooks/post-update .git/hooks/post-update
chmod u+x .git/hooks/post-update

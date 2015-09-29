#!/bin/bash
git diff --exit-code >> /dev/null && git diff --cached --exit-code >> /dev/null
SUCCESS=$?
HASH=$(git rev-parse --verify HEAD)
if [ $SUCCESS -eq 0 ]
then 
  echo "write src/betamc_version.hpp, all changes commited"
  MESSAGE="${HASH}"
else
  echo "write src/betamc_version.hpp, but not all changes are commited"
  MESSAGE="Not all changes commited, last commit: ${HASH}"
fi
echo "#define GITHASH_CONSUS \"$MESSAGE\"" > src/consus_version.hpp

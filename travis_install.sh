#!/bin/bash
set -v

if [ "${TRAVIS_BRANCH}" = "master" ]; then
	pip install git+https://github.com/sbrg/ecolime.git@master
else
    pip install git+https://github.com/sbrg/ecolime.git@devel
fi
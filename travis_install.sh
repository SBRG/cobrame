#!/bin/bash
set -v

pip install git+https://github.com/${TRAVIS_REPO_SLUG%/*}/ecolime.git@${TRAVIS_BRANCH}
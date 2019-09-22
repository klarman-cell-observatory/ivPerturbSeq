#!/bin/bash

script_name=$0
script_full_path=$(dirname "$0")

echo $(dirname "${BASH_SOURCE[0]}")
echo $BASH_SOURCE
echo $script_full_path

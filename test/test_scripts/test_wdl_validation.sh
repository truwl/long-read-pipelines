#!/bin/bash

find wdl/ -name '*.wdl' | xargs -t -L1 womtool validate

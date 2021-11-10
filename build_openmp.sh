#!/bin/bash
cd build && cmake .. -DUSE_OMP=true -DCMAKE_BUILD_TYPE=Debug && make && mv micro_aevol_cpu micro_aevol_cpu_openmp